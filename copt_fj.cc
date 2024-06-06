/*
 * XPress integration for the Feasibility Jump heuristic.
 */

#include "coptcpp_pch.h"

#include <cstdio>
#include <string>
#include <chrono>
#include <vector>
#include <cassert>
#include <mutex>
#include <thread>
#include <functional>
#include <atomic>
#include <cmath>
#include <climits>
#include <chrono>

#include "feasibilityjump.hh"

const int NUM_THREADS = 2;

std::atomic_size_t totalNumSolutionsFound(0);
std::atomic_size_t totalNumSolutionsAdded(0);
std::atomic_bool presolveFinished(false);
std::atomic_bool heuristicFinished(false);
std::chrono::steady_clock::time_point startTime;

struct Solution
{
	std::vector< double > assignment;
	bool includesContinuous;
};

std::vector< Solution > heuristicSolutions;
std::mutex heuristicSolutions_mutex;

std::mutex presolvedProblem_mutex;
std::mutex nonPresolvedProblem_mutex;

struct ProblemInstance
{
	int numCols;
	std::vector< char > varTypes;
	std::vector< double > lb;
	std::vector< double > ub;
	std::vector< double > objCoeffs;

	int numRows;
	int numNonZeros;
	std::vector< char > rowtypes;
	std::vector< double > rhs;
	std::vector< double > rhsrange;
	std::vector< int > rowStart;
	std::vector< int > colIdxs;
	std::vector< double > colCoeffs;
};

struct FJData
{
	std::vector< int > originalIntegerCols;
	Model originalProblemCopy = nullptr;
	ProblemInstance originalData;
	ProblemInstance presolvedData;
};

FJData gFJData;

// A function that receives the result of any solution added with XPRSaddmipsol.

std::string inputFilename;
std::string outDir;

ProblemInstance getCOPTProblemData(Model problem)
{
	ProblemInstance data;

	data.numCols = problem.GetIntAttr(COPT_INTATTR_COLS);
	data.varTypes = std::vector< char >(data.numCols);
	data.lb = std::vector< double >(data.numCols);
	data.ub = std::vector< double >(data.numCols);
	data.objCoeffs = std::vector< double >(data.numCols);
	for (int colIdx = 0; colIdx < data.numCols; colIdx += 1)
	{
		Var v = problem.GetVar(colIdx);
		data.varTypes[colIdx] = v.GetType();
		data.lb[colIdx] = v.Get(COPT_DBLINFO_LB);
		data.ub[colIdx] = v.Get(COPT_DBLINFO_UB);
		data.objCoeffs[colIdx] = v.Get(COPT_DBLINFO_OBJ);
	}

	data.numRows = problem.GetIntAttr(COPT_INTATTR_ROWS);
	data.rowtypes = std::vector< char >(data.numRows);
	data.rhs = std::vector< double >(data.numRows);
	data.rhsrange = std::vector< double >(data.numRows);
	for (int rowIdx = 0; rowIdx < data.numRows; rowIdx += 1)
	{
		Constraint c = problem.GetConstr(rowIdx);
		double lb = c.Get(COPT_DBLINFO_LB);
		double ub = c.Get(COPT_DBLINFO_UB);
		if (lb == ub)
		{
			data.rowtypes[rowIdx] = 'E';
			data.rhs[rowIdx] = lb;
		}
		else if (lb == -COPT_INFINITY)
		{
			data.rowtypes[rowIdx] = 'L';
			data.rhs[rowIdx] = ub;
		}
		else if (ub == COPT_INFINITY)
		{
			data.rowtypes[rowIdx] = 'G';
			data.rhs[rowIdx] = lb;
		}
		else
		{
			throw std::runtime_error("Range constraints not supported.");
			assert(lb < ub);
			data.rowtypes[rowIdx] = 'R';
			data.rhs[rowIdx] = lb;
			data.rhsrange[rowIdx] = ub - lb;  // FIXME?
		}
	}

	data.numNonZeros = problem.GetIntAttr(COPT_INTATTR_ELEMS);

	printf(FJ_LOG_PREFIX "copying %d x %d matrix with %d nonzeros.\n",
		data.numCols, data.numRows, data.numNonZeros);
	data.rowStart = std::vector< int >(data.numRows + 1);
	data.colIdxs = std::vector< int >(data.numNonZeros);
	data.colCoeffs = std::vector< double >(data.numNonZeros);

	// save it as CSR format
	for (int rowIdx = 0; rowIdx < data.numRows; rowIdx++)
	{
		Constraint c = problem.GetConstr(rowIdx);
		Expr r = problem.GetRow(c);
		size_t rowStart = data.rowStart[rowIdx];
		data.rowStart[rowIdx + 1] = rowStart + r.Size();
		for (int i = 0; i < r.Size(); i += 1)
		{
			data.colIdxs[rowStart + i] = r.GetVar(i).GetIdx();
			data.colCoeffs[rowStart + i] = r.GetCoeff(i);
		}
	}

	return data;
}

bool copyDataToHeuristicSolver(FeasibilityJumpSolver& solver, ProblemInstance& data, int relaxContinuous)
{
	printf("initializing FJ with %d vars %d constraints\n", data.numCols, data.numRows);
	for (int colIdx = 0; colIdx < data.numCols; colIdx += 1)
	{
		VarType vartype = VarType::Continuous;
		if (data.varTypes[colIdx] == 'C')
		{
			vartype = VarType::Continuous;
		}
		else if (data.varTypes[colIdx] == 'I')
		{
			vartype = VarType::Integer;
		}
		else if (data.varTypes[colIdx] == 'B')
		{
			vartype = VarType::Integer;
		}
		else
		{
			printf(FJ_LOG_PREFIX "unsupported variable type '%c' (%d).\n",
				data.varTypes[colIdx], data.varTypes[colIdx]);
			return false;
		}

		solver.addVar(vartype, data.lb[colIdx], data.ub[colIdx], data.objCoeffs[colIdx]);
	}

	for (int rowIdx = 0; rowIdx < data.numRows; rowIdx += 1)
	{
		RowType rowtype;
		if (data.rowtypes[rowIdx] == 'N')
		{
			continue;
		}
		else if (data.rowtypes[rowIdx] == 'L')
		{
			rowtype = RowType::Lte;
		}
		else if (data.rowtypes[rowIdx] == 'G')
		{
			rowtype = RowType::Gte;
		}
		else if (data.rowtypes[rowIdx] == 'E')
		{
			rowtype = RowType::Equal;
		}
		else if (data.rowtypes[rowIdx] == 'R')
		{
			// For the range constraint, we need two linear inequalities:
			// rhs - range <= lhs <= rhs

			if (data.rhsrange[rowIdx] < 0.0)
			{
				printf(FJ_LOG_PREFIX "unsupported negative range value '%g'.\n",
					data.rhsrange[rowIdx]);
				return false;
			}

			solver.addConstraint(RowType::Gte,
				data.rhs[rowIdx] - data.rhsrange[rowIdx],
				data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
				&data.colIdxs[data.rowStart[rowIdx]],
				&data.colCoeffs[data.rowStart[rowIdx]],
				relaxContinuous);
			solver.addConstraint(RowType::Lte,
				data.rhs[rowIdx],
				data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
				&data.colIdxs[data.rowStart[rowIdx]],
				&data.colCoeffs[data.rowStart[rowIdx]],
				relaxContinuous);
			continue;
		}
		else
		{
			printf(FJ_LOG_PREFIX "unsupported constraint type '%c'. Ignoring constraint.\n", data.rowtypes[rowIdx]);
			return false;
		}

		solver.addConstraint(rowtype,
			data.rhs[rowIdx],
			data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
			&data.colIdxs[data.rowStart[rowIdx]],
			&data.colCoeffs[data.rowStart[rowIdx]],
			relaxContinuous);
	}

	return true;
}

// An object containing a function to be executed when the object is destructed.
struct Defer
{
	std::function< void(void) > func;
	Defer(std::function< void(void) > pFunc) : func(pFunc)
	{
	};
	~Defer()
	{
		func();
	}
};

void mapHeuristicSolution(FJStatus& status, bool usePresolved)
{

	Solution s;
	bool conversionOk = false;
	{
		printf(FJ_LOG_PREFIX "received a solution from non-presolved instance.\n");
		s.assignment = std::vector< double >(status.solution, status.solution + status.numVars);
		s.includesContinuous = true;
		conversionOk = true;
	}

	if (conversionOk)
	{
		{
			std::lock_guard< std::mutex > guard(heuristicSolutions_mutex);
			heuristicSolutions.push_back(s);
			totalNumSolutionsFound += 1;
		}


		// FIXME: COPT set solution?
	}
}

const size_t maxEffort = 1 << 31;

// Starts background threads running the Feasibility Jump heuristic.
// Also installs the check-time callback to report any feasible solutions
// back to the MIP solver.
void start_feasibility_jump_heuristic(Model problem,
	size_t maxTotalSolutions,
	bool heuristicOnly,
	bool relaxContinuous = false,
	bool exponentialDecay = false,
	int verbose = 0)
{
	gFJData.originalProblemCopy = problem.Clone();
	{
		auto allThreadsTerminated = std::make_shared< Defer >([]()
		{ heuristicFinished = true; });

		for (int thread_idx = 0; thread_idx < NUM_THREADS; thread_idx += 1)
		{
			auto seed = thread_idx;
			bool usePresolved = thread_idx % 2 == 1;
			double decayFactor = (!exponentialDecay) ? 1.0 : 0.9999;

			std::thread(
				[verbose, maxTotalSolutions, usePresolved, seed,
					relaxContinuous, decayFactor, allThreadsTerminated]()
				{
				  // Prepare data for the non-presolved version.
				  {
					  std::lock_guard< std::mutex > guard(nonPresolvedProblem_mutex);
					  if (gFJData.originalData.numCols == 0)
					  {
						  gFJData.originalData = getCOPTProblemData(gFJData.originalProblemCopy);
					  }
				  }

				  ProblemInstance& data = gFJData.originalData;
				  FeasibilityJumpSolver solver(seed, verbose, decayFactor);
				  bool copyOk = copyDataToHeuristicSolver(solver, data, relaxContinuous);
				  if (!copyOk)
					  return;

				  solver.solve(
					  nullptr, [maxTotalSolutions, usePresolved](FJStatus status) -> CallbackControlFlow
					  {

						double time = std::chrono::duration_cast< std::chrono::milliseconds >(
							std::chrono::steady_clock::now() - startTime).count() / 1000.0;

						// If we received a solution, put it on the queue.
						if (status.solution != nullptr)
						{
							printf("FJSOL %g %g\n", time, status.solutionObjectiveValue);

							mapHeuristicSolution(status, usePresolved);
						}

						// If we have enough solutions or spent enough time, quit.
						auto quitNumSol = totalNumSolutionsFound >= maxTotalSolutions;
						if (quitNumSol)
							printf(FJ_LOG_PREFIX "quitting because number of solutions %zd >= %zd.\n",
								totalNumSolutionsFound.load(),
								maxTotalSolutions);
						auto quitEffort = status.effortSinceLastImprovement > maxEffort;
						if (quitEffort)
							printf(FJ_LOG_PREFIX "quitting because effort %ld > %ld.\n",
								status.effortSinceLastImprovement,
								maxEffort);

						auto quit = quitNumSol || quitEffort || heuristicFinished;
						if (quit)
							printf(FJ_LOG_PREFIX "effort rate: %g Mops/sec\n", status.totalEffort / time / 1.0e6);
						return quit ? CallbackControlFlow::Terminate : CallbackControlFlow::Continue;
					  });
				})
				.detach();
		}
	}

	while (!heuristicFinished)
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(50));
	}
	printf(FJ_LOG_PREFIX "all threads exited.\n");
}

int printUsage()
{
	printf(
		"Usage: fjcopt [--timeout|-t TIMEOUT] [--save-solutions|-s OUTDIR] [--verbose|-v] [--heuristic-only|-h] [--exponential-decay|-e] [--relax-continuous|-r] INFILE\n");
	return 1;
}

int main(int argc, char* argv[])
{
	int verbose = 0;
	bool heuristicOnly = false;
	bool relaxContinuous = false;
	bool exponentialDecay = false;
	int timeout = INT32_MAX / 2;
	size_t maxTotalSolutions = 10;

	std::string inputPath;
	for (int i = 1; i < argc; i += 1)
	{
		std::string argvi(argv[i]);
		if (argvi == "--save-solutions" || argvi == "-s")
		{
			if (i + 1 < argc)
				outDir = std::string(argv[i + 1]);
			else
				return printUsage();
			i += 1;
		}
		else if (argvi == "--timeout" || argvi == "-t")
		{
			if (i + 1 < argc)
				timeout = std::stoi(argv[i + 1]);
			else
				return printUsage();
			i += 1;
		}
		else if (argvi == "--verbose" || argvi == "-v")
			verbose += 1;
		else if (argvi == "--heuristic-only" || argvi == "-h")
			heuristicOnly = true;
		else if (argvi == "--relax-continuous" || argvi == "-r")
			relaxContinuous = true;
		else if (argvi == "--exponential-decay" || argvi == "-e")
			exponentialDecay = true;
		else if (!inputPath.empty())
			return printUsage();
		else
			inputPath = argvi;
	}

	if (inputPath.empty())
		return printUsage();

	inputFilename = inputPath.substr(inputPath.find_last_of("/\\") + 1);

	int returnCode = 0;

	startTime = std::chrono::steady_clock::now();

	Envr env;
	Model problem = env.CreateModel("COPT_FJ");
	problem.Read(inputPath.c_str());
//	problem.Solve();

	start_feasibility_jump_heuristic(problem, maxTotalSolutions, heuristicOnly, relaxContinuous, exponentialDecay, verbose);

	return returnCode;
}
