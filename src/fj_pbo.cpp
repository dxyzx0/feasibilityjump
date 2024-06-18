#include <cstdio>
#include <vector>
#include <atomic>
#include <chrono>
#include <mutex>
#include <tuple>
#include <string>
#include <thread>

#include "parser/PboCallback.h"
#include "parser/SimpleParser.h"
#include "fj_solver.h"
#include "fj_pbo.h"

namespace FJ_BIGINT
{

	struct ProblemInstance
	{
		size_t numCols{};
		std::vector< char > varTypes;
		std::vector< IntegerType > lb;
		std::vector< IntegerType > ub;
		std::vector< IntegerType > objCoeffs;

		size_t numRows{};
		size_t numNonZeros{};
		std::vector< char > rowtypes;
		std::vector< IntegerType > rhs;
		std::vector< size_t > rowStart;
		std::vector< size_t > colIdxs;
		std::vector< IntegerType > colCoeffs;
		std::vector< std::tuple< size_t, std::string > > rowRelOp;
	};

	struct FJData
	{
		std::vector< int > originalIntegerCols;
		ProblemInstance originalData;
		ProblemInstance presolvedData;
	};

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

	FJData gFJData;

	std::atomic_size_t global_thread_rank(0);
	std::atomic_size_t totalNumSolutionsFound(0);
	std::atomic_bool heuristicFinished(false);

	std::vector< Solution > heuristicSolutions;
	std::vector< IntegerType > heuristicSolutionValues;
	std::mutex heuristicSolutions_mutex;

	std::mutex nonPresolvedProblem_mutex;

	const size_t maxEffort = std::numeric_limits< size_t >::max();

	ProblemInstance getProblemData(PboCallback& pboCb)
	{
		ProblemInstance data;
		data.numCols = pboCb.getNVar();

		data.varTypes = std::vector< char >(data.numCols);
		std::fill(data.varTypes.begin(), data.varTypes.end(), 'B');

		data.lb = std::vector< IntegerType >(data.numCols, 0);
		data.ub = std::vector< IntegerType >(data.numCols, 1);

		data.objCoeffs = std::vector< IntegerType >(data.numCols);

		for (const auto& i : pboCb.getC())
		{
			size_t idx = std::get< 0 >(i);
			IntegerType val = std::get< 1 >(i);
			data.objCoeffs[idx] = val;
		}

		data.numRows = pboCb.getNCons();

		data.rowtypes = std::vector< char >(data.numRows);
		data.rowRelOp = pboCb.getRelOp();
		for (size_t i = 0; i < data.numRows; i++)
		{
			assert(i == std::get< 0 >(data.rowRelOp[i]));
			std::string relOp = std::get< 1 >(data.rowRelOp[i]);
			if (relOp == "=") data.rowtypes[i] = 'E';
			else if (relOp == ">=") data.rowtypes[i] = 'G';
			else
			{ throw std::runtime_error("Unsupported relation operator"); }
		}

		data.rhs = std::vector< IntegerType >(data.numRows);
		for (size_t i = 0; i < data.numRows; i++)
		{
			assert(i == std::get< 0 >(pboCb.getB()[i]));
			data.rhs[i] = std::get< 1 >(pboCb.getB()[i]);
		}

		data.numNonZeros = pboCb.getA().size();
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "copying %zu x %zu matrix with %zu nonzeros.\n",
			data.numCols, data.numRows, data.numNonZeros);

		data.rowStart = std::vector< size_t >(data.numRows + 1);
		data.colIdxs = std::vector< size_t >(data.numNonZeros);
		data.colCoeffs = std::vector< IntegerType >(data.numNonZeros);
		data.rowStart[0] = 0;

		// add assert to make sure PboCallback.getA() is sorted by row

		size_t i = 0;
		size_t row = 0;
		for (const auto& t : pboCb.getA())
		{
			data.colIdxs[i] = get< 1 >(t);
			data.colCoeffs[i] = get< 2 >(t);
			size_t row_i = get< 0 >(t);
			if (row_i != row)
			{
				data.rowStart[row + 1] = i;
				row = row_i;
			}
			i++;
		}
		data.rowStart[data.numRows] = data.numNonZeros;

		// printf(PBO_LOG_COMMENT_PREFIX "*****************\n");
		// printf(PBO_LOG_COMMENT_PREFIX "numCols= %zu\n", data.numCols);
		// printf(PBO_LOG_COMMENT_PREFIX "numRows= %zu\n", data.numRows);
		// printf(PBO_LOG_COMMENT_PREFIX "numNonZeros= %zu\n", data.numNonZeros);
		// string str = "varTypes: ";
		// for (const char& t : data.varTypes)
		// {
		// 	str += to_string(t) + " ";
		// }
		// printf(PBO_LOG_COMMENT_PREFIX "varTypes= %s\n", str.c_str());

		// str = "lb: ";
		// for (const IntegerType& t : data.lb)
		// {
		// 	str += to_string(t) + " ";
		// }
		// printf(PBO_LOG_COMMENT_PREFIX "lb= %s\n", str.c_str());

		// str = "ub: ";
		// for (const IntegerType& t : data.ub)
		// {
		// 	str += to_string(t) + " ";
		// }
		// printf(PBO_LOG_COMMENT_PREFIX "ub= %s\n", str.c_str());

		// str = "objCoeffs: ";
		// for (const IntegerType& t : data.objCoeffs)
		// {
		// 	str += to_string(t) + " ";
		// }
		// printf(PBO_LOG_COMMENT_PREFIX "objCoeffs= %s\n", str.c_str());

		// str = "rowtypes: ";
		// for (const char& t : data.rowtypes)
		// {
		// 	str += to_string(t) + " ";
		// }
		// printf(PBO_LOG_COMMENT_PREFIX "rowtypes= %s\n", str.c_str());

		// str = "rhs: ";
		// for (const IntegerType& t : data.rhs)
		// {
		// 	str += to_string(t) + " ";
		// }
		// printf(PBO_LOG_COMMENT_PREFIX "rhs= %s\n", str.c_str());

		// str = "rowStart: ";
		// for (const size_t& t : data.rowStart)
		// {
		// 	str += to_string(t) + " ";
		// }
		// printf(PBO_LOG_COMMENT_PREFIX "rowStart= %s\n", str.c_str());

		// str = "colIdxs: ";
		// for (const size_t& t : data.colIdxs)
		// {
		// 	str += to_string(t) + " ";
		// }
		// printf(PBO_LOG_COMMENT_PREFIX "colIdxs= %s\n", str.c_str());

		// str = "colCoeffs: ";
		// for (const IntegerType& t : data.colCoeffs)
		// {
		// 	str += to_string(t) + " ";
		// }
		// printf(PBO_LOG_COMMENT_PREFIX "colCoeffs= %s\n", str.c_str());
		// printf(PBO_LOG_COMMENT_PREFIX "*****************\n");

		return data;
	}

	bool copyDataToHeuristicSolver(FeasibilityJumpSolver& solver, ProblemInstance& data)
	{
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "initializing FJ with %zu vars %zu constraints\n",
			data.numCols,
			data.numRows);
		for (size_t colIdx = 0; colIdx < data.numCols; colIdx++)
		{
			VarType vartype = VarType::Continuous;
			if (data.varTypes[colIdx] == 'B')
			{
				vartype = VarType::Integer;
			}
			else
			{
				printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "unsupported variable type '%c' (%d).\n",
					data.varTypes[colIdx], data.varTypes[colIdx]);
				return false;
			}

			solver.addVar(vartype, data.lb[colIdx], data.ub[colIdx], data.objCoeffs[colIdx]);
		}

		for (size_t rowIdx = 0; rowIdx < data.numRows; rowIdx++)
		{
			assert(data.rowtypes[rowIdx] == 'G' || data.rowtypes[rowIdx] == 'E');
			RowType rowtype;
			if (data.rowtypes[rowIdx] == 'G')
			{
				rowtype = RowType::Gte;
			}
			else if (data.rowtypes[rowIdx] == 'E')
			{
				rowtype = RowType::Equal;
			}
			else
			{
				printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "unsupported constraint type '%c'. Ignoring constraint.\n",
					data.rowtypes[rowIdx]);
				return false;
			}

			solver.addConstraint(rowtype,
				data.rhs[rowIdx],
				data.rowStart[rowIdx + 1] - data.rowStart[rowIdx],
				&data.colIdxs[data.rowStart[rowIdx]],
				&data.colCoeffs[data.rowStart[rowIdx]]);
		}

		return true;
	}

	void mapHeuristicSolution(FJStatus& status)
	{

		Solution s;
		{
			printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "received a solution from non-presolved instance.\n");
			s.assignment = std::vector< IntegerType >(status.solution, status.solution + status.numVars);
		}
		{
			std::lock_guard< std::mutex > guard(heuristicSolutions_mutex);
			heuristicSolutions.push_back(s);
			heuristicSolutionValues.push_back(status.solutionObjectiveValue);
			totalNumSolutionsFound++;
		}
	}

	void postProcess()
	{
		assert(heuristicSolutions.size() == heuristicSolutionValues.size());
		// print the best solution
		IntegerType bestObj = PBOINTMAX;
		size_t bestIdx = 0;
		for (size_t i = 0; i < heuristicSolutions.size(); i++)
		{
			if (heuristicSolutionValues[i] < bestObj)
			{
				bestObj = heuristicSolutionValues[i];
				bestIdx = i;
			}
		}
		if (bestObj == PBOINTMAX)
		{
			printf(PBO_LOG_STATUS_PREFIX PBO_STATUS_UNKNOWN "\n");
		}
		else
		{
			printf(PBO_LOG_STATUS_PREFIX PBO_STATUS_SAT "\n");
			printf(PBO_LOG_OBJ_PREFIX "%s\n", to_string(bestObj).c_str());
			printFullSolution(heuristicSolutions[bestIdx]);
		}
	}

// Starts background threads running the Feasibility Jump heuristic.
// Also installs the check-time callback to report any feasible solutions
// back to the MIP solver.
	void startFeasibilityJumpHeuristic(PboCallback& pboCb,
		size_t maxTotalSolutions,
		size_t NUM_THREADS,
		std::chrono::steady_clock::time_point startTime,
		double timeout,
		int verbose
	)
	{
		{
			shared_ptr< Defer > allThreadsTerminated = std::make_shared< Defer >(
				[]()
				{ heuristicFinished = true; });

			printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "starting FJ with %zu threads.\n", NUM_THREADS);
			for (size_t thread_idx = 0; thread_idx < NUM_THREADS; thread_idx++)
			{
				int seed = thread_idx;
				IntegerType decayFactor = 1;
				auto f = [timeout, verbose, maxTotalSolutions, seed, startTime,
					decayFactor, allThreadsTerminated, &pboCb]()
				{
				  // Increment the thread rank.
				  size_t thread_rank = global_thread_rank++;
				  // Prepare data for the non-presolved version.
				  {
					  std::lock_guard< std::mutex > guard(nonPresolvedProblem_mutex);
					  if (gFJData.originalData.numCols == 0)
					  {
						  gFJData.originalData = getProblemData(pboCb);
					  }
				  }

				  ProblemInstance& data = gFJData.originalData;
				  FeasibilityJumpSolver solver(seed, verbose, decayFactor, thread_rank);
				  bool copyOk = copyDataToHeuristicSolver(solver, data);
				  if (!copyOk)
				  {
					  printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: failed to copy data to FJ solver.\n", thread_rank);
					  return;
				  }

				  solver.solve(
					  nullptr,
					  [maxTotalSolutions, thread_rank, timeout, startTime](FJStatus status) -> CallbackControlFlow
					  {

						double time = std::chrono::duration_cast< std::chrono::milliseconds >(
							std::chrono::steady_clock::now() - startTime).count() / 1000.0;

						// If we received a solution, put it on the queue.
						if (status.solution != nullptr)
						{
							string solutionObjectiveValueStr = to_string(status.solutionObjectiveValue);
							printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_SOL_PREFIX "%zu: time= %g obj= %s\n",
								thread_rank,
								time,
								solutionObjectiveValueStr.c_str());
							mapHeuristicSolution(status);
						}

						// If we have enough solutions or spent enough time, quit.
						bool quitNumSol = totalNumSolutionsFound >= maxTotalSolutions;
						if (quitNumSol)
							printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: quitting because number of solutions %zd >= %zd.\n",
								thread_rank,
								totalNumSolutionsFound.load(),
								maxTotalSolutions);
						bool quitEffort = status.effortSinceLastImprovement > maxEffort;
						if (quitEffort)
							printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: quitting because effort %ld > %ld.\n",
								thread_rank,
								status.effortSinceLastImprovement,
								maxEffort);
						time = std::chrono::duration_cast< std::chrono::milliseconds >(
							std::chrono::steady_clock::now() - startTime).count() / 1000.0;
						bool quitTime = time > timeout;
						if (quitTime)
							printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: quitting because time %g > %g.\n",
								thread_rank,
								time,
								timeout);
						bool quit = quitNumSol || quitEffort || quitTime || heuristicFinished;
						if (quit)
						{
							printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: effort rate: %g Mops/sec\n",
								thread_rank, status.totalEffort / time / 1.0e6);
						}
						return quit ? CallbackControlFlow::Terminate : CallbackControlFlow::Continue;
					  });
				};
				std::thread(f).detach();
			}
		}

		// Wait for all threads to finish.
		while (!heuristicFinished)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(50));
		}
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "all threads exited.\n");

		// post-process the solutions
		postProcess();
	}

} // namespace FJ_BIGINT

int printUsage()
{
	printf(PBO_LOG_COMMENT_PREFIX "Usage: pbo_fj [--jobs|-j JOBS] [--timeout|-t TIMEOUT] [--save-solutions|-s OUTDIR] [--verbose|-v] INFILE\n");
	return 1;
}

int runFeasibilityJumpHeuristic(int argc, const char* argv[])
{
	std::string inputFilename;
	std::string outDir;
	size_t maxTotalSolutions = 100;
	size_t NUM_THREADS = 1;
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	double timeout = 1e20;
	int verbose = 0;
	SimpleParser< PboCallback >* parser;

	try
	{
		std::string inputPath;
		for (int i = 1; i < argc; i++)
		{
			std::string argvi(argv[i]);
			if (argvi == "--save-solutions" || argvi == "-s")
			{
				if (i + 1 < argc)
					outDir = std::string(argv[i + 1]);
				else
					return printUsage();
				i++;
			}
			else if (argvi == "--timeout" || argvi == "-t")
			{
				if (i + 1 < argc)
					timeout = std::stod(argv[i + 1]);
				else
					return printUsage();
				i++;
			}
			else if (argvi == "--verbose" || argvi == "-v")
				verbose++;
			else if (argvi == "--jobs" || argvi == "-j")
			{
				if (i + 1 < argc)
					NUM_THREADS = std::stoi(argv[i + 1]);
				else
					return printUsage();
				i++;
			}
			else if (!inputPath.empty())
				return printUsage();
			else
				inputPath = argvi;
		}

		if (inputPath.empty())
			return printUsage();

		inputFilename = inputPath.substr(inputPath.find_last_of("/\\") + 1);

		const string& filename = inputPath;
		// assert filename ends with ".opb" or ".pb"
		printf(PBO_LOG_COMMENT_PREFIX "Parsing file: %s\n", filename.c_str());
		assert(filename.substr(filename.find_last_of('.') + 1) == "opb" ||
			filename.substr(filename.find_last_of('.') + 1) == "pb");
		parser = new SimpleParser< PboCallback >(filename.c_str());

		parser->setAutoLinearize(true);
		parser->parse();
	}
	catch (exception& e)
	{
		printf(PBO_LOG_COMMENT_PREFIX "Exception: %s\n", e.what());
		printf(PBO_LOG_STATUS_PREFIX PBO_STATUS_UNSUPPORTED "\n");
		return -1;
	}

	// start the Feasibility Jump heuristic
	try
	{
		int returnCode = 0;

		FJ_BIGINT::startFeasibilityJumpHeuristic(parser->cb,
			maxTotalSolutions,
			NUM_THREADS,
			startTime,
			timeout,
			verbose);

		return returnCode;
	}
	catch (exception& e)
	{
		printf(PBO_LOG_COMMENT_PREFIX "Exception: %s\n", e.what());
		printf(PBO_LOG_STATUS_PREFIX PBO_STATUS_UNKNOWN "\n");
		return -1;
	}
}
