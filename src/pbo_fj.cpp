#include "pbo_fj.h"
#include "parser/SimpleParser.h"

std::atomic_size_t global_thread_rank(0);
std::atomic_size_t totalNumSolutionsFound(0);
std::atomic_size_t totalNumSolutionsAdded(0);
std::atomic_bool presolveFinished(false);
std::atomic_bool heuristicFinished(false);
std::chrono::steady_clock::time_point startTime;

std::vector< Solution > heuristicSolutions;
std::vector< IntegerType > heuristicSolutionValues;
std::mutex heuristicSolutions_mutex;

std::mutex presolvedProblem_mutex;
std::mutex nonPresolvedProblem_mutex;

FJData gFJData;

std::string inputFilename;
std::string outDir;

const size_t maxEffort = 99999999999L;

ProblemInstance getProblemData(AbcCallback& abcCallback)
{
	ProblemInstance data;
	data.numCols = abcCallback.getNVar();

	data.varTypes = std::vector< char >(data.numCols);
	std::fill(data.varTypes.begin(), data.varTypes.end(), 'B');

	data.lb = std::vector< IntegerType >(data.numCols, 0);
	data.ub = std::vector< IntegerType >(data.numCols, 1);

	data.objCoeffs = std::vector< IntegerType >(data.numCols);

	for (const auto& i : abcCallback.getC())
	{
		size_t idx = std::get< 0 >(i);
		IntegerType val = std::get< 1 >(i);
		data.objCoeffs[idx] = val;
	}

	data.numRows = abcCallback.getNCons();

	data.rowtypes = std::vector< char >(data.numRows);
	data.rowRelOp = abcCallback.getRelOp();
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
		assert(i == std::get< 0 >(abcCallback.getB()[i]));
		data.rhs[i] = std::get< 1 >(abcCallback.getB()[i]);
	}

	data.numNonZeros = abcCallback.getA().size();
	printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "copying %zu x %zu matrix with %zu nonzeros.\n",
		data.numCols, data.numRows, data.numNonZeros);

	data.rowStart = std::vector< size_t >(data.numRows + 1);
	data.colIdxs = std::vector< size_t >(data.numNonZeros);
	data.colCoeffs = std::vector< IntegerType >(data.numNonZeros);
	data.rowStart[0] = 0;

	// add assert to make sure abcCallback.getA() is sorted by row


	size_t i = 0;
	size_t row = 0;
	for (const auto& t : abcCallback.getA())
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

//	cout << "*****************\n" << endl;
//	cout << "numCols= \n" << data.numCols << endl;
//	cout << "numRows= \n" << data.numRows << endl;
//	cout << "numNonZeros= \n" << data.numNonZeros << endl;
//	cout << " varTypes: " << endl;
//	for (const char& t : data.varTypes)
//	{
//		std::cout << t << " ";
//	}
//	std::cout << std::endl;
//	cout << " lb: " << endl;
//	for (const IntegerType& t : data.lb)
//	{
//		std::cout << t << " ";
//	}
//	std::cout << std::endl;
//	cout << " ub: " << endl;
//	for (const IntegerType& t : data.ub)
//	{
//		std::cout << t << " ";
//	}
//	std::cout << std::endl;
//	cout << " objCoeffs: " << endl;
//	for (const IntegerType& t : data.objCoeffs)
//	{
//		std::cout << t << " ";
//	}
//	std::cout << std::endl;
//	cout << " rowtypes: " << endl;
//	for (const char& t : data.rowtypes)
//	{
//		std::cout << t << " ";
//	}
//	std::cout << std::endl;
//	cout << " rhs:" << endl;
//	for (const IntegerType& t : data.rhs)
//	{
//		std::cout << t << " ";
//	}
//	std::cout << std::endl;
//	cout << " rowStart:" << endl;
//	for (const size_t& t : data.rowStart)
//	{
//		std::cout << t << " ";
//	}
//	std::cout << std::endl;
//	cout << " colIdxs:" << endl;
//	for (const size_t& t : data.colIdxs)
//	{
//		std::cout << t << " ";
//	}
//	std::cout << std::endl;
//	cout << " colCoeffs:" << endl;
//	for (const IntegerType& t : data.colCoeffs)
//	{
//		std::cout << t << " ";
//	}
//	std::cout << std::endl;
//
//	cout << "*****************\n" << endl;

	return data;
}

bool copyDataToHeuristicSolver(FeasibilityJumpSolver& solver, ProblemInstance& data, int relaxContinuous)
{
	printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "initializing FJ with %zu vars %zu constraints\n", data.numCols, data.numRows);
	for (size_t colIdx = 0; colIdx < data.numCols; colIdx += 1)
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

	for (size_t rowIdx = 0; rowIdx < data.numRows; rowIdx += 1)
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
			printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "unsupported constraint type '%c'. Ignoring constraint.\n", data.rowtypes[rowIdx]);
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

void mapHeuristicSolution(FJStatus& status)
{

	Solution s;
	bool conversionOk = false;
	{
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "received a solution from non-presolved instance.\n");
		s.assignment = std::vector< IntegerType >(status.solution, status.solution + status.numVars);
		s.includesContinuous = false;
		conversionOk = true;
	}

	if (conversionOk)
	{
		{
			std::lock_guard< std::mutex > guard(heuristicSolutions_mutex);
			heuristicSolutions.push_back(s);
            heuristicSolutionValues.push_back(status.solutionObjectiveValue);
			totalNumSolutionsFound += 1;
		}
	}
}

// Starts background threads running the Feasibility Jump heuristic.
// Also installs the check-time callback to report any feasible solutions
// back to the MIP solver.
void start_feasibility_jump_heuristic(AbcCallback& abcCallback,
	size_t maxTotalSolutions,
	bool heuristicOnly,
	size_t NUM_THREADS,
	bool relaxContinuous,
	bool exponentialDecay,
	int verbose)
{
	{
		shared_ptr< Defer > allThreadsTerminated = std::make_shared< Defer >(
			[]()
			{ heuristicFinished = true; });

        printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "starting FJ with %zu threads.\n", NUM_THREADS);
		for (size_t thread_idx = 0; thread_idx < NUM_THREADS; thread_idx += 1)
		{
			int seed = thread_idx;
			bool usePresolved = thread_idx % 2 == 1;
			IntegerType decayFactor = (!exponentialDecay) ? 1 : 1;
			auto f = [verbose, maxTotalSolutions, usePresolved, seed,
				relaxContinuous, decayFactor, allThreadsTerminated, &abcCallback]()
			{
			  // Increment the thread rank.
			  size_t thread_rank = global_thread_rank++;
			  // Prepare data for the non-presolved version.
			  {
				  std::lock_guard< std::mutex > guard(nonPresolvedProblem_mutex);
				  if (gFJData.originalData.numCols == 0)
				  {
					  gFJData.originalData = getProblemData(abcCallback);
				  }
			  }

			  ProblemInstance& data = gFJData.originalData;
			  FeasibilityJumpSolver solver(seed, verbose, decayFactor, thread_rank);
			  bool copyOk = copyDataToHeuristicSolver(solver, data, relaxContinuous);
			  if (!copyOk)
			  {
				  printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: failed to copy data to FJ solver.\n", thread_rank);
				  return;
			  }

			  solver.solve(
				  nullptr, [maxTotalSolutions, usePresolved, thread_rank](FJStatus status) -> CallbackControlFlow
				  {

					double time = std::chrono::duration_cast< std::chrono::milliseconds >(
						std::chrono::steady_clock::now() - startTime).count() / 1000.0;

					// If we received a solution, put it on the queue.
					if (status.solution != nullptr)
					{
#ifdef useGMP
						char* solutionObjectiveValueStr =
							mpz_get_str(nullptr, 10, status.solutionObjectiveValue.get_mpz_t());
#else
						char* solutionObjectiveValueStr = nullptr;
						asprintf(&solutionObjectiveValueStr, "%ld", status.solutionObjectiveValue);
#endif
                        printf(PBO_LOG_COMMENT_PREFIX "(FJSOL) %zu: %g %s\n", thread_rank, time, solutionObjectiveValueStr);
                        free(solutionObjectiveValueStr);
						mapHeuristicSolution(status);
					}

					// If we have enough solutions or spent enough time, quit.
					bool quitNumSol = totalNumSolutionsFound >= maxTotalSolutions;
					if (quitNumSol)
						printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: quitting because number of solutions %zd >= %zd.\n", thread_rank,
							totalNumSolutionsFound.load(), maxTotalSolutions);
					bool quitEffort = status.effortSinceLastImprovement > maxEffort;
					if (quitEffort)
						printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: quitting because effort %ld > %ld.\n", thread_rank,
							status.effortSinceLastImprovement, maxEffort);

					bool quit = quitNumSol || quitEffort || heuristicFinished;
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
    assert(heuristicSolutions.size() == heuristicSolutionValues.size());
    // print the best solution
    IntegerType bestObj = PBOINTMAX;
    size_t bestIdx = 0;
    for (size_t i = 0; i < heuristicSolutions.size(); i += 1)
    {
        if (heuristicSolutionValues[i] < bestObj)
        {
            bestObj = heuristicSolutionValues[i];
            bestIdx = i;
        }
    }
    if (bestObj == PBOINTMAX)
    {
        printf(PBO_LOG_STATUS_PREFIX "UNKNOWN\n");
    }
    else
    {
        printf(PBO_LOG_STATUS_PREFIX "SATISFIABLE\n");
        printf(PBO_LOG_OBJ_PREFIX "%ld\n", bestObj);
        printFullSolution(heuristicSolutions[bestIdx]);
    }
}

int printUsage()
{
	printf(PBO_LOG_COMMENT_PREFIX "Usage: pbo_fj [--jobs|-j JOBS] [--timeout|-t TIMEOUT] [--save-solutions|-s OUTDIR] [--verbose|-v] [--heuristic-only|-h] [--exponential-decay|-e] [--relax-continuous|-r] INFILE\n");
	return 1;
}

int run_feasibility_jump_heuristic(int argc, char* argv[])
{
	int verbose = 0;
	bool heuristicOnly = false;
	bool relaxContinuous = false;
	bool exponentialDecay = false;
	int timeout = INT32_MAX / 2;
	size_t maxTotalSolutions = 5;
	size_t NUM_THREADS = 1;

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
//        else if (argvi == "--relax-continuous" || argvi == "-r")
//            relaxContinuous = true;
//        else if (argvi == "--exponential-decay" || argvi == "-e")
//            exponentialDecay = true;
		else if (argvi == "--jobs" || argvi == "-j")
		{
			if (i + 1 < argc)
				NUM_THREADS = std::stoi(argv[i + 1]);
			else
				return printUsage();
			i += 1;
		}
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

	SimpleParser< AbcCallback >* parser;
	try
	{
		const string& filename = inputPath;
		// assert filename ends with ".opb" or ".pb"
        printf(PBO_LOG_COMMENT_PREFIX "Parsing file: %s\n", filename.c_str());
		assert(filename.substr(filename.find_last_of('.') + 1) == "opb" ||
			filename.substr(filename.find_last_of('.') + 1) == "pb");
		parser = new SimpleParser< AbcCallback >(filename.c_str());

		parser->setAutoLinearize(true);
		parser->parse();

		start_feasibility_jump_heuristic(parser->cb,
			maxTotalSolutions,
			heuristicOnly,
			NUM_THREADS,
			relaxContinuous,
			exponentialDecay,
			verbose);
	}
	catch (exception& e)
	{
		cerr << "ERROR: " << e.what() << endl;
		return -1;
	}

	return returnCode;
}
