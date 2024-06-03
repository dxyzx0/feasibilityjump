//
// Created by psw on 6/3/24.
//

#ifndef FEASIBILITYJUMP_HEAURISTIC__PBO_FJ_H_
#define FEASIBILITYJUMP_HEAURISTIC__PBO_FJ_H_

#include <cstdio>
#include <vector>
#include <atomic>
#include <chrono>
#include <mutex>
#include <tuple>
#include <string>
#include <thread>

#include "feasibilityjump.h"
#include "parser/AbcCallback.h"
#include "pbo_fj_pub.h"

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

ProblemInstance getProblemData(AbcCallback& abcCallback);

bool copyDataToHeuristicSolver(FeasibilityJumpSolver& solver, ProblemInstance& data, int relaxContinuous);

void mapHeuristicSolution(FJStatus& status);

void start_feasibility_jump_heuristic(AbcCallback& abcCallback, size_t maxTotalSolutions, bool heuristicOnly, size_t NUM_THREADS,
	bool relaxContinuous = false, bool exponentialDecay = false, int verbose = 0);

int printUsage();

#define CHECK_RETURN(call)                                  \
    do                                                      \
    {                                                       \
        int result_ = call;                                 \
        if (result_ != 0)                                   \
        {                                                   \
            fprintf(stderr, "Line %d: %s failed with %d\n", \
                    __LINE__, #call, result_);              \
            returnCode = result_;                           \
            goto cleanup;                                   \
        }                                                   \
    } while (0)

#endif //FEASIBILITYJUMP_HEAURISTIC__PBO_FJ_H_
