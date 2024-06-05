//
// Created by psw on 6/3/24.
//

#include "fj_solver.h"

namespace FJ_BIGINT
{

#ifdef useGMP
	string to_string(const IntegerType& x)
	{
		return x.get_str();
	}
#endif

	void printFullSolution(const Solution& s)
	{
		string str;
		for (size_t i = 0; i < s.assignment.size(); i++)
		{
			if (s.assignment[i] == 1)
				str += "x" + to_string(i + 1) + " ";
			else
				str += "-x" + to_string(i + 1) + " ";
		}
		printf(PBO_LOG_SOL_PREFIX "%s\n", str.c_str());
	}

	void printSolution(const Solution& s)
	{
		string str = "solution: ";
		for (const IntegerType& v : s.assignment)
		{
			str += to_string(v) + " ";
		}
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_SOL_PREFIX "%s\n", str.c_str());
	}

	void printIdxOfOneInSolution(const Solution& s)
	{
		string str = "solution: ";
		for (size_t i = 0; i < s.assignment.size(); i++)
		{
			if (s.assignment[i] == 1)
				str += to_string(i) + " ";
		}
		str += "\n";
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_SOL_PREFIX "%s", str.c_str());
	}

	void printSolution(IntegerType* s, size_t n)
	{
		string str = "solution: ";
		for (size_t i = 0; i < n; i++)
		{
			str += to_string(s[i]) + " ";
		}
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_SOL_PREFIX "%s\n", str.c_str());
	}

	void printIdxOfOneInSolution(const IntegerType* s, size_t n, size_t thread_rank)
	{
		string str;
		for (size_t i = 0; i < n; i++)
		{
			if (s[i] == 1)
				str += to_string(i) + " ";
		}
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_SOL_PREFIX "%zu: solution: %s\n", thread_rank, str.c_str());
	}

	void printVector(const std::vector< IntegerType >& v)
	{
		string str = "vector: ";
		for (const IntegerType& x : v)
		{
			str += to_string(x) + " ";
		}
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%s\n", str.c_str());
	}

	void printIdxOfOneInVector(const std::vector< IntegerType >& v)
	{
		string str = "vector: ";
		for (size_t i = 0; i < v.size(); i++)
		{
			if (v[i] == 1)
				str += to_string(i) + " ";
		}
		printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%s\n", str.c_str());
	}

// Measures if two Ints are equal within a tolerance of 1.0e-5.
	bool eq(IntegerType a, IntegerType b)
	{
		return a == b;
	}

	void modifyMove(LhsModification mod, Problem& problem, Move& move)
	{
		assert(move.value != -1);
		Constraint& c = problem.constraints[mod.constraintIdx];
		IntegerType incumbent = problem.incumbentAssignment[mod.varIdx];
		IntegerType oldModifiedLhs = mod.oldLhs + mod.coeff * (move.value - incumbent);
		IntegerType oldScoreTerm = c.weight * (c.score(oldModifiedLhs) - c.score(mod.oldLhs));
		IntegerType newModifiedLhs = mod.newLhs + mod.coeff * (move.value - incumbent);
		IntegerType newScoreTerm = c.weight * (c.score(newModifiedLhs) - c.score(mod.newLhs));
		move.score += newScoreTerm - oldScoreTerm;
	}

	bool checkSol(const IntegerType* solution, const Problem& problem)
	{
		// Check if the solution is feasible.
		for (const Constraint& c : problem.constraints)
		{
			IntegerType lhs = 0;
			for (const IdxCoeff& cell : c.coeffs)
				lhs += cell.coeff * solution[cell.idx];
			if (c.sense == RowType::Equal && !eq(lhs, c.rhs))
				return false;
			if (c.sense == RowType::Gte && lhs < c.rhs)
				return false;
			if (c.sense == RowType::Lte)
				throw std::runtime_error("Unsupported constraint type.");
		}
		return true;
	}

}