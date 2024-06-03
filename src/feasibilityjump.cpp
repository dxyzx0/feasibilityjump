//
// Created by psw on 6/3/24.
//

#include "feasibilityjump.h"

void printSolution(const Solution& s)
{
	cout << "solution: ";
	for (const IntegerType& v : s.assignment)
	{
		cout << v << " ";
	}
	cout << endl;
}

void printIdxOfOneInSolution(const Solution& s, size_t thread_rank)
{
	string str = "solution: ";
	for (size_t i = 0; i < s.assignment.size(); i += 1)
	{
		if (s.assignment[i] == 1)
			str += to_string(i) + " ";
	}
	str += "\n";
	printf(FJ_LOG_PREFIX "%s", str.c_str());
}

void printSolution(IntegerType* s, size_t n)
{
	cout << "solution: ";
	for (size_t i = 0; i < n; i += 1)
	{
		cout << s[i] << " ";
	}
	cout << endl;
}

void printIdxOfOneInSolution(IntegerType* s, size_t n, size_t thread_rank)
{
	string str = "solution: ";
	for (size_t i = 0; i < n; i += 1)
	{
		if (s[i] == 1)
			str += to_string(i) + " ";
	}
	str += "\n";
	printf(FJ_LOG_PREFIX "%zu: %s", thread_rank, str.c_str());

}

void printVector(const std::vector< IntegerType >& v)
{
	cout << "vector: ";
	for (const IntegerType& x : v)
	{
		cout << x << " ";
	}
	cout << endl;
}

void printIdxOfOneInVector(const std::vector< IntegerType >& v)
{
	cout << "vector: ";
	for (size_t i = 0; i < v.size(); i += 1)
	{
		if (v[i] == 1)
			cout << i << " ";
	}
	cout << endl;
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