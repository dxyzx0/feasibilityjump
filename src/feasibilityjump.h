#ifndef FEASIBILITYJUMP_HEAURISTIC__FEASIBILITYJUMP_H_
#define FEASIBILITYJUMP_HEAURISTIC__FEASIBILITYJUMP_H_

#include <algorithm>
#include <functional>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>
#include "type.h"

using namespace std;

#define FJ_LOG_PREFIX "(FJ) "
#define FJ_LOG_SOL_PREFIX "(FJSOL) "

#define PBO_LOG_COMMENT_PREFIX "c "
#define PBO_LOG_STATUS_PREFIX "s "
#define PBO_LOG_OBJ_PREFIX "o "
#define PBO_LOG_SOL_PREFIX "v "

const IntegerType violationTolerance = 0;
const IntegerType equalityTolerance = 0;

struct Solution
{
	std::vector< IntegerType > assignment;
	bool includesContinuous;
};

void printFullSolution(const Solution& s);

void printSolution(const Solution& s);

void printIdxOfOneInSolution(const Solution& s, size_t thread_rank);

void printSolution(IntegerType* s, size_t n);

void printIdxOfOneInSolution(const IntegerType* s, size_t n, size_t thread_rank);

void printVector(const std::vector< IntegerType >& v);

void printIdxOfOneInVector(const std::vector< IntegerType >& v);

// Measures if two Ints are equal within a tolerance of 1.0e-5.
bool eq(IntegerType a, IntegerType b);

enum RowType
{
	Equal,
	Lte,
	Gte,
};

enum VarType
{
	Continuous,
	Integer
};

enum CallbackControlFlow
{
	Terminate,
	Continue,
};

struct FJStatus
{
	size_t totalEffort;
	size_t effortSinceLastImprovement;
	size_t numVars;
	IntegerType solutionObjectiveValue;
	IntegerType* solution;
};

struct IdxCoeff
{
	size_t idx;
	IntegerType coeff;

	IdxCoeff(size_t idx, IntegerType coeff) : idx(idx), coeff(coeff)
	{
	}
};

struct Var
{
	VarType vartype;
	IntegerType lb;
	IntegerType ub;
	IntegerType objectiveCoeff;
	std::vector< IdxCoeff > coeffs;
};

struct Constraint
{
	RowType sense;
	IntegerType rhs;
	std::vector< IdxCoeff > coeffs;
	IntegerType weight;
	IntegerType incumbentLhs;
	size_t violatedIdx;
	IntegerType zero = 0;

	// Computes the constraint's contribution to the feasibility score:
	// If the constraint is satisfied by the given LHS value, returns 0.
	// If the constraint is violated by the given LHS value, returns -|lhs-rhs|.
	IntegerType score(const IntegerType& lhs) const
	{
		if (sense == RowType::Equal)
			return lhs > rhs ? rhs - lhs : lhs - rhs;
		else if (sense == RowType::Lte)
//            return -(std::max(0, lhs - rhs));
			return -(lhs > rhs ? lhs - rhs : zero);
		else
//            return -(std::max(0., rhs - lhs));
			return -(rhs > lhs ? rhs - lhs : zero);
	}
};

// A potential new value for a variable, including its score.
struct Move
{
	IntegerType value;  // FIXME: potential bug?
	IntegerType score;

	static Move undef()
	{
		Move move;
		move.value = -1; // potential value is 0 or 1, so -1 is a good choice
		move.score = PBOINTMIN;
		return move;
	}
};

// Represents a modification of the LHS in a constraint, for a specific
// variable/constraint combination.The `modifyMove` function below is used to
// update the score of a `Move` to reflect the LHS modification.
struct LhsModification
{
	size_t varIdx;
	size_t constraintIdx;
	IntegerType coeff;
	IntegerType oldLhs;
	IntegerType newLhs;
};

// Stores the MIP problem, an incumbent assignment, and the set of constraints
// that are violated in the current incumbent assignment. This set is maintained
// when changes are given to the incumbent assignment using `setValue`.
struct Problem
{
	std::vector< Var > vars;
	std::vector< Constraint > constraints;
	std::vector< IntegerType > incumbentAssignment;
	std::vector< size_t > violatedConstraints;
	bool usedRelaxContinuous = false;

	size_t nNonzeros = 0;
	IntegerType incumbentObjective = 0;  // FIXME

	size_t addVar(VarType vartype, IntegerType lb, IntegerType ub, IntegerType objCoeff)
	{
		size_t idx = vars.size();
		Var var;
		var.vartype = vartype;
		var.lb = lb;
		var.ub = ub;
		var.objectiveCoeff = objCoeff;
		vars.push_back(var);
		incumbentAssignment.push_back(lb);
		return idx;
	}

	size_t addConstraint(RowType sense, IntegerType rhs, size_t numCoeffs, size_t* rowVarIdxs, IntegerType* rowCoeffs,
		int relax_continuous)
	{
		if (relax_continuous)
			usedRelaxContinuous = true;

		// If we are relaxing continuous variables, an equality needs to be split into Gte and Lte.
		assert(relax_continuous == 0);
//        if (relax_continuous > 0 && sense == RowType::Equal)
//            if (std::any_of(rowVarIdxs, rowVarIdxs + numCoeffs,
//                            [&](size_t varIdx) { return vars[varIdx].vartype == VarType::Continuous; })) {
//                addConstraint(RowType::Gte, rhs, numCoeffs, rowVarIdxs, rowCoeffs, relax_continuous);
//                addConstraint(RowType::Lte, rhs, numCoeffs, rowVarIdxs, rowCoeffs, relax_continuous);
//                return INT_MAX;
//            }

		std::vector< IdxCoeff > coeffs;
		for (size_t i = 0; i < numCoeffs; i++)
		{
			if (relax_continuous > 0 && vars[rowVarIdxs[i]].vartype == VarType::Continuous)
			{
				if (sense == RowType::Lte)
				{
					if (rowCoeffs[i] >= 0)
						rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].lb;
					else
						rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].ub;
				}
				else if (sense == RowType::Gte)
				{
					if (rowCoeffs[i] >= 0)
						rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].ub;
					else
						rhs -= rowCoeffs[i] * vars[rowVarIdxs[i]].lb;
				}
				else
					return PBOINTMIN;
			}
			else
				coeffs.emplace_back(rowVarIdxs[i], rowCoeffs[i]);
		}

		if (coeffs.empty())
		{
			bool ok;
			if (sense == RowType::Lte)
				ok = 0 <= rhs + equalityTolerance;
			else if (sense == RowType::Gte)
				ok = 0 + equalityTolerance >= rhs;
			else
				ok = eq(0, rhs);

			return ok ? PBOINTMAX : PBOINTMIN;
		}

		size_t newConstraintIdx = constraints.size();
		for (IdxCoeff& c : coeffs)
		{
			vars[c.idx].coeffs.emplace_back(newConstraintIdx, c.coeff);
		}

		nNonzeros += coeffs.size();
		Constraint newConstraint;
		newConstraint.coeffs = coeffs;
		newConstraint.incumbentLhs = 0;  // FIXME
		newConstraint.violatedIdx = -1;
		newConstraint.rhs = rhs;
		newConstraint.sense = sense;
		newConstraint.weight = 1;

		constraints.push_back(newConstraint);
		return newConstraintIdx;
	}

	void resetIncumbent(IntegerType* initialValues)
	{
		// Set the initial values, if given.

		if (initialValues)
			for (size_t i = 0; i < vars.size(); i++)
				incumbentAssignment[i] = initialValues[i];
		// std::copy(initialValues, initialValues + vars.size(), incumbentAssignment);

		// Reset the incumbent objective.
		incumbentObjective = 0;
		for (size_t i = 0; i < vars.size(); i++)
			incumbentObjective += vars[i].objectiveCoeff * incumbentAssignment[i];

		// Reset the constraint LHSs and the violatedConstraints list.
		violatedConstraints.clear();
		for (size_t cIdx = 0; cIdx < constraints.size(); cIdx++)
		{
			Constraint& cstr = constraints[cIdx];

			cstr.incumbentLhs = 0;
			for (IdxCoeff& vc : cstr.coeffs)
				cstr.incumbentLhs += vc.coeff * incumbentAssignment[vc.idx];

			if (initialValues == nullptr)
				assert(cstr.incumbentLhs == 0);

			if (cstr.score(cstr.incumbentLhs) < -violationTolerance)
			{
				cstr.violatedIdx = violatedConstraints.size();
				violatedConstraints.push_back(cIdx);
			}
			else
				cstr.violatedIdx = -1;
		}
	}

	// Updates a variable assignment for `varIdx` to `newValue`.
	// Takes a function parameter f that receives a LhsModification
	// for every variable/constraint combination (except for `varIdx` itself)
	// where the LHS of the constraint has changed.
	template< typename F >
	size_t setValue(size_t varIdx, IntegerType newValue, F f)
	{
		size_t dt = 0;
		IntegerType oldValue = incumbentAssignment[varIdx];
		IntegerType delta = (newValue - oldValue);
		incumbentAssignment[varIdx] = newValue;
		incumbentObjective += vars[varIdx].objectiveCoeff * delta;
		// printf("Setting v%d to from %g to value %g\n", varIdx, oldValue, newValue);

		// Update the LHSs of all involved constraints.
		for (IdxCoeff& cstrCoeff : vars[varIdx].coeffs)
		{
			IntegerType oldLhs = constraints[cstrCoeff.idx].incumbentLhs;
			IntegerType newLhs = oldLhs + cstrCoeff.coeff * delta;
			constraints[cstrCoeff.idx].incumbentLhs = newLhs;
			IntegerType newCost = constraints[cstrCoeff.idx].score(newLhs);

			// Add/remove from the violatedConstraints list.
			if (newCost < -violationTolerance && constraints[cstrCoeff.idx].violatedIdx == -1)
			{
				// Became violated.
				constraints[cstrCoeff.idx].violatedIdx = violatedConstraints.size();
				violatedConstraints.push_back(cstrCoeff.idx);
			}
			if (newCost >= -violationTolerance && constraints[cstrCoeff.idx].violatedIdx != -1)
			{
				// Became satisfied.
				size_t lastViolatedIdx = violatedConstraints.size() - 1;
				size_t lastConstraintIdx = violatedConstraints[lastViolatedIdx];
				size_t thisViolatedIdx = constraints[cstrCoeff.idx].violatedIdx;
				std::swap(violatedConstraints[thisViolatedIdx], violatedConstraints[lastViolatedIdx]);
				constraints[lastConstraintIdx].violatedIdx = thisViolatedIdx;
				constraints[cstrCoeff.idx].violatedIdx = -1;
				violatedConstraints.pop_back();
			}

			// Now, report the changes in LHS for other variables.
			dt += constraints[cstrCoeff.idx].coeffs.size();
			for (IdxCoeff& varCoeff : constraints[cstrCoeff.idx].coeffs)
			{
				if (varCoeff.idx != varIdx)
				{
					LhsModification m;
					m.varIdx = varCoeff.idx;
					m.constraintIdx = cstrCoeff.idx;
					m.coeff = varCoeff.coeff;
					m.oldLhs = oldLhs;
					m.newLhs = newLhs;
					f(m);
				}
			}
		}

		return dt;
	}
};

void modifyMove(LhsModification mod, Problem& problem, Move& move);

bool check_feasibility(const IntegerType* solution, const Problem& problem);

// Stores current moves and computes updated jump values for
// the "Jump" move type.
class JumpMove
{
	std::vector< Move > moves;
	std::vector< std::pair< IntegerType, IntegerType>> bestShiftBuffer;

 public:
	void init(Problem& problem)
	{
		moves.resize(problem.vars.size());
	}

	template< typename F >
	void forEachVarMove(size_t varIdx, F f)
	{
		f(moves[varIdx]);
	}

	void updateValue(Problem& problem, size_t varIdx)
	{
		bestShiftBuffer.clear();
		IntegerType varIncumbentValue = problem.incumbentAssignment[varIdx];
		IntegerType currentValue = problem.vars[varIdx].lb;
		IntegerType currentScore = 0;
		IntegerType currentSlope = 0;

		// printf(" updatevalue lb %g ub %g numcells %d\n",
		//        problem.vars[varIdx].lb,
		//        problem.vars[varIdx].ub, problem.vars[varIdx].coeffs.size());

		for (IdxCoeff& cell : problem.vars[varIdx].coeffs)
		{
			Constraint& constraint = problem.constraints[cell.idx];

			std::vector< std::pair< IntegerType, IntegerType > > constraintBounds;
			assert(constraint.sense != RowType::Lte);
			if (constraint.sense == RowType::Lte)
				constraintBounds.emplace_back(PBOINTMIN, constraint.rhs);
			else if (constraint.sense == RowType::Gte)
				constraintBounds.emplace_back(constraint.rhs, PBOINTMAX);
			else
			{
				constraintBounds.emplace_back(PBOINTMIN, constraint.rhs);
				constraintBounds.emplace_back(constraint.rhs, constraint.rhs);  // FIXME: ???
				constraintBounds.emplace_back(constraint.rhs, PBOINTMAX);
			}

			for (const std::pair< IntegerType, IntegerType >& bound : constraintBounds)
			{
				IntegerType residualIncumbent = constraint.incumbentLhs - cell.coeff * varIncumbentValue;

				std::pair< IntegerType, IntegerType > validRange = {
					(bound.first - residualIncumbent + cell.coeff - 1) / cell.coeff,  // Round up  FIXME: correct?
					(bound.second - residualIncumbent) / cell.coeff,  // Round down
				};

				if (problem.vars[varIdx].vartype == VarType::Integer)
					validRange = {
						validRange.first - equalityTolerance,
						validRange.second + equalityTolerance,
					};

				if (validRange.first > validRange.second)
					continue;

				if (validRange.first > currentValue)
				{
					currentScore += constraint.weight * (validRange.first - currentValue);
					currentSlope -= constraint.weight;
					if (validRange.first < problem.vars[varIdx].ub)
						bestShiftBuffer.emplace_back(validRange.first, constraint.weight);
				}

				if (validRange.second <= currentValue)
				{
					currentScore += constraint.weight * (validRange.second - currentValue);
					currentSlope += constraint.weight;
				}
				else if (validRange.second < problem.vars[varIdx].ub)
					bestShiftBuffer.emplace_back(validRange.second, constraint.weight);
			}
		}

		bestShiftBuffer.emplace_back(problem.vars[varIdx].lb, 0);
		bestShiftBuffer.emplace_back(problem.vars[varIdx].ub, 0);
		std::sort(bestShiftBuffer.begin(), bestShiftBuffer.end());

		IntegerType bestScore = currentScore;
		IntegerType bestValue = currentValue;
		// printf("evaluating best shift buffer size %d \n", bestShiftBuffer.size());
		for (const std::pair< IntegerType, IntegerType >& item : bestShiftBuffer)
		{
			currentScore += (item.first - currentValue) * currentSlope;
			currentSlope += item.second;
			currentValue = item.first;

			// printf("bestshift cscore %g cslope %g cval %g bestval %g bestscore %g\n",
			// currentScore,currentSlope, currentValue, bestScore, bestValue
			// );

			if (eq(bestValue, problem.incumbentAssignment[varIdx]) ||
				(!eq(currentValue, problem.incumbentAssignment[varIdx]) && currentScore < bestScore))
			{
				bestScore = currentScore;
				bestValue = currentValue;
			}

			// Slope is always increasing, so if we have a valid value, we can quit
			// as soon as the slope turns nonnegative, since we must already have
			// visited the minimum.
			if (!eq(bestValue, problem.incumbentAssignment[varIdx]) && currentSlope >= 0)
				break;
		}

		// printf("Setting jump for %d to from %g to %g\n", varIdx, problem.incumbentAssignment[varIdx], moves[varIdx].value);
		moves[varIdx].value = bestValue;
		assert(bestValue != -1);
	}
};

class FeasibilityJumpSolver
{
	int verbosity;
	Problem problem;
	JumpMove jumpMove;

	std::vector< size_t > goodVarsSet;
	std::vector< size_t > goodVarsSetIdx;

	std::mt19937 rng;

	IntegerType bestObjective = PBOINTMAX;
	IntegerType objectiveWeight = 0;
	size_t bestViolationScore = SIZE_MAX;
	size_t effortAtLastCallback = 0;
	size_t effortAtLastImprovement = 0;
	size_t totalEffort = 0;

	IntegerType weightUpdateDecay;
	IntegerType weightUpdateIncrement = 1;

	size_t nBumps;

	int seed;

	size_t thread_rank;

	// The probability of choosing a random positive-score variable.
	const double randomVarProbability = 0.001;

	// The probability of choosing a variable using a random constraint's
	// non-zero coefficient after updating weights.

	const double randomCellProbability = 0.01;

	// The number of moves to evaluate, if there are many positive-score
	// variables available.
	const size_t maxMovesToEvaluate = 25;

 public:
	explicit FeasibilityJumpSolver(int _seed = 0,
		int _verbosity = 0,
		IntegerType _weightUpdateDecay = 1,
		size_t _thread_rank = 0)
	{
		verbosity = _verbosity;
		weightUpdateDecay = _weightUpdateDecay;
		thread_rank = _thread_rank;
		seed = _seed;
		rng = std::mt19937(seed);
		nBumps = 0;
	}

	size_t addVar(VarType vartype, IntegerType lb, IntegerType ub, IntegerType objCoeff)
	{
		goodVarsSetIdx.push_back(-1);
		return problem.addVar(vartype, lb, ub, objCoeff);
	}

	size_t addConstraint(RowType sense, IntegerType rhs, size_t numCoeffs, size_t* rowVarIdxs, IntegerType* rowCoeffs,
		int relax_continuous)
	{
		return problem.addConstraint(sense, rhs, numCoeffs, rowVarIdxs, rowCoeffs, relax_continuous);
	}

	int solve(IntegerType* initialValues, const std::function< CallbackControlFlow(FJStatus) >& callback)
	{
		assert(callback);
		if (verbosity >= 1)
#ifdef useGMP
			printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: starting solve. weightUpdateDecay=%s, relaxContinuous=%d, seed=%d \n",
				thread_rank, weightUpdateDecay.get_str().c_str(), problem.usedRelaxContinuous, seed);
#else
			printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: starting solve. weightUpdateDecay=%ld, relaxContinuous=%d, seed=%d \n",
				thread_rank, weightUpdateDecay, problem.usedRelaxContinuous, seed);
#endif

		init(initialValues);

		for (size_t step = 0; step < PBOINTMAX; step++)
		{
			if (user_terminate(callback, nullptr))
				break;

			if (step % 1000 == 0)
			{
				if (verbosity >= 1)
					printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: step %zu viol %zd good %zd bumps %zd\n", thread_rank, step,
						problem.violatedConstraints.size(), goodVarsSet.size(), nBumps);
			}

			if (problem.violatedConstraints.size() < bestViolationScore)
			{
				effortAtLastImprovement = totalEffort;
				bestViolationScore = problem.violatedConstraints.size();
			}

			if (problem.violatedConstraints.empty() && problem.incumbentObjective < bestObjective)
			{
				effortAtLastImprovement = totalEffort;
				bestObjective = problem.incumbentObjective;
				if (user_terminate(callback, problem.incumbentAssignment.data()))
					break;
			}

			if (!problem.violatedConstraints.empty())
			{
				assert(!check_feasibility(problem.incumbentAssignment.data(), problem));
			}

			if (problem.vars.empty())
				break;

			size_t var = selectVariable();
			doVariableMove(var);
		}

		return 0;
	}

 private:
	void init(IntegerType* initialValues)
	{
		problem.resetIncumbent(initialValues);
		jumpMove.init(problem);
		totalEffort += problem.nNonzeros;

		// Reset the variable scores.
		goodVarsSet.clear();
		for (size_t i = 0; i < problem.vars.size(); i++)
			resetMoves(i);
	}

	size_t selectVariable()
	{
		if (!goodVarsSet.empty())
		{
			if (std::uniform_real_distribution< double >(0., 1.)(rng) < randomVarProbability)
				return goodVarsSet[rng() % goodVarsSet.size()];

			size_t sampleSize = std::min(maxMovesToEvaluate, goodVarsSet.size());
			totalEffort += sampleSize;

			IntegerType bestScore = PBOINTMIN;
			size_t bestVar = PBOINTMAX;
			for (size_t i = 0; i < sampleSize; i++)
			{
				size_t setidx = rng() % goodVarsSet.size();
				size_t varIdx = goodVarsSet[setidx];
				// assert(goodVarsSetIdx[varIdx] >= 0 && goodVarsSetIdx[varIdx] == setidx);
				Move move = bestMove(varIdx);
				// assert(move.score > equalityTolerance);
				if (move.score > bestScore)
				{
					bestScore = move.score;
					bestVar = varIdx;
				}
			}
			assert(bestVar != PBOINTMAX);
			return bestVar;
		}

		// Local minimum, update weights.
		updateWeights();

		if (!problem.violatedConstraints.empty())
		{
			size_t cstrIdx = problem.violatedConstraints[rng() % problem.violatedConstraints.size()];
			Constraint& constraint = problem.constraints[cstrIdx];

			if (std::uniform_real_distribution< double >(0., 1.)(rng) < randomCellProbability)
				return constraint.coeffs[rng() % constraint.coeffs.size()].idx;

			IntegerType bestScore = PBOINTMIN;
			size_t bestVarIdx = PBOINTMAX;

			for (IdxCoeff& cell : constraint.coeffs)
			{
				Move move = bestMove(cell.idx);
				if (move.score > bestScore)
				{
					bestScore = move.score;
					bestVarIdx = cell.idx;
				}
			}
			return bestVarIdx;
		}

		// Fallback to random choice.
		return rng() % problem.vars.size();
	}

	void updateWeights()
	{
		if (verbosity >= 2)
			printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: Reached a local minimum.\n", thread_rank);
		nBumps++;
		bool rescaleAllWeights = false;
		size_t dt = 0;

		if (problem.violatedConstraints.empty())
		{
			objectiveWeight += weightUpdateIncrement;
			if (objectiveWeight > (1L << 40))
				rescaleAllWeights = true;

			dt += problem.vars.size();
			for (size_t varIdx = 0; varIdx < problem.vars.size(); varIdx++)
				forEachMove(
					varIdx, [&](Move& move)
					{
					  assert(move.value != -1);
					  move.score += weightUpdateIncrement *
						  problem.vars[varIdx].objectiveCoeff *
						  (move.value - problem.incumbentAssignment[varIdx]);
					});
		}
		else
		{
			for (size_t& cIdx : problem.violatedConstraints)
			{
				Constraint& constraint = problem.constraints[cIdx];
				constraint.weight += weightUpdateIncrement;
				if (constraint.weight > (1L << 40))
					rescaleAllWeights = true;

				dt += constraint.coeffs.size();
				for (IdxCoeff& cell : constraint.coeffs)
				{
					forEachMove(
						cell.idx, [&](Move& move)
						{
						  assert(move.value != -1);
						  IntegerType candidateLhs = constraint.incumbentLhs + cell.coeff * (move.value -
							  problem.incumbentAssignment[cell.idx]);
						  IntegerType diff = weightUpdateIncrement * (constraint.score(candidateLhs) -
							  constraint.score(constraint.incumbentLhs));
						  move.score += diff;
						});

					updateGoodMoves(cell.idx);
				}
			}
		}

		weightUpdateIncrement /= weightUpdateDecay;
		if (rescaleAllWeights)
		{
			weightUpdateIncrement = 1;  // FIXME: 0 or 1?
			objectiveWeight = 0;

			for (Constraint& c : problem.constraints)
				c.weight = 1;
			dt += problem.constraints.size();

			for (size_t i = 0; i < problem.vars.size(); i++)
				resetMoves(i);
		}

		totalEffort += dt;
	}

	Move bestMove(size_t varIdx)
	{
		Move best = Move::undef();
		forEachMove(varIdx, [&](Move& move)
		{
		  if (move.score > best.score)
			  best = move;
		});
		assert(best.value != -1);
		return best;
	}

	void doVariableMove(size_t varIdx)
	{
		// First, we get the best move for the variable;
		Move m = bestMove(varIdx);
		IntegerType newValue = m.value;
		assert(newValue != -1);
		// assert(!isnan(newValue));

		// Update the incumbent solution.
		// printf("Setting var %d from %g to %g for a score of %g\n", varIdx, oldValue, newValue, m.score);

		totalEffort += problem.setValue(
			varIdx, newValue, [&](LhsModification mod)
			{
			  forEachMove(mod.varIdx, [&](Move& m)
			  { modifyMove(mod, problem, m); });
			  updateGoodMoves(mod.varIdx);
			});

		resetMoves(varIdx);
	}

	void updateGoodMoves(size_t varIdx)
	{
		bool anyGoodMoves = bestMove(varIdx).score > 0;
		if (anyGoodMoves && goodVarsSetIdx[varIdx] == -1)
		{
			// Became good, add to good set.
			goodVarsSetIdx[varIdx] = goodVarsSet.size();
			goodVarsSet.push_back(varIdx);
		}
		else if (!anyGoodMoves && goodVarsSetIdx[varIdx] != -1)
		{
			// Became bad, remove from good set.
			size_t lastSetIdx = goodVarsSet.size() - 1;
			size_t lastVarIdx = goodVarsSet[lastSetIdx];
			size_t thisSetIdx = goodVarsSetIdx[varIdx];
			std::swap(goodVarsSet[thisSetIdx], goodVarsSet[lastSetIdx]);
			goodVarsSetIdx[lastVarIdx] = thisSetIdx;
			goodVarsSetIdx[varIdx] = -1;
			goodVarsSet.pop_back();
		}
	}

	template< typename F >
	void forEachMove(size_t varIdx, F f)
	{
		jumpMove.forEachVarMove(varIdx, f);

		// TODO: here, we can add more move types.
		// upDownMove.forEachVarMove(varIdx, f);
	}

	void resetMoves(size_t varIdx)
	{
		totalEffort += problem.vars[varIdx].coeffs.size();
		jumpMove.updateValue(problem, varIdx);

		forEachMove(
			varIdx, [&](Move& move)
			{
			  assert(move.value != -1);
			  move.score = 0;
			  move.score += objectiveWeight *
				  problem.vars[varIdx].objectiveCoeff *
				  (move.value - problem.incumbentAssignment[varIdx]);

			  for (IdxCoeff& cell : problem.vars[varIdx].coeffs)
			  {
				  Constraint& constraint = problem.constraints[cell.idx];
				  IntegerType candidateLhs = constraint.incumbentLhs +
					  cell.coeff *
						  (move.value - problem.incumbentAssignment[varIdx]);
				  move.score += constraint.weight *
					  (constraint.score(candidateLhs) - constraint.score(constraint.incumbentLhs));
			  }
			});

		updateGoodMoves(varIdx);
	}

	bool user_terminate(const std::function< CallbackControlFlow(FJStatus) >& callback, IntegerType* solution)
	{
		const size_t CALLBACK_EFFORT = 500000;
		if (solution != nullptr || totalEffort - effortAtLastCallback > CALLBACK_EFFORT)
		{
			if (verbosity >= 2)
				printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: calling user termination.\n", thread_rank);
			effortAtLastCallback = totalEffort;

			FJStatus status{};
			status.totalEffort = totalEffort;
			status.effortSinceLastImprovement = totalEffort - effortAtLastImprovement;

			status.solution = solution;
			if (solution != nullptr)
			{
				//check if the solution is feasible
				if (!check_feasibility(solution, problem))
				{
					throw std::runtime_error("Solution is not feasible.");
				}
				printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: solution is feasible.\n", thread_rank);
				printIdxOfOneInSolution(solution, problem.vars.size(), thread_rank);
			}
			else
				printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: no solution.\n", thread_rank);
			status.numVars = problem.vars.size();
			status.solutionObjectiveValue = problem.incumbentObjective;

			CallbackControlFlow result = callback(status);
			if (result == CallbackControlFlow::Terminate)
			{
				if (verbosity >= 2)
					printf(PBO_LOG_COMMENT_PREFIX FJ_LOG_PREFIX "%zu: quitting.\n", thread_rank);
				return true;
			}
		}
		return false;
	}
};

#endif //FEASIBILITYJUMP_HEAURISTIC__FEASIBILITYJUMP_H_