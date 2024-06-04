//
// Created by psw on 5/19/24.
//

#ifndef FEASIBILITYJUMP_HEURISTIC_DEFAULTCALLBACK_H_
#define FEASIBILITYJUMP_HEURISTIC_DEFAULTCALLBACK_H_

#include <vector>
#include <string>
#include "../fj_type_int.h"

using std::vector;

class DefaultCallback
{
 public:
	virtual void metaData(int nbvar, int nbconstr);
	virtual void beginObjective();
	virtual void endObjective();
	virtual void objectiveTerm(IntegerType coeff, int idVar);
	virtual void objectiveProduct(IntegerType coeff, std::vector< int > list);
	virtual void beginConstraint();
	virtual void endConstraint();
	virtual void constraintTerm(IntegerType coeff, int idVar);
	virtual void constraintProduct(IntegerType coeff, std::vector< int > list);
	virtual void constraintRelOp(std::string relop);
	virtual void constraintRightTerm(IntegerType val);
	virtual void linearizeProduct(int newSymbol, std::vector< int > product);
};

#endif //FEASIBILITYJUMP_HEURISTIC_DEFAULTCALLBACK_H_
