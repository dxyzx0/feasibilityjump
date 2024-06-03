//
// Created by psw on 5/19/24.
//

#include <iostream>
#include <stdexcept>
#include <vector>
#include "DefaultCallback.h"

using namespace std;

/**
 * defines callback that are used by the parser to tranmit information
 * to the solver
 *
 * These functions must be redefined in a subclass.
 */

/**
 * callback called when we get the number of variables and the
 * expected number of constraints
 * @param nbvar: the number of variables
 * @param nbconstr: the number of contraints
 */
void DefaultCallback::metaData(int nbvar, int nbconstr)
{
//	cout << "[nbvar=" << nbvar << "]" << endl;
//	cout << "[nbconstr=" << nbconstr << "]" << endl;
}

/**
 * callback called before we read the objective function
 */
void DefaultCallback::beginObjective()
{
//	cout << "objective:  ";
}

/**
 * callback called after we've read the objective function
 */
void DefaultCallback::endObjective()
{
//	cout << endl;
}

/**
 * callback called when we read a term of the objective function
 *
 * @param coeff: the coefficient of the term
 * @param idVar: the numerical identifier of the variable
 */
void DefaultCallback::objectiveTerm(IntegerType coeff, int idVar)
{
//	cout << "[" << showpos << coeff << noshowpos << " x" << idVar << "] ";
}

/**
 * callback called when we read a term of the objective function
 * which is a product of literals
 *
 * @param coeff: the coefficient of the term
 * @param list: list of literals which appear in the product
 */
void DefaultCallback::objectiveProduct(IntegerType coeff, vector< int > list)
{
//	cout << "[" << showpos << coeff << noshowpos << " ";
//	for (int i = 0; i < list.size(); ++i)
//	{
//		if (list[i] < 0)
//			cout << "~x" << -list[i] << ' ';
//		else
//			cout << "x" << list[i] << ' ';
//	}
//	cout << "] ";
}

/**
 * callback called before we read a constraint
 */
void DefaultCallback::beginConstraint()
{
//	cout << "constraint: ";
}

/**
 * callback called after we've read a constraint
 */
void DefaultCallback::endConstraint()
{
//	cout << endl;
}

/**
 * callback called when we read a term of a constraint
 *
 * @param coeff: the coefficient of the term
 * @param idVar: the numerical identifier of the variable
 */
void DefaultCallback::constraintTerm(IntegerType coeff, int idVar)
{
//	cout << "[" << showpos << coeff << noshowpos << " x" << idVar << "] ";
}

/**
 * callback called when we read a term of a constraint which is a
 * product of literals
 *
 * @param coeff: the coefficient of the term
 * @param list: list of literals which appear in the product
 */
void DefaultCallback::constraintProduct(IntegerType coeff, vector< int > list)
{
//	cout << "[" << showpos << coeff << noshowpos << " ";
//	for (int i = 0; i < list.size(); ++i)
//	{
//		if (list[i] < 0)
//			cout << "~x" << -list[i] << ' ';
//		else
//			cout << "x" << list[i] << ' ';
//	}
//	cout << "] ";
}

/**
 * callback called when we read the relational operator of a constraint
 *
 * @param relop: the relational oerator (>= or =)
 */
void DefaultCallback::constraintRelOp(string relop)
{
//	cout << "[" << relop << "] ";
}

/**
 * callback called when we read the right term of a constraint (also
 * known as the degree)
 *
 * @param val: the degree of the constraint
 */
void DefaultCallback::constraintRightTerm(IntegerType val)
{
//	cout << "[" << val << "]";
}

/**
 * add the necessary constraints to define newSymbol as equivalent
 * to the product (conjunction) of literals in product.
 */
void DefaultCallback::linearizeProduct(int newSymbol, vector< int > product)
{
	IntegerType r;

	// product => newSymbol (this is a clause)
	// not x0 or not x1 or ... or not xn or newSymbol
	r = 1;
	beginConstraint();
	constraintTerm(1, newSymbol);
	for (int i = 0; i < product.size(); ++i)
		if (product[i] > 0)
		{
			constraintTerm(-1, product[i]);
			r -= 1;
		}
		else
			constraintTerm(1, -product[i]);
	constraintRelOp(">=");
	constraintRightTerm(r);
	endConstraint();

#ifdef ONLYCLAUSES
	// newSymbol => product translated as
	// not newSymbol of xi (for all i)
	for(int i=0;i<product.size();++i)
	{
	  r=0;
	  beginConstraint();
	  constraintTerm(-1,newSymbol);
	  if (product[i]>0)
	constraintTerm(1,product[i]);
	  else
	  {
	constraintTerm(-1,-product[i]);
	r-=1;
	  }
	  constraintRelOp(">=");
	  constraintRightTerm(r);
	  endConstraint();
	}
#else
	// newSymbol => product translated as
	// x0+x1+x3...+xn-n*newSymbol>=0
	r = 0;
	beginConstraint();
	constraintTerm(-(int)product.size(), newSymbol);
	for (int i = 0; i < product.size(); ++i)
		if (product[i] > 0)
			constraintTerm(1, product[i]);
		else
		{
			constraintTerm(-1, -product[i]);
			r -= 1;
		}
	constraintRelOp(">=");
	constraintRightTerm(r);
	endConstraint();
#endif
}
