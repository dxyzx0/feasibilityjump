//
// Created by psw on 5/19/24.
//

#include "AbcCallback.h"
#include <iostream>
#include <cassert>

using namespace std;

void AbcCallback::objectiveTerm(IntegerType coeff, int idVar)
{
	DefaultCallback::objectiveTerm(coeff, idVar);
	c.emplace_back(idVar - 1, coeff);
}

void AbcCallback::metaData(int nbvar, int nbconstr)
{
	DefaultCallback::metaData(nbvar, nbconstr);
	nVar = nbvar;
	nCons = nbconstr;
}

void AbcCallback::constraintTerm(IntegerType coeff, int idVar)
{
	DefaultCallback::constraintTerm(coeff, idVar);
	A.emplace_back(iCons, idVar - 1, coeff);  // idVar is 1-based, so we need to subtract 1
	tmp_trip.emplace_back(idVar - 1, coeff);
}

void AbcCallback::endConstraint()
{
	DefaultCallback::endConstraint();

	if (iRelOp == "=")
	{
		b_eq.emplace_back(iConsEq, get< 1 >(b.back()));
		iConsEq++;
	}
	else
	{
		assert(iRelOp == ">=");  // only support >= for now
		b_ineq.emplace_back(iConsIneq, get< 1 >(b.back()));
		iConsIneq++;
	}

//	cout << "RelOp: " << iRelOp << endl;
//	cout << "b: " << get< 1 >(b.back()) << endl;

	tmp_trip.clear();
	iCons++;

}

void AbcCallback::constraintRelOp(std::string relop)
{
	DefaultCallback::constraintRelOp(relop);
	relOp.emplace_back(iCons, relop);
	iRelOp = relop;
}

void AbcCallback::constraintRightTerm(IntegerType val)
{
	DefaultCallback::constraintRightTerm(val);
	b.emplace_back(iCons, val);
}

void AbcCallback::linearizeProduct(int newSymbol, std::vector< int > product)
{
	DefaultCallback::linearizeProduct(newSymbol, product);
	throw std::runtime_error("Not implemented");
}

const std::vector< std::tuple< size_t, size_t, IntegerType > >& AbcCallback::getA()
{
	return A;
}

const std::vector< std::tuple< size_t, string > >& AbcCallback::getRelOp()
{
	return relOp;
}

const std::vector< std::tuple< size_t, IntegerType > >& AbcCallback::getB()
{
	return b;
}

const std::vector< std::tuple< size_t, IntegerType > >& AbcCallback::getB_eq()
{
	return b_eq;
}

const std::vector< std::tuple< size_t, IntegerType > >& AbcCallback::getB_ineq()
{
	return b_ineq;
}

const std::vector< std::tuple< size_t, IntegerType > >& AbcCallback::getC()
{
	return c;
}
