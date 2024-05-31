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

	vector< Eigen::Triplet< IntegerType > > tmp_triplet;

	if (iRelOp == "=")
	{
		for (auto& t : tmp_trip)
		{
			tmp_triplet.emplace_back(iConsEq, get<0>(t), get<1>(t));
		}
		A_eq_trip.insert(A_eq_trip.end(), tmp_triplet.begin(), tmp_triplet.end());
		b_eq.emplace_back(iConsEq, get<1>(b.back()));
		iConsEq++;
	}
	else
	{
		assert(iRelOp == ">=");  // only support >= for now
		for (auto& t : tmp_trip)
		{
			tmp_triplet.emplace_back(iConsIneq, get<0>(t), get<1>(t));
		}
		A_ineq_trip.insert(A_ineq_trip.end(), tmp_triplet.begin(), tmp_triplet.end());
		b_ineq.emplace_back(iConsIneq, get<1>(b.back()));
		iConsIneq++;
	}

	for (auto& t : tmp_triplet)
	{
		cout << t.row() << " " << t.col() << " " << t.value() << endl;
	}
	cout << "RelOp: " << iRelOp << endl;
	cout << "b: " << get<1>(b.back()) << endl;

	tmp_trip.clear();
	iCons ++;

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

const std::vector< std::tuple< long, long, IntegerType > >& AbcCallback::getA()
{
	return A;
}

const std::vector< Eigen::Triplet < IntegerType > >& AbcCallback::getA_eq()
{
	return A_eq_trip;
}

const std::vector< Eigen::Triplet < IntegerType > >& AbcCallback::getA_ineq()
{
	return A_ineq_trip;
}


const std::vector< std::tuple< long, string > >& AbcCallback::getRelOp()
{
	return relOp;
}

const std::vector< std::tuple< long, IntegerType > >& AbcCallback::getB()
{
	return b;
}

const std::vector< std::tuple< long, IntegerType > >& AbcCallback::getB_eq()
{
	return b_eq;
}

const std::vector< std::tuple< long, IntegerType > >& AbcCallback::getB_ineq()
{
	return b_ineq;
}

const std::vector< std::tuple< long, IntegerType > >& AbcCallback::getC()
{
	return c;
}
