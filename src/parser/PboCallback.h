//
// Created by psw on 5/19/24.
//

#ifndef PBO_HEURISTICS__ABCCALLBACK_H_
#define PBO_HEURISTICS__ABCCALLBACK_H_

#include <vector>
#include <memory>
#include <string>
#include "DefaultCallback.h"

using namespace std;

class PboCallback : public DefaultCallback
{
 private:
	vector< tuple< size_t, size_t, IntegerType > > A;
	vector< tuple< size_t, string > > relOp;
	vector< tuple< size_t, IntegerType > > b;
	vector< tuple< size_t, IntegerType > > b_eq, b_ineq;
	vector< tuple< size_t, IntegerType > > c;
	size_t nVar;
	size_t nCons;
	size_t iCons;
	size_t iConsEq, iConsIneq;
	string iRelOp;
	vector< tuple< size_t, IntegerType > > tmp_trip;
 public:
	void metaData(int nbvar, int nbconstr) override;
//    void beginObjective() override;
//    void endObjective() override;
	void objectiveTerm(IntegerType coeff, int idVar) override;
//    void beginConstraint() override;
	void endConstraint() override;
	void constraintTerm(IntegerType coeff, int idVar) override;
	void constraintRelOp(std::string relop) override;
	void constraintRightTerm(IntegerType val) override;
	void linearizeProduct(int newSymbol, std::vector< int > product) override;
	const std::vector< std::tuple< size_t, size_t, IntegerType > >& getA();
	const std::vector< std::tuple< size_t, string > >& getRelOp();
	const std::vector< std::tuple< size_t, IntegerType > >& getB();
	const std::vector< std::tuple< size_t, IntegerType > >& getB_eq();
	const std::vector< std::tuple< size_t, IntegerType > >& getB_ineq();
	const std::vector< std::tuple< size_t, IntegerType > >& getC();
	size_t getNVar() const
	{
		return nVar;
	}
	size_t getNCons() const
	{
		return nCons;
	}
	size_t getNConsEq() const
	{
		return iConsEq;
	}
	size_t getNConsIneq() const
	{
		return iConsIneq;
	}

};

#endif //PBO_HEURISTICS__ABCCALLBACK_H_
