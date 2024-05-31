//
// Created by psw on 5/19/24.
//

#ifndef PBO_HEURISTICS__ABCCALLBACK_H_
#define PBO_HEURISTICS__ABCCALLBACK_H_

#include <vector>
#include <memory>
#include <string>
#include "DefaultCallback.h"
#include "../type_int.h"
#include "../Eigen/Sparse"

using namespace std;

class AbcCallback : public DefaultCallback
{
 private:
	vector< tuple< long, long, IntegerType > > A;
	vector< Eigen::Triplet< IntegerType > > A_eq_trip, A_ineq_trip;
	vector< tuple< long, string>> relOp;
	vector< tuple< long, IntegerType > > b;
	vector< tuple< long, IntegerType > > b_eq, b_ineq;
	vector< tuple< long, IntegerType>> c;
	long nVar;
	long nCons;
	long iCons;
	long iConsEq, iConsIneq;
	string iRelOp;
	vector< tuple< long, IntegerType > > tmp_trip;
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
	const std::vector< std::tuple< long, long, IntegerType > >& getA();
	const std::vector< Eigen::Triplet< IntegerType > >& getA_eq();
	const std::vector< Eigen::Triplet< IntegerType > >& getA_ineq();
	const std::vector< std::tuple< long, string > >& getRelOp();
	const std::vector< std::tuple< long, IntegerType > >& getB();
	const std::vector< std::tuple< long, IntegerType > >& getB_eq();
	const std::vector< std::tuple< long, IntegerType > >& getB_ineq();
	const std::vector< std::tuple< long, IntegerType > >& getC();
	long getNVar() const
	{
		return nVar;
	}
	long getNCons() const
	{
		return nCons;
	}
	long getNConsEq() const
	{
		return iConsEq;
	}
	long getNConsIneq() const
	{
		return iConsIneq;
	}

};

#endif //PBO_HEURISTICS__ABCCALLBACK_H_
