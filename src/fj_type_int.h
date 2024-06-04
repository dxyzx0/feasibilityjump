//
// Created by psw on 5/19/24.
//

#ifndef PBO_HEURISTICS__TYPE_INT_H_
#define PBO_HEURISTICS__TYPE_INT_H_

// define IntegerType
#ifdef useGMP
#include <gmpxx.h>

typedef mpz_class IntegerType;

const IntegerType PBOINTMAX = []{
	IntegerType result;
	mpz_ui_pow_ui(result.get_mpz_t(), 2, 1000);
	return result;
}();

const IntegerType PBOINTMIN = -PBOINTMAX;

#else
//#warning this IntegerType may not be suitable for some input file with int size > 64. Consider using GMP
typedef long IntegerType;

const IntegerType PBOINTMAX = (1L << 62);
const IntegerType PBOINTMIN = -PBOINTMAX;

#endif

#endif //PBO_HEURISTICS__TYPE_INT_H_
