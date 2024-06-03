//
// Created by psw on 5/19/24.
//

#ifndef PBO_HEURISTICS__TYPE_INT_H_
#define PBO_HEURISTICS__TYPE_INT_H_

// define IntegerType
#ifdef useGMP
#include <gmpxx.h>

typedef mpz_class IntegerType;

#define PBOINTMAX 1L << 62 // FIXME: numeric_limits< mpz_class >::max() has bug, and 1LL will fail
#define PBOINTMIN -1L << 62

#else
//#warning this IntegerType may not be suitable for some input file with int size > 64. Consider using GMP
typedef long IntegerType;

#define PBOINTMAX 1L << 62
#define PBOINTMIN -1L << 62

#endif

#endif //PBO_HEURISTICS__TYPE_INT_H_
