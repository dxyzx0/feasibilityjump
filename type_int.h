//
// Created by psw on 5/19/24.
//

#ifndef PBO_HEURISTICS__TYPE_INT_H_
#define PBO_HEURISTICS__TYPE_INT_H_

// define IntegerType
#ifdef useGMP
#include <gmpxx.h>

typedef mpz_class IntegerType;

#define PBOINTMAX 2e64 // FIXME: numeric_limits< mpz_class >::max() has bug
#define PBOINTMIN -2e64
#else

#warning this IntegerType may not be suitable for some input file. Consider using GMP
typedef long IntegerType;

#define PBOINTMAX numeric_limits< IntegerType >::max()

#endif

#endif //PBO_HEURISTICS__TYPE_INT_H_
