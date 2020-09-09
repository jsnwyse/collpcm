
#ifndef _REQUIRED_LIBS_H_
#define _REQUIRED_LIBS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

#include "atable.h"

// define the constants
#define TRUE 1
#define FALSE 0
#define log_2_pi 1.83787706640934533908
#define log_pi 1.1447298858494001639
// courtesty of BayesLogit pakage
#define __PI 3.141592653589793238462643383279502884197
#define HALFPISQ 0.5 * __PI * __PI
#define FOURPISQ 4 * __PI * __PI
#define __TRUNC 0.64
#define __TRUNC_RECIP 1.0 / __TRUNC


#endif
