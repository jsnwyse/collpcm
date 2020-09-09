#ifndef __POLYA_GAMMA_LIB__
#define __POLYA_GAMMA_LIB__

#include "required_libs.h"

#define __MAX_IT_POLYA_GAMMA 1000

double Polya_gamma_a(int n, double x);

double Polya_gamma_rtigauss(double z);

double Polya_gamma_mass_t_exp(double z);

double r_Polya_gamma_Devroye(double z);

double r_Polya_gamma(int n, double z);

#endif