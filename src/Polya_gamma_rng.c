
#include "Polya_gamma_rng.h"

/*these methods are translated directly from the C++ code in the BayesLogit R package*/

double Polya_gamma_a(int n, double x)
{
  double K = (n + 0.5) * __PI;
  double y = 0;
  if (x > __TRUNC) 
  {
    y = K * exp( -0.5 * K*K * x );
  }else if (x > 0) {
    double expnt = -1.5 * (log(0.5 * __PI)  + log(x)) + log(K) - 2.0 * (n+0.5)*(n+0.5) / x;
    y = exp(expnt);
  }
  return y;
}


double Polya_gamma_rtigauss(double z)
{
  z = fabs(z);
  double t = __TRUNC;
  double x = t + 1.0;
  if (__TRUNC_RECIP > z) 
  { // mu > t
    double alpha = 0.0;
    while(runif(0.0,1.0) > alpha) 
    {
      // Slightly faster to use truncated normal.
      double E1 = rexp(1.0); double E2 = rexp(1.0);
      while ( E1*E1 > 2.0 * E2 / t) 
      {
        E1 = rexp(1.0); E2 = rexp(1.0);
      }
      x = 1 + E1 * t;
      x = t / (x * x);
      alpha = exp(-0.5 * z*z * x);
    }
  }else{
    double mu = 1.0 / z;
    while(x > t) 
    {
      double y = rnorm(0., 1.); 
      y = y * y;
      double half_mu = 0.5 * mu;
      double mu_y = mu * y;
      x = mu + half_mu * mu_y - half_mu * sqrt(4 * mu_y + mu_y * mu_y);
      if( runif(0.0,1.0) > mu / (mu + x)) x = mu * mu / x;
    }
  }
  return x;
}

double Polya_gamma_mass_t_exp(double z)
{
  double t = __TRUNC;
  double fz = 0.125 * __PI*__PI + 0.5 * z*z;
  double b = sqrt(1.0 / t) * (t * z - 1);
  double a = sqrt(1.0 / t) * (t * z + 1) * -1.0;
  
  double x0 = log(fz) + fz * t;
  double xb = x0 - z + pnorm(b, 0.0, 1.0, TRUE, TRUE); 
  double xa = x0 + z + pnorm(a, 0.0, 1.0, TRUE, TRUE);
  
  double qdivp = 4.0 / __PI * ( exp(xb) + exp(xa) );
  
  return(1.0 / (1.0 + qdivp));  
}

double r_Polya_gamma_Devroye(double z)
{
  z = 0.5 * fabs(z);
  
  double fz = 0.125 * __PI * __PI + 0.5 * z*z;
  double x=0.0, s=1.0, y=0.0;
  int n, go;
  
  while(1)
  {
    if( runif(0.0,1.0) < Polya_gamma_mass_t_exp(z) ) x = __TRUNC + rexp(1.0) / fz ; else x = Polya_gamma_rtigauss(z) ;
    s = Polya_gamma_a(0,x);
    y = runif(0.0,1.0) * s;
    n = 0;
    go = TRUE;
    while( go && n < __MAX_IT_POLYA_GAMMA )
    {
      n += 1;
      if( n%2 == 1 )
      {
        //odd
        s -= Polya_gamma_a(n,x);
        if(y<=s) return(0.25*x);
      }else{
        //even
        s += Polya_gamma_a(n,x);
        if(y>s) go = FALSE;
      }
    }
    //can only break out here if max its reached
    return(-1.0);
  }
}

double r_Polya_gamma(int n, double z)
{
  if( n < 1 ) n = 1;
  int i;
  double s = 0.0, r;
  for( i=0; i<n; i++ )
  {
    r = r_Polya_gamma_Devroye(z);
    if( r < 0.0 ) return(r);
    s += r;
  }
  return(s);
}
