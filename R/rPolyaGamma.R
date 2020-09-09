rPolyaGamma <- function(n,par)
{
  x <- numeric(n)
  if( length(par) == 1 ) 
  {
    par <- rep(par,n) 
  }else{ 
    if( length(par) < n ){ 
      warning("length of par is less than n, only first entry will be used")
      par <- rep(par[1],n)
    }
  }
  stop <- 0
  w <- .C("simulate_random_polya_gamma", as.integer(n),
          as.double(par), x = as.double(x), stop = as.integer(stop),
          PACKAGE="collpcm") 
  if( w$stop > 0 ) stop("\n something went wrong when generating Polya-Gamma RVs")
  return(w$x)
}