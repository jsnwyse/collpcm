collpcm.control <- function( x = list() , n, d )
{
  
  vals <- list(
    # hyperparameters and priors
    G = sample(2:5, size=1), Gmax = floor(n/2), Gprior = NULL, xi = 0, psi = sqrt(2), 
    gamma = 0.103, delta = 2, alpha = 3, kappa = 0.1,
    
    # initialization
    betainit = rnorm( 0, sd=.01 ), Xinit = rnorm(d*n,0,1),
    
    # mcmc settings
    sample = 5000, burn = 5000, interval = 10, model.search = TRUE, pilot = 0,
    
    # proposals
    sd.beta.prop = sqrt(0.5), sd.X.prop = 1,
    
    # updates and storage
    gamma.update = TRUE, adapt = TRUE, adapt.interval = 200, store.sparse = FALSE, MKL = TRUE,
    
    # progress report
    verbose = FALSE
  ) 
  
  nmx <- names(x)
  m <- match( nmx, names(vals) ) 
  
  if( any(is.na(m)) ) stop( "Argument", nmx[ which(is.na(m)) ],"not recognized as a control argument." )
  
  for( k in seq_along(nmx) )
  {
    vals[[ nmx[k] ]] <- x[[ nmx[k] ]]
  }
  
  # set Gprior
  if( is.null( vals[["Gprior"]] ) ) vals[["Gprior"]] <- dpois( 0:vals[["Gmax"]], lambda=1., log=TRUE )
  
  return( vals )
}
