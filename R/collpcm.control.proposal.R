collpcm.control.proposal <- function( x )
{
	members <- c( "variance.beta.pr", "variance.latent.position.pr", "variance.theta.pr" )
	missing <- setdiff( members, names(x) )
	for(i in missing)
	{
		if( i == "variance.beta.pr" )
			x$variance.beta.pr <- 0.5
		if( i == "variance.latent.position.pr" )
			x$variance.latent.position.pr <- 1
		if( i == "variance.theta.pr" )
			x$variance.theta.pr <- 0.1
	}
	#check all the parameters
	if( x$variance.beta.pr < 0 )
		stop("The proposal variance for the abundance parameter must be positive.")
	if( x$variance.latent.position.pr < 0 )
		stop("The proposal variance for the latent position update must be positive.")
	if( x$variance.theta.pr < 0 )
		stop("The proposal variance for the covariate update must be positive.")
	return( x )
}

