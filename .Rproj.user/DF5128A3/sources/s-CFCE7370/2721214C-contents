
collpcm.control.initial <- function( x, nnodes )
{
	members <- c("beta","theta","latent.positions","ngroups")
	missing <- setdiff(members,names(x))
	for(i in missing)
	{
		if(i == "beta")
			x$beta <- -1.
		if(i == "theta")
			x$theta <- .001
		if(i == "latent.positions")
			x$latent.positions <- rnorm(2*nnodes,0,1)
		if(i == "ngroups")
			x$ngroups <- 2
	}
	#one check
	if( x$ngroups < 1 || floor(x$ngroups) != x$ngroups )
		stop("The initial number of groups must be positive integer.")
	#x$latent.positions[1:2] <- c(0,0)
	return( x )	
}

