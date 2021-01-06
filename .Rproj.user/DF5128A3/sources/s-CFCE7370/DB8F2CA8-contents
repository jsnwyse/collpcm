collpcm.control.hyper <- function( x, nnodes )
{
	members <- c( "min.groups", "max.groups", "prior.groups", "mean.beta", "variance.beta", "gamma", "delta", "alpha", "kappa" )
	if( is.null(x) ){
		missing <- members
	} else {
		missing <- setdiff( members, names(x) )
	}
	for(i in missing)
	{
		if( i == "min.groups" )
			x$min.groups <- 1
		if( i == "max.groups" )
			x$max.groups <- floor(nnodes/2)
		if( i == "prior.groups" )
			x$prior.groups <- 1
		if( i == "mean.beta" )
			x$mean.beta <- 0
		if( i == "variance.beta" )
			x$variance.beta <- 2.
		if( i == "gamma" )	
			x$gamma <- 0.103
		if( i == "delta" )
			x$delta <- 2
		if( i == "alpha" )
			x$alpha <- 3
		if( i == "kappa" )
			x$kappa <- .1
	}	
		
	#check the entries to make sure ok
	if( x$min.groups >= x$max.groups )
		stop("The minimum number of groups must be less than the allowed maximum. Otherwise select a fixed number of groups.")
	if( x$variance.beta < 0 )
		stop("The prior variance for the abundance parameter must be positive.")
	if( x$prior.groups > 1 )
		stop("The prior for the number of groups can be either 0 (Richardson and Green) or 1 (Nobile and Fearnside)")
		#maybe this can be modified to let the user place their own prior on the number of groups?
	if( x$gamma < 0 )
		stop("The prior rate for the cluster precisions must be positive.")
	if( x$delta < 0 )
		stop("The prior shape for the cluster precisions must be positive.")
	if( x$alpha < 0 )
		stop("The symmetric Dirichlet prior parameter must be positive.")
	if( x$kappa < 0 )
		stop("The prior scaled precision (omega squared) must be positive.")	
	return( x )
}

