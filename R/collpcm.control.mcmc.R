collpcm.control.mcmc <- function( x )
{
	#member names- check to see if present and if not set to default
	# check members for anomolies
	members <- c( "sample", "burn", "interval", "fixed.groups", "pilot" )
	if( is.null(x) ){
		missing <- members
	}else{
		missing <-  setdiff( members, names(x) )
	}
	for( i in missing )
	{
		if(i == "sample")
			x$sample <- 4000
		if(i == "burn")
			x$burn <- 2000
		if(i == "interval")
			x$interval <- 1
		if(i == "fixed.groups")
			x$fixed.groups <- FALSE
		#a pilot scheme without updating the no. of groups gives initialization a chance
		if( i== "pilot" )
			x$pilot <- 0
	}
	#if( x$niterations < x$nburn )
	#	stop("Please set both the number of iterations and burn in, or use default for both.")
	#if( x$thinby > 100  )
	#	warning("The MCMC thinning parameter seems quite large. Please check.") 
	return( x )
}
