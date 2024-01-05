#### The following functions are taken from latentnet package only modified for 
####  our specific set up and variable naming conventions (for clarity)
#### files in latentnet ergmm.probs.R, ergmm.statseval.R, ergmm.families.R, ergmm.latent.effects.R

#Based on 'statnet' project software (http://statnet.org).
#  For license and citation information see http://statnet.org
 
#Krivitsky P and Handcock M (2015). latentnet: Latent Position and Cluster Models 
		# for Statistical Networks. The Statnet Project (http://	www.statnet.org). R package version 2.7.1, http://CRAN.R-project.org/package=latentnet. 

collpcm.get.MKL.latent.positions <- function( nw, beta, latentpos )
{	
	#get the minimum Kullback-Liebler divergences as explained 
	# in Appendix A of Handcock et al (2007)
	
	#this usually gives a nicer representation of the network positions
	
	sample <- length( beta )
	
	#k <- dim( latentpos )
	n <- nw$call$Y$gal$n #k[1]
	d <- nw$call$d #k[2]
	
	# get the expected value of links
	ExY <- matrix( 0 , nrow=n, ncol=n )
	
	for( l in 1:sample )
	{
		Xl <- latentpos[,,l] #else Xl <- latentpos[,l]
 		eta <- beta[l] - as.matrix( dist(Xl) )
		ExY <- ExY + 1 / ( 1 + exp(-eta) )
	}
	ExY <- ExY/sample

	nw$EofY <- ExY
	
	#optim to optimize 
	
	optim.control <- list( fnscale=-1, maxit=1000 )
	
	start.val <- c( mean(beta), as.vector( nw$sample$Xref ) )
	
	optim.ret <- optim( par = start.val, fn = collpcm.get.llike, gr = collpcm.get.grad.llike, nw, method = "L-BFGS-B", lower=rep(-Inf,n*d+1),   control=optim.control )
	
	if( inherits( optim.ret, "try-error" ) ){
		warning("\n BFGS did not converge to find MKL positions\n")
		return( NULL )
	}
	
	#otherwise return the positions
	
	retlist <- list()
	retlist[[ "beta" ]] <- optim.ret$par[1]
	if( d > 1 ) 
		retlist[[ "X" ]] <- matrix( optim.ret$par[2:(n*d+1)], nrow=n, ncol=d )
	else
		retlist[[ "X" ]] <- optim.ret$par[2:(n+1)]
	
	#centre the MKL positions at the origin
	if( d > 1 )
	{
		m <- colMeans( retlist$X )
		retlist$X <- t( apply( retlist$X, 1, function(z) {z-m} ) ) 
	}else{
		retlist$X <- retlist$X - mean( retlist$X )
	}

	return( retlist )
}


  

