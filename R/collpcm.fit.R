collpcm.fit <- function( Y , d = 2, G=NULL, Gmax = NULL, control = list(), Xref = NA )
{
	#arguments to the function:
	#	Y:: a network object 
	#	directed:: is the network to be considered as directed?
	#	d:: the assumed dimension of the latent space (usually 2)
	#	control :: a list giving the settings

	#do checks first to make sure everything in order
  
  if( d < 1 || d != floor(d) ) stop("Argument d gives the dimension of the latent space and must be a positive integer")
  if( d > 2 ) warning("The package can visualize the latent space only when d=1 or d=2")
	
	A <- as.sociomatrix( Y ) 
	dir <- is.directed( Y )
	
	if( nrow(A) != ncol(A) )
		stop("Adjacency matrix must be a square matrix")

	n <- nrow( A )

	#check args and set defaults
	if( TRUE ){
		ncovariates <- 0	
		covariates <- NULL
	}else{
		ncovariates <- ncol(covariates)
	}
	
	control <- collpcm.control( control, n, d )
	
	if( !is.null(G) ) control$G <- G 
	if( !is.null(Gmax) ) control$Gmax <- Gmax
	
	hyparams <- c( control$xi, (control$psi)^2, control$gamma, control$delta, control$alpha, 1/control$kappa )
	
	prparams <- c( control$sd.beta.prop, control$sd.X.prop )
	
	#if( control$bradley.terry ) control$betainit <- 0 
	
	initparams <- c( control$betainit, 0. )
	
	len <- control$sample 
	
	if( !control$store.sparse )
	{
		ret.allocations <- numeric( len*n )
		ret.X <- numeric(len*n*d)
		ret.theta <- numeric(ncovariates*len)
		ret.abundance <- numeric(len)
		ret.ngroups <- numeric(len)
		ret.llike <- numeric(len)
		ret.gamma <- numeric(len)
		ret.kappa <- numeric(len)
	}else{
		ret.allocations <- NULL
		ret.X <- NULL
		ret.theta <- NULL
		ret.abundance <- NULL
		ret.llike <- NULL
		ret.gamma <- NULL
		ret.kappa <- NULL
	}

	ret.ngroups <- numeric(len)
	
	loo <- FALSE
	acc.rt.latent.positions <- 0
	acc.rt.beta <- 0
	acc.rt.metmoves <- numeric(3)
	acc.rt.ejab <- numeric(2)
	acc.rt.theta <- numeric(ncovariates)
	
	tstart.mcmc <- proc.time()
	
	mcmcout <- .C( "collapsed_lpcm", 						as.integer(A),
					 as.integer(n), 								as.integer(d), 
					 as.integer(ncovariates),					as.double(covariates), 
					 as.integer(dir), 							as.integer(control$Gmax), 
					 as.integer(control$G),						as.integer(control$sample), 
					 as.integer(control$burn), 				as.integer(control$interval), 
					 as.integer(control$model.search),	 	as.double(hyparams), 
					 as.double(control$Gprior),
					 prop.vars = as.double(prparams), 		as.double(initparams),
					 as.double(control$Xinit),		 			allocations = as.integer(ret.allocations), 
					 Xsamp = as.double(ret.X),				 	betasamp = as.double(ret.abundance), 
					 Gsamp = as.integer(ret.ngroups), 		llike = as.double(ret.llike), 
					 thetasamp = as.double(ret.theta), 		gammasamp = as.double(ret.gamma), 
					 kappasamp = as.double(ret.kappa), 		accrtlp = as.double(acc.rt.latent.positions), 
					 accrtab = as.double(acc.rt.beta),		accrtmet = as.double(acc.rt.metmoves), 
					 accrtejab = as.double(acc.rt.ejab), 	accrtth = as.double(acc.rt.theta), 
					 as.integer( 1 ), 					 		as.integer(control$gamma.update),
					 as.integer( FALSE ),
					 as.integer(control$pilot),				as.integer(control$store.sparse),
					 as.integer(control$adapt),				as.integer(control$adapt.interval),			
					 as.integer(FALSE),							as.integer(control$verbose),
					 PACKAGE = "collpcm" )
	
	tend.mcmc <- proc.time()
	
	zz <- list()
	
	#put the algorithm settings into results to make life easier later
	zz$call <- list( Y=Y, d = d, ncovariates =  ncovariates, control = control )
	
	#arrange the latent positions into a 3 indexed array with index 1 giving the iteration number...

	Gsamp <- mcmcout$Gsamp

	betasamp <- mcmcout$betasamp

	if( !control$store.sparse ) theta <- matrix( mcmcout$thetasamp, nrow = len, ncol = ncovariates, byrow=TRUE )
	
	uG <- sort( unique(Gsamp) )
	gdist <- numeric( length(uG) )
	j <- 1
	for( k in seq_along(uG) ) gdist[k] <- sum( Gsamp == uG[k] )
	gdist <- gdist/len
	Gdist <- matrix( nrow = length(gdist), ncol = 2 )
	colnames(Gdist) <- c( "Number groups", "Posterior probability" )
	rownames(Gdist) <- rep( " ", length(gdist) )
	Gdist[,1] <- uG
	Gdist[,2] <- gdist
	
	zz$sample <- list( G = Gsamp, beta = betasamp, llike = mcmcout$llike )
	zz$Gpost <- Gdist
	
	#do the label switching
	
	if( !control$store.sparse )
	{
	
		tstart.labels <- proc.time()
		
		if( control$verbose ) cat("\n\t Starting label switching post processing...")
		
		allocations <- matrix( mcmcout$allocations, nrow=len, ncol=n, byrow=TRUE )
		allocations <- allocations + 1
		labels <- collpcm.undo.label.switching( allocations, Gsamp )
		zz$sample$labels <- labels$relab
		zz$sample$Gslot <- labels$numcomponents
		zz$sample$label.probs <- labels$label.probs
		zz$sample$label.probs.idx <- labels$item.tags
		
		if( control$verbose ) cat("\n\t Finished label switching post processing...")
	
		tend.labels <- proc.time()
	
		#do the Procrustes
		
		if( control$verbose ) cat("\n\t Starting Procrustes matching...")
		
		tstart.procrustes <- proc.time()
		
		Xsamp <- array( mcmcout$Xsamp, c(d,n,len) )
		Xsamp <- aperm( Xsamp, c(2,1,3) )
		
		# centre all sampled positions around origin
		Xsamp <- apply( Xsamp, c(2,3), function(x) scale(x, scale=FALSE) )
		
		# reference positions
		idx <- which.max( mcmcout$llike )[1] # in case of ties
		if( is.na(Xref) ) Xref <- Xsamp[,,idx]
		zz$sample$Xref <- Xref
		
		# compute the crossprod
		crpr <- apply( Xsamp, c(2,3) , function(x) crossprod( Xref, x) )
		# get the svd
		if( d > 1 ) svald <- apply( crpr, 3, svd ) else svald <- apply( crpr, 2, svd)
		# get the rotation matrices
		rot <- lapply( svald, function(x) x$v %*% t(x$u) )
		# rotate the samples
		for( k in 1:len ) Xsamp[,,k] <- Xsamp[,,k] %*% rot[[k]]
		
		if( control$verbose ) cat("\n\t Finished Procrustes matching...")	
		tend.procrustes <- proc.time()
	
		numgr <- length(uG)
		
		nol <- paste0("G=",uG)
		zz$sample$X <- zz$Xpostmean <- vector( numgr, mode="list" )
		names(zz$sample$X) <- names(zz$Xpostmean) <- nol
		
		for( k in seq_along(uG) )
		{
		  idxs <- which( Gsamp == uG[k] )
		  zz$sample$X[[k]] <- Xsamp[,,idxs, drop=FALSE]
		  zz$Xpostmean[[k]] <- apply( zz$sample$X[[k]], c(1,2), mean ) 
		}
		
	}
	
	if( control$gamma.update ) zz$sample$gamma <- mcmcout$gammasamp
	
	if( !control$store.sparse )
	{
		if( control$verbose ) cat("\n\t Finding MKL positions...")
	
		#do the maximum Kullback-Liebler calc in here as need return intact for this
		if( control$MKL ) 
		{
		  zz$XpostMKL  <- vector( numgr, mode="list" )
		  names( zz$XpostMKL ) <- nol
		  
		  for( k in seq_along(uG) )
		  {
		    idxs <- which( Gsamp == uG[k])
		    zz$XpostMKL[[k]] <- collpcm.get.MKL.latent.positions( zz, betasamp[ idxs ] , zz$sample$X[[k]]   )[["X"]]
		  }
		  
		}
		
		if( control$verbose ) cat("\n\t MKL positions found...")
	}
	
	
	zz$acceptance.rates <- list( X = mcmcout$accrtlp, beta=mcmcout$accrtab, move1=mcmcout$accrtmet[1], move2=mcmcout$accrtmet[2], move3=mcmcout$accrtmet[3] )
	if( control$model.search ) zz$acceptance.rates <- c( zz$acceptance.rates, eject=mcmcout$accrtejab[1], absorb=mcmcout$accrtejab[2] )
	
	zz$adapted.sd.prop = list( beta = sqrt(mcmcout$prop.vars[1]), X = sqrt(mcmcout$prop.vars[2]) )
	
	
	times <- list()
	times$total <- system.time(1)
	times$mcmc <- tend.mcmc - tstart.mcmc
	times$total <- times$total + times$mcmc
	if( !control$store.sparse )
	{
		times$label.switching <- tend.labels - tstart.labels
		times$procrustes <- tend.procrustes - tstart.procrustes
		times$total <- times$total + times$label.switching + times$procrustes
	}
	zz$timings <- times

	if( control$verbose ) cat("\n\t Returning...\n")

	class( zz ) <- "collpcm"
	
	return( zz )
}

