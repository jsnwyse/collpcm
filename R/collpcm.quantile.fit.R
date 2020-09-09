collpcm.quantile.fit <- function( Y , d = 2, q = NULL, G=NULL, Gmax = NULL, control = list(), Xref = NA )
{
  #arguments to the function:
  #	Y:: a network object 
  #	directed:: is the network to be considered as directed?
  #	d:: the assumed dimension of the latent space (usually 2)
  # q:: the quantile value for the asymmetric Laplace distribution
  #	control :: a list giving the settings
  
  #do checks first to make sure everything in order
  
  A <- as.sociomatrix( Y ) 
  dir <- is.directed( Y )
  
  if( nrow(A) != ncol(A) )
    stop("Adjacency matrix must be a square matrix.")
  
  n <- nrow( A )
  
  # if no quantile passed determine from data 
  if( is.null(q) ) q <- sum(A) / ( (n-1)*n )
  
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
  
  mcmcout <- .C( "collapsed_quantile_lpcm", 		as.integer(A),
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
                 as.double(q),
                 PACKAGE = "collpcm" )
  
  tend.mcmc <- proc.time()
  
  zz <- list()
  
  #put the algorithm settings into results to make life easier later
  zz$call <- list( Y=Y, d = d, ncovariates =  ncovariates, control = control )
  
  #arrange the latent positions into a 3 indexed array with index 1 giving the iteration number...
  
  if( !control$store.sparse ) 
  {
    if( d > 1 ) 
    {
      Xsamp <- matrix( mcmcout$Xsamp, nrow = n*len, ncol = d, byrow=FALSE )
    }else{
      Xsamp <- mcmcout$Xsamp 
    }
    
  }
  
  Gsamp <- mcmcout$Gsamp
  
  betasamp <- mcmcout$betasamp
  
  if( !control$store.sparse ) theta <- matrix( mcmcout$thetasamp, nrow = len, ncol = ncovariates, byrow=TRUE )
  
  gdist <- numeric( max(Gsamp) - min(Gsamp) + 1 )
  j <- 1
  for( i in min(Gsamp):max(Gsamp) ){
    gdist[j] <- sum( Gsamp == i )
    j <- j + 1
  }
  gdist <- gdist/len
  Gdist <- matrix( nrow = length(gdist), ncol = 2 )
  colnames(Gdist) <- c( "Number groups", "Posterior probability" )
  rownames(Gdist) <- rep( " ", length(gdist) )
  Gdist[,1] <- min(Gsamp):max(Gsamp)
  Gdist[,2] <- gdist
  
  zz$sample <- list( G = Gsamp, beta = betasamp, llike = mcmcout$llike )
  zz$Gpost <- Gdist
  
  #do the label switching
  
  if( !control$store.sparse )
  {
    
    tstart.labels <- proc.time()
    
    allocations <- matrix( mcmcout$allocations, nrow=len, ncol=n, byrow=TRUE )
    allocations <- allocations + 1
    labels <- collpcm.undo.label.switching( allocations, Gsamp )
    zz$sample$labels <- labels$relab
    zz$sample$Gslot <- labels$numcomponents
    zz$sample$label.probs <- labels$label.probs
    
    
    #zz$labels <- labels
    
    #zz$labels.loo <- mcmcout$allocations
    
    #zz$loo <- loo
    
    
    tend.labels <- proc.time()
    
    #do the Procrustes
    
    tstart.procrustes <- proc.time()
    
    if( d > 1 )
    {
      t <- array( Xsamp, c( d, n, len ) )
      Xsamp <- aperm( t, c(2,1,3) )
    }else{
      Xsamp <- matrix( Xsamp, nrow=n, ncol=len  )
    }
    
    #latent.positions <- aperm( t, c(2,1,3) )
    
    idx <- which.max( mcmcout$llike )
    
    if( length(idx) > 1 ) idx <- idx[1]
    
    if( d > 1 )
    {
      if( is.na(Xref) ) Xref <- Xsamp[ , , idx]
    }else{
      if( is.na(Xref) ) Xref <- Xsamp[ , idx]
    }
    
    zz$sample$Xref <- Xref
    
    #match positions
    Xproc <- array( numeric(), c( n, d, len ) )
    
    if( control$verbose ) cat("\n\t Starting Procrustes matching...")
    
    if( d > 1 )
    {
      for(t in 1:len)
      {
        procr <- procrustes( Xref, Xsamp[,,t], truemean=FALSE, scale=FALSE )
        Xproc[,,t] <- procr$Yrot
      }
    }else{
      for( t in 1:len ) Xsamp[, t] <- Xsamp[,t] - mean( Xsamp[,t] - Xref )
    }
    
    if( control$verbose ) cat("\n\t Finished Procrustes matching...")	
    
    tend.procrustes <- proc.time()
    
    #special case for leave one out method
    
    #make a list with the latent positions for each number of groups
    latentpos <- list()
    Xpostmean <- list()
    
    j <- 1
    ug <- unique( Gsamp )
    for(i in ug ){
      idxs <- which( Gsamp == i )
      if( d > 1 ) 
      {
        latentpos[[j]] <- Xproc[,,idxs] 
        if( length(idxs) > 1 )
        {
          Xpostmean[[j]] <- apply( latentpos[[j]], c(1,2), mean  )
        }else{
          Xpostmean[[j]] <- latentpos[[j]]
        }
      }else{ 
        latentpos[[j]] <- Xsamp[,idxs]
        if( length(idxs) > 1 )
        {
          Xpostmean[[j]] <- apply( latentpos[[j]], 1 , mean  )
        }else{
          Xpostmean[[j]] <- latentpos[[j]]
        }
      }
      j <- j + 1
    }
    
    #sort the entries of latentpos and label
    
    ord <- sort( ug, index.return=TRUE )$ix
    
    zz$sample$X <- zz$Xpostmean <- list()
    
    for( k in 1:length(ord) ) 
    {
      zz$sample$X[[k]] <- latentpos[[ ord[k] ]]
      zz$Xpostmean[[k]] <- Xpostmean[[ ord[k] ]]
    }
    
    nol <- paste0("G=",sort(ug))
    
    names( zz$sample$X ) <- names( zz$Xpostmean ) <- nol
  }
  
  if( control$gamma.update ) zz$sample$gamma <- mcmcout$gammasamp
  
  
  if( !control$store.sparse )
  {
    if( control$verbose ) cat("\n\t Finding MKL positions...")
    
    #do the maximum Kullback-Liebler calc in here as need return intact for this
    if( control$MKL ) 
    {
      latentpos.mkl <- list()
      beta.mkl <- list()
      j <- 1
      ug <- unique( Gsamp )
      for(i in ug){
        idxs <- which( Gsamp == i)
        if( length(idxs) == 1 ) 
        {
          stre <- latentpos[[j]]
          if( d > 1 )
          {
            latentpos[[j]] <- array( dim=c( dim(stre) , 2 ) )
            latentpos[[j]][,,1] <- stre
          }else{
            latentpos[[j]] <- array( dim=c( n, 2 ) )
            latentpos[[j]][,1] <- stre
          }
        }
        if( length(idxs) > 1 ) 
        {
          latentpos.mkl[[j]] <- collpcm.get.MKL.latent.positions( zz, betasamp[ idxs ] , latentpos[[j]]   )[["X"]]
        }
        if( length(idxs) == 1 ) latentpos.mkl[[j]] <- stre
        j <- j + 1
      }
      
      zz$XpostMKL <- list()
      for( k in 1:length(ord) ) zz$XpostMKL[[ k ]] <- latentpos.mkl[[ ord[k] ]] 
      names( zz$XpostMKL ) <- nol
      
      zz$Gslot <- sort( ug )
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

