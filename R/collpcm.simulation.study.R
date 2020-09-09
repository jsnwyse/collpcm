collpcm.simulate.networks <- function( N, n, beta = 0, mu, tau, w, seed = NA )
{
	
	if( !is.na(seed) ) set.seed( seed )
	
	#simulate N networks	

	if( 1 ) 
	{
		Ymat <- array( dim=c( N, n, n) )	
		Xmat <- array( dim=c( N, n, 2) )
		labels <- matrix( nrow=N, ncol=n )
	}else{
	   Ymat <- matrix( n, n )	
	   labels  <- numeric(n)
	}
	
	X <- matrix( nrow=n, ncol= 2 )
	
	sig <- 1/sqrt(tau)

	for( k in 1:N )
	{
		#simulate the labels
		
		z <- sample( 1:2, size=n, replace=T, prob=w )
		i1 <- which( z==1 )
		i2 <- which( z==2 )
		
		X[ i1, ] <- matrix( rnorm( length(i1)*2, mean=mu, sd=sig ), ncol=2 )
		X[ i2, ] <- matrix( rnorm( length(i2)*2, mean=-mu, sd=sig ), ncol=2 )
	
		d <- as.matrix( dist( X ) ) 
		
		eta <- beta - d
		
		Pr <- ( 1 + exp( -eta ) )^{-1}
		
		p <- as.vector( Pr )
		y <- rbinom( n^2, size=1, prob=p )
		
		Y <- matrix( y, nrow=n, ncol=n )
		Y <- Y * lower.tri( Y )
		Y <- Y + t(Y)
		
		Ymat[ k, , ] <- Y
		Xmat[ k, , ] <- X
		labels[ k, ] <- z
		
	}

	return( list( N=N, n=n, Y=Ymat, X=Xmat, z=labels) )
}

collpcm.plot.simulated.network <- function( x, k=1, str=" ", vertex.cex = 1, object.scale = formals(plot.network.default)[["object.scale"]] )
{
	
	if( k > x$N ) stop("\nIndex provided is greater than no. of networks simulated")
	
	network <- as.network( x$Y[k,,], directed = FALSE )
	
	plot( network, coord = x$X[k,,], suppress.axes=FALSE, xlab="",  ylab="", edge.col="grey", main=str )
	
	
	pad <- .2
	xlims <- range( x$X[,1] )
	xlims <- xlims + c(-pad, pad)
	ylims <- range( x$X[,2] )
	ylims <- ylims + c(-pad,pad)
	
	probs = cbind( x$z[k,]==1, x$z[k,]==2 )
	pie.rad <- collpcm.get.pie.radius( vertex.cex, xlims, ylims, object.scale )
	for(i in 1:x$n){
		pr <- probs[i,]
		ergmm.drawpie( x$X[k,i,], radius=pie.rad,probs[i,],n=50,colours=c(2,3) )
	}
	
}


collpcm.set.up.simulation.study <- function( ez, beta = 0, N = 10^4 )
{
	
	#simple two cluster simulation study described in the paper
	# ez is a vector of "problem easiness" values (consider 1, 5, 10)
	# where 10 is easier than 1 etc

	#make a lookup table for the expected probability of a link
	#  within and between clusters

	mu <- seq( .1, 2, by=.1 )
	tau <- seq( 1, 20, by=1 )

	n <- length(mu)
	m <- length(tau)

	d <- 2

	# value of abundance
	b <- beta

	P.in <- numeric( m )
	P.bw <- matrix( nrow=n , ncol=m )

	p <- numeric(N)

	for( s in c(-1,1) )
	{
		for( l in 1:m )
		{
			#  simulate N realizations and compute 
			sdev = 1/sqrt(tau[l])
			if( s == 1 ) ner <- 1 else ner <- n
		
			for( k in 1:ner )
			{
				for( t in 1:N )
				{
					x1 <- rnorm( d, mean=mu[k], sd=sdev )
					x2 <- rnorm( d, mean= s*mu[k], sd=sdev )
					eta <- b - sqrt( sum( (x1-x2)^2 ) )
					p[t] <- 1/( 1 + exp( -eta ) )  
				}
				if( s == 1 ){
				 P.in[ l ] = mean(p) 
				}else{
				  P.bw[ k, l] = mean(p)
				}
			}
		}
	}

	r <-t( apply(P.bw, 1, function(x) P.in/x) )

	n.ez <- length( ez )
	
	vals <- list()
	
	for( k in 1:n.ez )
	{
		idx <- arrayInd( which.min( abs( r-ez[k] ) ) , dim(r)  )
		vals[[k]] <- list( ez = ez[k], actual.ez = r[ idx[1], idx[2] ], par = c(mu[ idx[1] ] , tau[ idx[2] ]) )
	}
	
	vals$P.in =  P.in
	vals$P.bw = P.bw

	#image( Matrix(r) , colorkey=T )
	return(vals)
}

