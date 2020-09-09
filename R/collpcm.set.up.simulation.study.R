collpcm.set.up.simulation.study <- function( ez, beta = 1, N = 10^4 )
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
	
	# probability within cluster
	P.in <- numeric( m )
	# probability link between cluster 
	P.bw <- matrix( nrow=n , ncol=m )

	p <- numeric(N)

	for( s in c(-1,1) )
	{
		for( l in 1:m )
		{
			# simulate N realizations and compute 
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

	#image( Matrix(r) , colorkey=T )
	return(vals)
}

