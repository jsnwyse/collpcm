collpcm.summaryplot <- function( x )
{
	
	#plot a summary of the results from a run!
	
	par(mfrow = c(2,2))
	
	par(mar = c(4.1,4.1,3.1,2.1))
	
	Gsamp <- x$sample$G
	
	plot( x$sample$G, type="l", col="blue", xlab="iterate", ylab="G", yaxt="n" )
	axis( side=2, at=min(Gsamp):max(Gsamp), labels=min(Gsamp):max(Gsamp) )
	
	plot( x$sample$llike, type="l", col="blue", xlab="iterate", ylab="log-likelihood" )
	
	plot( x$sample$beta, type="l", col="blue", xlab="iterate", ylab=expression(beta) )
	
	id <- which.max( x$Gpost[,2] )
	
	G <- x$Gpost[id,1]
	
	plot( x, G=G, pie=FALSE, vertex.cex=2 )
	
	par(mfrow=c(1,1))

}

