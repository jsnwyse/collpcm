#Plotting function has parts heavily based on code from latentnet package
# Attribution:
#Based on 'statnet' project software (http://statnet.org).
#  For license and citation information see http://statnet.org

plot.collpcm <- function( x, ..., G = NULL, label.nodes = NULL, pie = TRUE, vertex.col = c("red","green","blue","cyan","magenta",
                        "orange","yellow","purple"), vertex.cex = 1, object.scale = formals(plot.network.default)[["object.scale"]] )
{
	#input: a collpcm object

	if( class( x ) != "collpcm" ) stop("\n Argument is not of type collpcm" ) 
	
	if( x$call$d > 2 ) stop("Plotting only avilable for d <= 2")

	p <- x$Gpost
	if( is.null(G) ) G <- p[ which.max( p[  , 2] )  , 1 ]  
	
	idx <- which( x$sample$Gslot == G )
	if( length( idx ) == 0 )
		stop("There was no visit to ",G," group(s) in the sampler run: cannot plot")
	#if small posterior probability print a warning
	v <- min( x$sample$Gslot ): max( x$sample$Gslot )
	jx <- which( v == G )
	if( x$Gpost[ jx, 2 ] * x$call$control$sample < 100 )
		warning("The posterior probability for ",G," groups is small: this plot is based on less than 100 visits to this model")
	
	if( !x$call$control$MKL ) lpos <-  x$Xpostmean[[idx]] else lpos<- x$XpostMKL[[idx]] 
	
	nw <- x$call$Y #as.network( x$settings$network, directed = x$settings$directed )
	
	idx1 = which( x$sample$Gslot == G )
	
	if( G > 1 ) labels <- apply( x$sample$label.probs[[idx1]], 1, which.max ) else labels <- rep( 1, x$call$Y$gal$n )
	
	if( is.null(label.nodes) ) labs <- rep(" ", x$call$Y$gal$n ) else labs <- label.nodes
	
	if( x$call$d == 1 ) 
	{
		plot1D <- TRUE
		# idea taken from latentnet package
		normdist <- as.matrix( dist( lpos ) )
		normdist <- normdist / max( normdist )
		lpos <- cbind( lpos, rep(0,length(lpos) ) )
	}else{
		plot1D <- FALSE
		normdist <- NULL
	}
	
	#plotting ideas from latentnet package
	#Based on 'statnet' project software (http://statnet.org).
   #  For license and citation information see http://statnet.org
	pad <- .2
	xlims <- range( lpos[,1] )
	xlims <- xlims + c(-pad, pad)
	ylims <- range( lpos[,2] )
	ylims <- ylims + c(-pad,pad)
	if( plot1D ) ylims <- xlims
	
	plot.network( nw, coord = lpos, suppress.axes=FALSE, label= labs , arrowhead.cex=1.5, vertex.col= vertex.col[ labels ], edge.col = 8 , xlab="x1",  ylab="x2", xlim=xlims, ylim=ylims, object.scale=object.scale, usecurve = plot1D , edge.curve= normdist  )
	
	if(pie){
		#better way to get the pie radius 
		#plotting ideas from latentnet package
	#Based on 'statnet' project software (http://statnet.org).
   #  For license and citation information see http://statnet.org. 
		pie.rad <- collpcm.get.pie.radius( vertex.cex, xlims, ylims, object.scale )
		# limits
		#pie.radius <- min( diff(xlims), diff(ylims) ) * .05  
		probs = x$sample$label.probs[[idx1]]
		for(i in 1:x$call$Y$gal$n ){
			#pr <- probs[i,]
		  if( G > 1 ) pr <- probs[i,] else pr <- 1
			ergmm.drawpie( lpos[i,], radius=pie.rad, pr, n=50, cols=vertex.col[1:length(pr)] )
		}
	}

}

#Taken from latentnet package almost verbatim
	#Based on 'statnet' project software (http://statnet.org).
   #  For license and citation information see http://statnet.org
collpcm.get.pie.radius <- function( vertex.cex, xlims, ylims, object.scale )
{
	baserad<-min(diff(xlims),diff(ylims))*2.1*object.scale
	vertex.cex*baserad
}
