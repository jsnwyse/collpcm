print.collpcm <- function( x, ... )
{
	#print the MCMC run
	cat("\nSummary of the collapsed LPCM MCMC run:\n")
	if( x$call$control$model.search ){
		cat("\nPosterior for the number of groups:\n")
		xx <- x$Gpost
		print(x$Gpost,type="html")
	}else{
		cat("\nFitted a model with ", x$call$control$G , " groups\n")
	}
	
	cat("\n")

}


