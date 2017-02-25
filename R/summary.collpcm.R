summary.collpcm <- function( object, ... )
{

	x <- object
	#print a neat summary of the MCMC run
	cat("\nSummary of the collapsed LPCM run:\n")
	if( x$call$control$model.search ){
		cat("\nPosterior for the number of groups:\n")
		print(x$Gpost,type="html")
	}
	its <- x$call$control$burn + x$call$control$sample * x$call$control$interval
	cat("\nMCMC inferences based on : \n\t Samples: \t", x$call$control$sample,"\n\t Burn-in: \t", x$call$control$burn, "\n\t Interval: \t",x$call$control$interval,
	"\n\t \t\t ------- \n\t Iterations: \t", its )
	cat("\n\nAcceptance rates for the various moves (% accepted): \n\t Latent postions: \t",round(x$acceptance.rates$X,digits=2),"\n\t Beta (intercept): \t",round(x$acceptance.rates$beta,digits=2), "\n\t Move 1: \t\t",round(x$acceptance.rates$move1,digits=2),"\n\t Move 2: \t\t",round(x$acceptance.rates$move2,digits=2),"\n\t Move 3: \t\t",round(x$acceptance.rates$move3,digits=2))
	if( x$call$control$model.search ) cat("\n\t Eject move: \t\t",round(x$acceptance.rates$eject,digits=2),"\n\t Absorb move: \t\t",round(x$acceptance.rates$absorb,digits=2))
	#if(x$settings$ncovariates>0) cat("\n\t Covariates: ",round(x$acceptance.rates$covs,digits=2))
	
	cat("\n\n")
	
	cat("Run times for various parts of analysis in seconds:")
	cat("\n\tMCMC for samples : \t\t",round(x$timings$mcmc[3],digits=2),"\n\tPost-processing procedures:\n\tLabel switching : \t\t",round(x$timings$label.switching[3],digits=2),"\n\tProcrustes matching : \t\t",round(x$timings$procrustes[3],digits=2))
	cat("\n")
}


