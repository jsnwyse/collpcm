collpcm.run.monks.example <- function( run=1 )
{
	# runs the example from Ryan, Wyse and Friel
	
	#load the data
	data( monks )
	
	if( run == 1 )
	{
		control.mcmc <- list( sample = 10000, burn = 10000, interval = 10 )

		control.proposal <- list(variance.beta.pr = 0.5,variance.latent.position.pr = 0.7)

		control.initial <- list( beta= rnorm( 0, sd= 1) , ngroups= sample(2:5)[1] )

		#run the MCMC
		fit <- collpcm.fit( monks, control.initial=control.initial, control.proposal=control.proposal, control.mcmc=control.mcmc, update.gamma = TRUE )
		
	}else{
		
		control.mcmc <- list( sample = 10^4, burn = 10^5, interval = 10^2 )
		
		control.initial <- list( beta= rnorm( 0, sd= 1) , ngroups= sample(2:5)[1] )
		fit <- collpcm.fit( monks, control.initial=control.initial, control.mcmc=control.mcmc, update.gamma = TRUE, adapt = TRUE  )
	
	}

	#plot a summary of the main results of interest
	collpcm.summaryplot(fit)

	return( fit )
}

