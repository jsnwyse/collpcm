collpcm.run.dolphins.example <- function( )
{

	# runs the example from Ryan, Wyse and Friel

	#load the data 
	data( dolphins )
	
	control.mcmc <- list( sample = 10^4, burn = 10^5, interval = 10^2 )
	
	#control.mcmc = list(niterations=1100000, nburn=100000, thinby=100)

	control.proposal <- list(variance.beta.pr = 0.2, variance.latent.position.pr = 3)

	control.initial <- list( beta= rnorm( 0, sd= 1) , ngroups= sample(2:5)[1] )

	#run the MCMC
	fit <- collpcm.fit( dolphins,  control.initial=control.initial, control.proposal=control.proposal, control.mcmc=control.mcmc, update.gamma = TRUE )
	
	#plot a summary
	collpcm.summaryplot( fit )
	
	return( fit )
}

