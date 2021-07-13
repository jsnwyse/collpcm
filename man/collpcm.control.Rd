\name{collpcm.control}
\alias{collpcm.control}
\title{Specify parameters determining the collapsed LPCM model and MCMC fitting run}
\description{
  Specify the number of samples to be collected, burn in to be used, sub-sampling interval, whether variable model jumps are allowed, and whether to run a pilot sample in the initial model.
}
\usage{
 collpcm.control( x = list() , n, d )
}
\arguments{
 \item{x}{An optional list setting the set up parameters of the model. Any parameters not set in the list will default to the values described below. }
 \item{n}{The number of nodes in the network.}
 \item{d}{The dimension of the latent space for model fitting.}
}

\value{

\code{\link{collpcm.control}} returns a list giving the set up of the problem containing the following items:

 \item{G}{Initial value of \code{G} for the chain.}
 \item{Gmax}{The maximum allowed value of \code{G} if doing model search.}
 \item{Gprior}{Log of the prior mass on the number of components \code{G}.}
 \item{xi}{Mean of the prior on the model intercept.}
 \item{psi}{Standard deviation of the prior on the model intercept.}
 \item{gamma}{Twice the rate of the Gamma prior on the cluster precision.}
 \item{delta}{Twice the shape of the Gamma prior on the cluster precision.}
 \item{alpha}{The parameter of the Dirichlet prior on group weights.}
 \item{kappa}{The scaling of the prior mean for the cluster centre (in units of cluster precision).}
 \item{betainit}{Initial value given to the intercept for the MCMC run.}
 \item{Xinit}{Initial configuration of latent positions for the MCMC run.}
 \item{sample}{Number of MCMC samples to be stored.}
 \item{burn}{Number of MCMC iterations to discard as burn-in.}
 \item{interval}{Number of iterations at which to sub-sample the chain and store i.e. total iterations post burn-in is sample*interval.}
 \item{model.search}{Logical; if TRUE (default) the model space for G is searched.}
 \item{pilot}{Number of iterations to run as a pilot to adapt the proposal standard deviations for the MCMC chains (in addition to adaptation during burn-in).}
 \item{sd.beta.prop}{Standard deviation of the random walk proposal updating the intercept.}
 \item{sd.X.prop}{Standard deviation of the (possibly multivariate) random walk proposal for an actor's latent position.}
 \item{gamma.update}{Logical; if TRUE (default) then the \code{gamma} hyperparameter is updated as part of the MCMC run.}
 \item{store.sparse}{Logical; do a sparse form of storage and don't return or store some of the MCMC run and only keep summary values.}
 \item{adapt}{Logical; if TRUE (default) use an adaptive phase during burn-in to tune the standard deviation of the proposals to get an "optimal" acceptance rate.}
 \item{adapt.interval}{The number of iterations between tweaks of the proposal standard deviations in the adaptation phase.}
 \item{MKL}{Logical; if TRUE (default) compute the maximum Kullback-Liebler configuration of the latent positons from Handcock, Raftery & Tantrum (2007) }
 \item{verbose}{Logical; if TRUE (default) print out progression messages througout the MCMC run and stages of fitting.}

}

\author{
Jason Wyse
}

\references{
Ryan, C., Wyse, J. and Friel, N. (2017)
 \emph{Bayesian model selection for the latent position cluster model for Social Networks.} 
 Network Science, volume 5, 70-91.
}



