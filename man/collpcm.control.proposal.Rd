\name{collpcm.control.proposal}
\alias{collpcm.control.proposal}
\title{Set the parameters of the proposal densities for the collapsed sampler.}
\description{
  \code{\link{collpcm.control.proposal}} sets the parameters of the proposal densities for the various proposals used in the collapsed sampler. 
}
\usage{
collpcm.control.proposal( x )
}
\arguments{
 \item{x}{A list which gives values for the parameters of the proposal densities. If empty, all values are given default values. This list can include members \code{variance.beta.pr}, \code{variance.latent.position.pr}. This is random walk Metropolis updating.}
    
}

\value{

\code{\link{collpcm.control.proposal}} returns a list with (possibly) default values of the parameters to be used in the proposal densities.
 \item{variance.beta.pr}{Proposal variance of the random walk proposal for the intercept.}
 \item{variance.latent.position.pr}{Proposal variance of the random walk proposal for the latent positions.}
}


\author{
Jason Wyse
}

\references{
 Caitriona Ryan, Jason Wyse and Nial Friel. 
 \emph{Bayesian model selection for the latent position cluster model for Social Networks.}
}




