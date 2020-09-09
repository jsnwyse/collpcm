\name{collpcm.control.initial}
\alias{collpcm.control.initial}
\title{Set the initial values of the collapsed sampler.}
\description{
  \code{\link{collpcm.control.initial}} returns a list giving initial values for the parameters to be sampled using the collapsed sampler. 
}
\usage{
collpcm.control.initial( x, nnodes )
}
\arguments{
 \item{x}{A list which may already set some of the initial values. Those not given in the list are set to default values. The list can contain members \code{beta}, \code{latent.positions} and  \code{ngroups}. }
 \item{nnodes}{The number of nodes in the network.}
}

\value{

\code{\link{collpcm.control.initial}} returns a list containing initial values of the parameters to be sampled by the collapsed sampler with slots 
 \item{beta}{Initial value of the intercept parameter.}
 \item{latent.positions}{Initial values of the latent positions corresponding to each node.}
 \item{ngroups}{Initial number of clusters for the sampler.}
}
\seealso{ }

\references{
 
 Caitriona Ryan, Jason Wyse and Nial Friel. 
 \emph{Bayesian model selection for the latent position cluster model for Social Networks.}
 
}



