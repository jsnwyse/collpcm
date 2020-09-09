\name{collpcm.control.mcmc}
\alias{collpcm.control.mcmc}
\title{Specify parameters determining the MCMC run}
\description{
  Specify the number of samples to be collected, burn in to be used, sub-sampling interval, whether variable model jumps are aloud, and whether to run a pilot sample in the initial model.
}
\usage{
 collpcm.control.mcmc( x )
}
\arguments{
 \item{x}{An optional list setting the parameters of the MCMC run which can contain elements  \code{sample}, \code{burn}, \code{interval}, \code{fixed.groups} . }
}

\value{

\code{\link{collpcm.control.mcmc}} returns a list giving the specifications of the MCMC run. 

 \item{sample}{The number of MCMC samples that should be collected and stored.}
 \item{burn}{The number of initial samples to be discarded as burn-in.}
 \item{interval}{The sub-sampling interval for storing iterations from the chain.}
 \item{fixed.groups}{Keep the number of groups fixed for entire run (logical).}

}
\seealso{ }

\author{
Jason Wyse
}
\references{
 
 Caitriona Ryan, Jason Wyse and Nial Friel. 
 \emph{ Bayesian model selection for the latent position cluster model for Social Networks.}
 
}



