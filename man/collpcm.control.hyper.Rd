\name{collpcm.control.hyper}
\alias{collpcm.control.hyper}
\title{Set the hyperparameters for the collapsed latent position cluster model}
\description{
  \code{\link{collpcm.control.hyper}} creates a list containing the hyperparameters for the collapsed latent position cluster model. Hyperparameters can be passed if wished and these will be appended to default values for those not specified. If none of the hyperparameters are specified, then a list of default values is created. 
}
\usage{
collpcm.control.hyper( x, nnodes )
}
\arguments{
 \item{x}{An optional list giving named hyperparameters and their 
 values with possible entries \code{min.groups}, \code{max.groups}, \code{prior.groups},	\code{mean.beta}, \code{variance.beta}, \code{gamma}, 
 \code{delta}, \code{alpha}, \code{kappa}.}
    \item{nnodes}{The number of nodes in the network.}
}

\value{

\code{\link{collpcm.control.hyper}} returns a list giving the hyperparameters of the collapsed model (see the elements of \code{x} above). 

 \item{min.groups}{The minimum allowed number of groups.}
 \item{max.groups}{The maximum allowed number of groups.}
 \item{prior.groups}{The type of prior to use on the number of groups.}
 \item{mean.beta}{The prior mean on the intercept parameter.}
 \item{variance.beta}{The prior variance on the intercept parameter.}
 \item{gamma}{Twice the rate of the Gamma prior on the cluster precision.}
 \item{delta}{Twice the shape of the Gamma prior on the cluster precision.}
 \item{alpha}{The parameter of the Dirichlet prior on group weights.}
 \item{kappa}{The scaling of the prior mean for the cluster centre (in units of cluster precision).}
}

\author{
Jason Wyse
}

\references{
 
 Caitriona Ryan, Jason Wyse and Nial Friel. 
 \emph{ Bayesian model selection for the latent position cluster model for Social Networks.}
 
}



