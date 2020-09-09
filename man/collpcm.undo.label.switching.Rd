\name{collpcm.undo.label.switching}
\alias{collpcm.undo.label.switching}
\title{Correct samples of label vectors for label switching.}
\description{
  \code{\link{collpcm.undo.label.switching}} is used to correct sampled label vectors for label switching using the method proposed by Nobile and Fearnside (2007) which relies on the assignment algorithm of Carpaneto and Toth (1980).
}
\usage{
collpcm.undo.label.switching( Z, Gsamp = NULL )
}
\arguments{
 \item{Z}{A matrix of dimensions (num samples) by n giving the sampled label vectors for each iteration of the MCMC run.}
 
 \item{Gsamp}{ A vector of length (num samples) giving the corresponding number of groups for each iteration of the MCMC run. }

}

\value{

\code{\link{collpcm.undo.label.switching}} returns a list with the following slots.
  
 \item{call}{The function call.}

 \item{relab}{The post processed Z matrix after label switching has been corrected for}
 
 \item{label.probs}{List giving the probability of belonging to each group for each item, each entry of the list corresponding to a given number of components.} 

 \item{permutation}{The permutation applied to each row to correct for label switching.}
 
}

\author{
Jason Wyse
}

\references{
Nobile, A. and Fearnside A. T. (2007). \emph{Bayesian finite mixtures with an unknown number of components: The allocation sampler}
 Statistics and Computing, Vol. 17, 147-162 <doi:10.1007/s11222-006-9014-7>
 
 Carpaneto, G. and Toth, P. (1980). \emph{Algorithm 548: Solution of the Assignment Problem [H]}
 ACM Transactions on Mathematical Software,  Vol.  6, 104-111 <doi:10.1145/355873.355883>
}


