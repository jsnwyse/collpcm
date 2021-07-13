\name{plot.collpcm}
\alias{plot.collpcm}
\title{Plotting a collpcm object}
\description{
  Plot the posterior mean latent positions for G groups.
}

\usage{
\method{plot}{collpcm}( x, ..., G = NULL, label.nodes = NULL, pie = TRUE, 

vertex.col = c( "red", "green", "blue", "cyan", "magenta", "orange", "yellow", "purple"), 

vertex.cex = 1, object.scale = formals(plot.network.default)[["object.scale"]] )
	
}
\arguments{
 \item{x}{An object of class \code{collpcm}.}
 \item{...}{Additional arguments including.}
 \item{G}{The number of groups in the model to be plotted. Defaults to most visited in MCMC run.}
 \item{label.nodes}{A vector of labels to print beside corresponding nodes on the plot.}
 \item{pie}{Logical; Draw small pie charts to indicate group membership probabilities.}
 \item{vertex.col}{The colour for the slices of pie (previous).}
 \item{vertex.cex}{Magnify the vertex}
 \item{object.scale}{Scale up/down the size of the plotting of vertex and arrows.}
}

\details{
This function gives a plot of the latent positions for a given number of groups (assuming the model with the specified number of groups has been visited during the run of the sampler). If argument \code{pie} is set to \code{TRUE}, membership probabilities of the nodes are indicated by pie charts with each colour corresponding to a different group in the model. Some of the code to implement this function draws heavily on code contained in the \code{latentnet} package (Krivitsky & Handcock, 2015).
}

\author{
Jason Wyse
}


\references{
 
 Ryan, C., Wyse, J. and Friel, N. (2017).
 \emph{Bayesian model selection for the latent position cluster model for Social Networks.} 
 Network Science, volume 5, 70-91.

Krivitsky P and Handcock M (2015). latentnet: Latent Position and Cluster Models for Statistical Networks. The Statnet Project (http://www.statnet.org). R package version 2.7.1, http://CRAN.R-project.org/package=latentnet. 
 
}



