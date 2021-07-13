#=========================================================
#	Implementation of the label switching algorithm
#	of Nobile and Fearnside (2007) Statistics & Computing.
#
#	Author of this implementation:
#		Jason Wyse,
#		Discipline of Statistics,
#		School of Computer Science and Statistics,
#		Trinity College Dublin,
#		Dublin 2, Ireland.
#
#	Last update:
#		Fri 24 Feb 2017 21:42:08 GMT  
#=========================================================


collpcm.undo.label.switching <- function( Z, Gsamp = NULL )
  #undo label switching using Nobile and Fearnside approach
{
  
  # Z is a matrix of labels with row i indexing classifications
  # to groups from 1 : ngroups[i] or 1:ngroups if ngroups is an integer 
  
  if( is.null(Gsamp) ) stop("argument ngroups must be specified as either a fixed number of groups or a vector of the number of groups corresponding to each row of Z")
  
  if( length(Gsamp) == 1 ) ngroups <- rep( Gsamp, nrow(Z) )
  
  #different no.'s groups
  
  Zrelab <- matrix( nrow = nrow(Z), ncol = ncol(Z) )
  
  Perm <- matrix( nrow = nrow(Z), ncol = max(Gsamp) )
  
  nobs <- ncol(Z)
  
  G <- sort( unique(Gsamp) )
  
  ret <- list()
  ret$groups <- Gsamp
  #ret$ncomponents <- numeric(length(G))
  ret$item.tags <- vector( length(G), mode="list" ) #list()
  
  groups.not.done <- NULL
  
  for( j in seq_along(G) )
  {
    
    idx.k <- which( Gsamp == G[j] )
    labels.k <- Z[idx.k,]
    
    nsamp.k <- length(idx.k)
    ngrp.k <- G[j]
    
    permutation <- numeric( nsamp.k* G[j])
    
    if( length(idx.k) == 1 ) groups.not.done <- c( groups.not.done, G[j] )
    
    if( ( G[j] != 1 ) & ( length(idx.k) > 1 ) )
    {
      
      nonempty.k <- apply( labels.k, 1, function(x) length(unique(x)) )
      t.k <- sort( nonempty.k, index.return=TRUE )
      
      labels.arranged.k <- labels.k[ t.k$ix, ]
      
      item.tags.k <- idx.k[ t.k$ix ] #this is the actual index of each row in the original matrix Z
      
      labels.out.k <- numeric( nsamp.k * nobs )
      
      w <- .C(	"Relabel",								as.integer(nobs),
               as.integer(nsamp.k),					as.integer(ngrp.k),
               as.integer(labels.arranged.k),	x=as.integer(labels.out.k),
               xx = as.integer(permutation),		PACKAGE = "collpcm" )
      
      #ret$ncomponents[j] <- G[j]
      ret$memberships[[j]] <- matrix( w$x, nrow=nsamp.k, ncol=nobs, byrow=FALSE)
      ret$permutation[[j]] <- matrix( w$xx, nrow=nsamp.k, ncol=G[j], byrow=FALSE)
      
      #compute membership probabilities for each data pt
      
      probs <- matrix( nrow = nobs, ncol=G[j] )
      for(id in 1:nobs)
      {
        for(c in 1:G[j]) probs[id,c] <- length( which( ret$memberships[[j]][,id] == c ) )
      }
      
      probs <- probs/nsamp.k
      
      ret$membership.probabilities[[j]] <- probs
      
      #for variable indicator
      ret$item.tags[[j]] <- item.tags.k  
      
      #store in the new Z matrix Zrelab
      
      Zrelab[ item.tags.k, ] <- ret$memberships[[j]]
      
      Perm[ item.tags.k, 1:G[j] ] <- ret$permutation[[j]]
      
    }else{
      
      ret$ncomponents[j] <- G[j]
      ret$memberships[[j]] <- labels.k
      ret$membership.probabilities[[j]] <- matrix( 0 , nrow=nobs, ncol=G[j] )
      for( id in 1:nobs ) ret$membership.probabilities[[j]][ id, labels.k[id] ] <- 1
      idx.k <- which( Gsamp == G[j] )
      ret$item.tags[[j]] <- idx.k #only in this case
      
    }
    
  }
  
  # do a last step to match the labels overall with one another
  # start at smallest G and work up
  
  M <- ret$membership.probabilities[[1]]
  # do a permutation on the columns of this to put largest group first
  o <- order( colSums(M), decreasing=TRUE )
  ret$membership.probabilities[[1]] <- ret$membership.probabilities[[1]][,o] # reorder the columns
  for( k in ret$item.tags[[1]] ) Zrelab[k,] <- pmatch( Zrelab[k,], o, duplicates.ok = TRUE ) # relabel the memberships
  
  # now match decreasing groups to this
  for( j in 2:length(G) ) 
  {
    p <- permutations( G[j], G[j] )
    cst <- numeric( nrow(p) ) 
    prev.pr <- cbind( ret$membership.probabilities[[j-1]], matrix( 0, nrow=nobs, ncol= G[j]-G[j-1] ))
    for( k in 1:nrow(p) ) cst[k] <-  sum( prev.pr * ret$membership.probabilities[[j]][,p[k,]] )
    o <- p[ which.max(cst), ]
    ret$membership.probabilities[[j]] <- ret$membership.probabilities[[j]][,o]
    for( k in ret$item.tags[[j]] ) Zrelab[k,] <- pmatch( Zrelab[k,], o, duplicates.ok = TRUE )
  }
  
  names( ret$membership.probabilities ) <- paste0("G=", G ) 
  
  x <- list()
  
  x$call <- match.call()
  
  x$relab <- Zrelab
  
  x$numcomponents <- G
  
  x$components <- G
  
  x$label.probs <- ret$membership.probabilities
  
  x$permutation <- Perm
  
  return(x)
  
  
  
}

