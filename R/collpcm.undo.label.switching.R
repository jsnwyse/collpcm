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


collpcm.undo.label.switching <- function( Z, ngroups = NULL )
#undo label switching using Nobile and Fearnside approach
{

	# Z is a matrix of labels with row i indexing classifications
  	# to groups from 1 : ngroups[i] or 1:ngroups if ngroups is an integer 
  	
  	if( is.null(ngroups) ) stop("\t argument ngroups must be specified as either a fixed number of groups or a vector of the number of groups corresponding to each row of Z")
    
    if( length(ngroups) == 1 ) ngroups <- rep(ngroups, nrow(Z))
    
    #different no.'s groups
    
    Zrelab <- matrix( nrow = nrow(Z), ncol = ncol(Z) )
    
    Perm <- matrix( nrow = nrow(Z), ncol = max(ngroups) )
    
    nobs <- ncol(Z)
    
    G <- unique(ngroups)
    
    ret <- list()
    ret$groups <- ngroups
    ret$ncomponents <- numeric(length(G))
    ret$item.tags <- list()
    
    groups.not.done <- NULL
    
    j <- 1
    
    for(k in G)
    {
      
      idx.k <- which(ngroups == k)
      labels.k <- Z[idx.k,]
      
      nsamp.k <- length(idx.k)
      ngrp.k <- k
      
      permutation <- numeric(nsamp.k*k)
      
      if( length(idx.k) == 1 ) groups.not.done <- c( groups.not.done, k )

      if((k!=1) && (length(idx.k) > 1))
      {
      
        nonempty.k <- apply( labels.k, 1, function(x) length(unique(x)) )
        t.k <- sort( nonempty.k, index.return=TRUE )
        
        labels.arranged.k <- labels.k[ t.k$ix, ]
        
        item.tags.k <- idx.k[ t.k$ix ] #this is the actual index of each row in the original matrix Z
        
        labels.out.k <- numeric(nsamp.k*nobs)
        
        w <- .C(	"Relabel",								as.integer(nobs),
        				as.integer(nsamp.k),					as.integer(ngrp.k),
        				as.integer(labels.arranged.k),	x=as.integer(labels.out.k),
        			   xx = as.integer(permutation),		PACKAGE = "collpcm" )
        
        ret$ncomponents[j] <- k
        ret$memberships[[j]] <- matrix(w$x,nrow = nsamp.k,ncol=nobs,byrow=FALSE)
        ret$permutation[[j]] <- matrix(w$xx,nrow=nsamp.k,ncol=k,byrow=FALSE)
        
        #compute membership probabilities for each data pt
        
        probs <- matrix( nrow = nobs, ncol=k)
        for(id in 1:nobs){
        	for(c in 1:k){
        		probs[id,c] <- length(which(ret$memberships[[j]][,id] == c))
        	}
        }
        
        probs <- probs/nsamp.k
        
        ret$membership.probabilities[[j]] <- probs
 
        #for variable indicator
		ret$item.tags[[j]] <- item.tags.k  
		
		#store in the new Z matrix Zrelab
		
		Zrelab[ item.tags.k, ] <- ret$memberships[[j]]
		
		Perm[ item.tags.k, 1:k ] <- ret$permutation[[j]]
        
      }else{
        
        ret$ncomponents[j] <- k
        ret$memberships[[j]] <- labels.k
        ret$membership.probabilities[[j]] <- matrix( 0 , nrow=nobs, ncol=k )
        for( id in 1:nobs ) ret$membership.probabilities[[j]][ id, labels.k[id] ] <- 1
        idx.k <- which(ngroups == k)
        ret$item.tags[[j]] <- idx.k #only in this case
        
      }
      
      j = j+1
      
    }
    
    
    ord <- sort( ret$ncomponents, index.return=TRUE )$ix
    
    membprob <- list()
    
   for( k in 1:length(ord) ) membprob[[k]] <- ret$membership.probabilities[[ ord[k] ]]
    
    names( membprob ) <- paste0("G=", sort(ret$ncomponents) ) 
    
    x <- list()
    
    x$call <- match.call()
    
    x$relab <- Zrelab
    
    x$numcomponents <- sort( ret$ncomponents )
    
    x$components <- sort( ret$ncomponents )
    
    x$label.probs <- membprob
    
    x$permutation <- Perm
    
    return(x)
    
  
  
}

