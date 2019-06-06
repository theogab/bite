#' @import ape

# input: n - number of species, n.p - which parameter has been updated, pars - c(sig1, ..., sigN, the0), tree and map
# does: calculate log-likelihood; 
update_bm <-function(n, n.p, pars, tree, map, ...){
  
  mat <- list(e = list(F),
              det = list(F),
              inv = list(F))
  sig <- pars[-length(pars)] # sigma(s)

	## calculate matricies
	
	# theta has been updated: change e
	if (any(length(pars) %in% n.p)) {
	  e  <- matrix(1, n, 1)
	  e[,] <- pars[length(pars)] # ancestral mean
	  mat$e[[1]] <- T
	  mat$e[[2]] <- e
	}
	
	# sigma(s) have been updated: change v
	if (any(1:(length(pars)-1) %in% n.p))
	{
	  V <- v_bm(tree, map, n, sig)[tree$tip.label,tree$tip.label]
	  # determinant
	  det <- as.numeric(determinant(V)$modulus)
	  mat$det[[1]] <- T
	  mat$det[[2]] <- det
	  # inverse
	  inv <- solve(V)
	  mat$inv[[1]] <- T
	  mat$inv[[2]] <- inv
	}
	
	return(mat)

}


# input: tree, map, n, T.len, alp
# does: calculates the vcv matrix according to the different regimes
v_bm <- function(tree, map, n, sig){
  
  if (is.null(tree$edge.length)){
    stop("the tree has no branch lengths")
  }
  
  pp <- prop.part(tree)
  tree <- reorder(tree, "postorder")
  e1 <- tree$edge[, 1]
  e2 <- tree$edge[, 2]
  el <- map[paste(e1, e2, sep = ","),]
  if (dim(map)[2] == 1) el <- matrix(el) # a subset into a one column matrix becomes a vector (thanks R!)
  xx <- numeric(n + tree$Nnode)
  vcv <- matrix(0, n, n)
  
  # for each edge, calculate the variance accumulated rom the root
  for(i in length(e1):1) { #loop ascending from the root to the first tip
    var.cur.node <- xx[e1[i]]
    xx[e2[i]] <- var.cur.node + sum(el[i,] * sig) # branch length under each regime * sig[regime]
    
    # the covariance of the descendant to the right and to the left is the variance accumulated from the root to the node
    j <- i - 1L
    while(e1[j] == e1[i] && j > 0) {
      
      if(e2[j] > n){ # wether it is a species or a node
        left <- pp[[e2[j] - n]]
      } else {
        left <- e2[j]
      }
      
      if(e2[i] > n){ # wether it is a species or a node
        right <- pp[[e2[i] - n]]
      } else {
        right <- e2[i]
      }
      
      vcv[left, right] <- vcv[right, left] <- var.cur.node
      j <- j - 1L
    }
  }
  
  # compute the diagonal
  diag.elts <- 1 + 0:(n - 1) * (n + 1)
  vcv[diag.elts] <- xx[1:n]
  dimnames(vcv)[1:2] <- list(tree$tip.label)
  
  return(vcv)
  
}
