#' @import ape MASS


# input: pars - c(sig1, ..., sigN, the0), x - var or mean of trait by sp, tree and map
# does: calculate log-likelihood; 
lik_bm <-function(pars, x, tree, map, ...){
  
  #extract variables
	Y <- as.matrix(x)	
	sig <- pars[-length(pars)] # sigma(s)
	n  <- length(tree$tip.label)
	m  <- matrix(1, n, 1)
	m[,] <- pars[2] # ancestral mean
	
	#calculate matricies
	V <- v_reg(tree, map, n, sig)[tree$tip.label,tree$tip.label]
	
	log.lik.BM <- try((-n/2 * log(2 * pi) - (as.numeric(determinant(V)$modulus))/2 - 1/2 * (t(Y - m)%*%ginv(V)%*%(Y - m))), silent=T)
	
	if (is.na(log.lik.BM) | (class(log.lik.BM) == "try-error" )) {
			return(-Inf)
	} else {
		return(log.lik.BM)
	}

}


# input: tree, map, n, T.len, alp
# does: calculates the vcv matrix according to the different regimes
v_reg <- function(tree, map, n, sig){
  
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
