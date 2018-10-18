# input: pars - c(sig1, ..., sigN, the0), x - var or mean of trait by sp, tree and map
# does: calculate log-likelihood; 
lik_wn <- function(pars, x, tree, map, ...){
	
  #extract variables
  Y <- as.matrix(x)	
  sig <- pars[-length(pars)] # sigma(s)
  n  <- length(tree$tip.label)
  m  <- matrix(1, n, 1)
  m[,] <- pars[2] # ancestral mean
  
  #calculate matricies
  V <- v_wn(tree, map, n, sig)[tree$tip.label,tree$tip.label]
  
  log.lik.WN <- try((-n/2 * log(2 * pi) - (as.numeric(determinant(V)$modulus))/2 - 1/2 * (t(Y - m)%*%ginv(V)%*%(Y - m))), silent=T)
  
  if (is.na(log.lik.WN) | (class(log.lik.WN) == "try-error" )) {
    return(-Inf)
  } else {
    return(log.lik.WN)
  }
  
}

### default tuning for initial window size. The output can be manually modified using control.mcmc()
# input: x - sd (VWN) or mean (MWN) of each species trait value, nreg - number of regimes 
ws_wn <- function(x, nreg){
  
  ws <- list()
  ws$sig.wn <- rep(0.5, nreg)
  ws$the.wn <- sd(x) # 2 in the previous version?
  return(ws)
  
}


### default tuning for initial parameter values. The output can be manually modified using control.mcmc()
# input: x - sd (VWN) or mean (MWN) of each species trait value, nreg - number of regimes 
init_wn <- function(x, nreg){
  
  pv <- list()
  pv$sig.wn <- runif(nreg, 0.5, 3)
  pv$the.wn <- mean(x)
  return(pv)
  
}



# input: tree, map, n, T.len, alp
# does: calculates the vcv matrix according to the different regimes
v_wn <- function(tree, map, n, sig){
  
  if (is.null(tree$edge.length)){
    stop("the tree has no branch lengths")
  }
  
  pp <- prop.part(tree)
  tree <- reorder(tree, "postorder")
  e1 <- tree$edge[, 1]
  e2 <- tree$edge[, 2]
  el <- map[paste(e1, e2, sep = ","),]
  xx <- numeric(n + tree$Nnode)
  vcv <- matrix(0, n, n)
  
  # for each edge, calculate the variance accumulated rom the root
  for(i in length(e1):1) { #loop ascending from the root to the first tip
    var.cur.node <- xx[e1[i]]
    xx[e2[i]] <- var.cur.node + sum(el[i,] * sig) # branch length under each regime * sig[regime]
  }
  
  # compute the diagonal
  diag.elts <- 1 + 0:(n - 1) * (n + 1)
  vcv[diag.elts] <- xx[1:n]
  dimnames(vcv)[1:2] <- list(phy$tip.label)
  
  return(vcv)
  
}
