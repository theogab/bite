#' @import ape

# input: n - number of species, n.p - which parameter has been updated, pars - c(sig1, ..., sigN, the0), tree and map
# does: calculate log-likelihood; 
update_wnt <-function(n, n.p, pars, tree, map, ...){
  
  mat <- list(e = list(F),
              det = list(F),
              inv = list(F))
  sig <- pars[1] # sigma(s)
  
  ## calculate matricies
  
  # theta has been updated: change e
  if (any(length(pars) %in% n.p)) {
    e  <- matrix(1, n, 1)
    e <- pars[1 + as.numeric(as.factor(sapply(map[which(tree$edge[,2] <= n)], function(x) names(x)[length(x)])))] # ancestral mean
    mat$e[[1]] <- T
    mat$e[[2]] <- e
  }
  
  # sigma(s) have been updated: change v
  if (any(1:(length(pars)-1) %in% n.p))
  {
    V <- v_wnt(tree, map, n, sig)[tree$tip.label,tree$tip.label]
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
v_wnt <- function(tree, map, n, sig){
  
  if (is.null(tree$edge.length)){
    stop("the tree has no branch lengths")
  }
  
  pp <- prop.part(tree)
  tree <- reorder(tree, "postorder")
  e1 <- tree$edge[, 1]
  e2 <- tree$edge[, 2]
  nreg <- max(do.call(cbind,map)[1,])
  xx <- numeric(n + tree$Nnode)
  vcv <- matrix(0, n, n)
  
  # for each edge, calculate the variance accumulated from the root
  for(i in length(e1):1) { #loop ascending from the root to the first tip
    var.cur.node <- xx[e1[i]]
    reg.dur <- sapply(1:nreg, function(r){
      sum(map[[e2[i]]][3,map[[e2[i]]][1,] == r] - map[[e2[i]]][2,map[[e2[i]]][1,] == r])
    })
    xx[e2[i]] <- var.cur.node + sum(reg.dur * sig) # branch length under each regime * sig[regime]
  }
  
  # compute the diagonal
  diag.elts <- 1 + 0:(n - 1) * (n + 1)
  vcv[diag.elts] <- xx[1:n]
  dimnames(vcv)[1:2] <- list(tree$tip.label)
  
  return(vcv)
  
}
