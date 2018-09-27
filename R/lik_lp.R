## conditional prior on sd/mean of species-specific normal likelihoods under LP model, 
# pars - vector of parameters
# x - vector of species means/sd
# tree - phylogenetic tree
# map - map of the jumps (of length equal to the number of edges)

lik_lp <- function(pars, x, tree, map){
  
  #@ parameters
  sigsq.bm <- pars[1]
  sigsq.jp <- pars[2]
  th.lp <- pars[3]
  # lam <- pars[4]
  Y <- as.matrix(x)
  
  #@ variance covariance matrix
  bm.vcv <- vcv(tree)*sigsq.bm
  lp.vcv <- njm(tree, map)*sigsa.jp
  vcv <- bm.vcv + lp.vcv
  
  #@ ancestral state vector
  n <- dim(vcv)[1]
  m      <- matrix(1, n, 1)
  m[,]  <- th.lp
  
  #@ calculations
  DET <- as.numeric(determinant(vcv, logarithm = T)$modulus)
  INV <- ginv(vcv)
  
  #@ prior
  log.lik <- try(-1/2 * (n * log(2 * pi) + DET + t(Y - m) %*% INV %*% (Y - m)), silent=T)
  
  if (is.na(log.lik) | (class(log.lik) == "try-error" )) {
    return(-Inf)
  } else {
    return(log.lik)
  }
  
}


## matrix that counts the number of jumps between species

njm <- function(tree, map, n){
  
  njm <- mrca(tree)
  anc <- unique(unlist(apply(njm, 1, unique)))
  nbj <- unlist(lapply(anc ,function(x){
    edges <- which(tree$edge[,1] %in% Ancestors(tree,x))
    return(sum(map[edges])/2)
  }))
  anc.loc <- lapply(anc, function(x) which(njm %in% x))
  for(i in 1:length(anc)){
    njm[anc.loc[[i]]] <- nbj[i]
  }
  return(njm)
}


