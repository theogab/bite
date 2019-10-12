#' @import ape

# input: n - number of species, n.p - which parameter has been updated, pars - c(alp,sig,the0,the1...theN), tree, map, vcv and vector of size 3 giving the number of regime every parameter is undergoing for nr 
# does: calculate log-likelihood; see Butler and King 2004, appendix eq. A8 and A9 and Beaulieu et al. 2012
update_ou <- function(n, n.p, pars, tree, map, t.vcv, nr){
  
  mat <- list(e = list(F),
              det = list(F),
              inv = list(F))
  ## extract variables
  sv_i <- 1:nr[1]
  sig_i <- nr[1] + 1:nr[2]
  the_i <- nr[1] + nr[2] + 1:nr[3]
  alp  <- pars[sig_i]/(2*pars[sv_i])
  sig <- pars[sig_i]
  the  <- pars[the_i]
  T.len  <- t.vcv[1, 1]

  # update expectation or variance matrices depending on the updated parameter
  pp <- prop.part(tree)
  tree <- reorder(tree, "postorder")
  e1 <- tree$edge[, 1]
  e2 <- tree$edge[, 2]
  nreg <- max(do.call(cbind,map)[1,])
  if(length(sig) == 1) sig <- rep(sig, nreg)
  if(length(alp) == 1) alp <- rep(alp , nreg)
  
  ## Matrices used for the calculations of e and v
  if(any(c(sv_i, the_i) %in% n.p)){
    var.cur.reg <- matrix(0,nreg,n)
  }
  
  var.reg <- matrix(0,nreg,n)
  
  if(any(c(sig_i, sv_i) %in% n.p)){
    covar.reg <- matrix(0,n,n)
  }
  
  # loop over all edges
  for(i in length(e1):1){
    
    var.node <- sapply(1:nreg, function(r){
      -alp[r] * sum(map[[e2[i]]][3,map[[e2[i]]][1,] == r] - map[[e2[i]]][2,map[[e2[i]]][1,] == r])
    })
    
    if(any(c(sig_i, sv_i) %in% n.p)){
      covar.node <- sum(sapply(1:nreg, function(r){
        sig[r]/(2*alp[r]) * sum(exp(2 * alp[r] * map[[e2[i]]][3,map[[e2[i]]][1,] == r]) - exp(2 * alp[r] * map[[e2[i]]][2,map[[e2[i]]][1,] == r]))
      }))
    }
    
    if(any(c(the_i, sv_i) %in% n.p)){
      var.cur.node <- sapply(1:nreg, function(r){
        sum(exp(alp[r] * map[[e2[i]]][3,map[[e2[i]]][1,] == r]) - exp(alp[r] * map[[e2[i]]][2,map[[e2[i]]][1,] == r]))
      }) 
    }
    
    if(e2[i] > n){ # wether it is a species or a node
      desc <- pp[[e2[i] - n]]
    } else {
      desc <- e2[i]
    } 
    
    # increment each matrix with the information of the ith edge
    var.reg[,desc] <- var.reg[,desc] + var.node
    
    if(any(c(sig_i, sv_i) %in% n.p)){
      covar.reg[desc, desc] <- covar.reg[desc, desc] + covar.node
    }
    
    if(any(c(the_i ,sv_i) %in% n.p)){
      var.cur.reg[,desc] <- var.cur.reg[,desc] + var.cur.node
    }
  }
  
  
  if(any(c(sv_i, the_i) %in% n.p)){
    # calculate weight matrix
    if(any(is.infinite(var.cur.reg))){stop("inf values", call. = F)}
    w <- exp(t(var.reg))*t(var.cur.reg)
    if(nr[3] > nreg) w <- cbind(exp(colSums(var.reg)), w)
    w <- w/rowSums(w)
  
    # calculate expectation
    e <- w%*%the
    mat$e[[1]] <- T
    mat$e[[2]] <- e
  }
  
  # calculate variance covariance
  if (any(c(sv_i, sig_i) %in% n.p))
  {
    if(any(is.infinite(covar.reg))){stop("inf values", call. = F)}
    v <- exp(t(matrix(0, n, n) + colSums(var.reg)) + colSums(var.reg)) * covar.reg
    # determinant
    det <- as.numeric(determinant(v)$modulus)
    mat$det[[1]] <- T
    mat$det[[2]] <- det
    # inverse
    inv <- solve(v)
    mat$inv[[1]] <- T
    mat$inv[[2]] <- inv
  }
  
  return(mat)

}
