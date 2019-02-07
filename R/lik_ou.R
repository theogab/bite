#' @import ape

# input: n - number of species, n.p - which parameter has been updated, pars - c(alp,sig,the0,the1...theN), tree, map, vcv and logical for root.station 
# does: calculate log-likelihood; see Butler and King 2004, appendix eq. A8 and A9 
update_ou <- function(n, n.p, pars, tree, map, t.vcv, root.station){
  
  mat <- list(e = list(F),
              det = list(F),
              inv = list(F))
  
  ## extract variables
  alp  <- pars[2]/2*pars[1]
  sig <- pars[2]
  the  <- pars[3:length(pars)]
  T.len  <- t.vcv[1, 1]

  # alpha or theta(s) have been updated: change e
  if (any(c(1, 3:length(pars)) %in% n.p)) {
    if(root.station){
      w <- w_reg(tree, map, n, T.len, alp)
    } else {
      w <- cbind(rep(exp(-alp * T.len), n), w_reg(tree, map, n, T.len, alp))
    }  
    e <- w%*%the
    mat$e[[1]] <- T
    mat$e[[2]] <- e
  }
  
  # alpha or sigma have been updated: change v
  if (any(c(1,2) %in% n.p))
  {
    V <- sig/(2 * alp) * (exp(-2 * alp * (T.len-t.vcv)) * (1-exp(-2 * alp * t.vcv)))
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
# does: calculates the difference exp(-alpha * t[gamma]) - exp(-alpha * t[gamma-1]) for each regime according to map and tree
w_reg <- function(tree, map, n, T.len, alp){
  
  pp <- prop.part(tree)
  tree <- reorder(tree, "postorder")
  e1 <- tree$edge[, 1]
  e2 <- tree$edge[, 2]
  el <- map[paste(e1, e2, sep = ","),]
  if (dim(map)[2] == 1) el <- matrix(el) # a subset into a one column matrix becomes a vector (thanks R!)
  nodage <- T.len - branching.times(tree)
  w.reg <- matrix(0, dim(map)[2], n)
  
  for(i in length(e1):1){
    var.cur.node <- exp(alp * ifelse(e2[i] > n,  nodage[e2[i] - n], T.len)) - exp(alp * nodage[e1[i] - n])
    
    if(e2[i] > n){ # wether it is a species or a node
      desc <- pp[[e2[i] - n]]
    } else {
      desc <- e2[i]
    } 
    
    w.reg[,desc] <- w.reg[,desc] + ((el[i,]/sum(el[i,])) * var.cur.node)
  }
  
  w.reg <- exp(-alp * T.len) * t(w.reg) 
  return(w.reg)
  
}
