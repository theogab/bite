#' @import ape MASS

# input: pars - c(alp,sig,the0,the1...theN), x - var or mean of trait by sp, tree and map
# does: calculate log-likelihood; see Butler and King 2004, appendix eq. A8 and A9 
lik_ou <- function(pars, x, tree, map, root.station){
  
  # extract variables
  alp  <- pars[1]
  sig <- pars[2]
  the  <- pars[3:length(pars)]
  t.vcv  <- vcv(tree)
  T.len  <- t.vcv[1, 1]
  n <- dim(t.vcv)[1]
  
  # calculate matricies
  if(root.station){
    w <- w_reg(tree, map, n, T.len, alp)
  } else {
    w <- cbind(rep(exp(-alp * T.len), n), w_reg(tree, map, n, T.len, alp))
  }
  E <- w%*%the
  V <- sig/(2 * alp) * (exp(-2 * alp * (T.len-t.vcv)) * (1-exp(-2 * alp * t.vcv)))
  
  # log likelihood
  log.lik.OU <- try((-n * log(2 * pi)/2 - (as.numeric(determinant(V)$modulus))/2 - (t(x - E)%*%ginv(V)%*%(x - E))/2),silent=T)
  
  if (is.na(log.lik.OU) | (class( log.lik.OU) == "try-error" )) {
    return(-Inf)
  } else {
    return(log.lik.OU)
  }
  
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
