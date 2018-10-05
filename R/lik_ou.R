# input: pars - c(alp,sig,the0,the1...theN), x - var or mean of trait by sp, tree and map
# does: calculate log-likelihood; see Butler and King 2004, appendix eq. A8 and A9 
lik_ou <- function(pars, x, tree, map){
  
  # extract variables
  alp  <- pars[1]
  sig <- pars[2]
  the  <- pars[3:length(pars)]
  t.vcv  <- vcv(tree)
  T.len  <- t.vcv[1, 1]
  n <- dim(t.vcv)[1]
  
  # calculate matricies
  w <- cbind(rep(exp(-alp * T.len), n), w.reg(tree, map, n, T.len, alp))
  E <- w%*%the
  V <- sig.sq/(2 * alpha) * (exp(-2 * alpha * (T.len-t.vcv)) * (1-exp(-2 * alpha * t.vcv)))
  
  log.lik.OU <- try((-n * log(2 * pi)/2 - (as.numeric(determinant(V)$modulus))/2 - (t(x - E)%*%ginv(V)%*%(x - E))/2),silent=T)
  
  #print(log.lik)
  if (is.na(log.lik.OU) | (class( log.lik.OU) == "try-error" )) {
    return(-Inf)
  } else {
    return(log.lik.OU)
  }
  
}

### default tuning for initial window size. The output can be manually modified using control.mcmc()
# input: x - sd (VOU) or mean (MOU) of each species trait value, nreg - number of regimes 
ws_ou <- function(x, nreg){
  
  ws <- list()
  ws$alp.ou <- 0.5
  ws$sig.ou <- 2
  ws$the.ou <- rep(sd(x), nreg+1) # 2 in the previous version?
  return(ws)
  
}
  

### default tuning for initial parameter values. The output can be manually modified using control.mcmc()
# input: x - sd (VOU) or mean (MOU) of each species trait value, nreg - number of regimes 
init_ou <- function(x, nreg){
  
  pv <- list()
  pv$alp.ou <- runif(1, 0.1, 1)
  pv$sig.ou <- runif(1, 0.5, 3)
  pv$the.ou <- rep(mean(x), nreg+1)
  return(pv)
  
}


# input: tree, map, n, T.len, alp
# does: calculates the difference exp(-alpha * t[gamma]) - exp(-alpha * t[gamma-1]) for each regime according to map and tree
w.reg <- function(tree, map, n, T.len, alp){
    
    #apply calc to each species  
    w.reg <- lapply(1:n, function(sp){
    
      nodage <- T.len - branching.times(tree)
      # finds all direct ancestors of sp i until the root
      root <- n + 1
      branch <- which(tree$edge[,2] == sp) # branch linking sp and its most recent ancestor
      anc <- tree$edge[branch,1] # most recent ancestor
      exp.time <- exp(alp * T.len) - exp(alp * nodage[as.character(anc)]) # see Butler and King 2004, appendix eq. A7 
      
      while(anc != root){
        desc <- anc
        next.branch <- which(tree$edge[,2] == desc)
        anc <- tree$edge[next.branch, 1]
        
        branch <- c(branch, next.branch)
        exp.time <- c(exp.time, exp(alp * nodage[as.character(desc)]) - exp(alp * nodage[as.character(anc)]))
      }
      
      ## calculates the time spent in each regime
      beta.sp <- map[branch,]/rowSums(map[branch,])
      
      ## calculates wik
      wik <- exp(-alp * T.len) * colSums(beta.sp*exp.time)
      
      
      return(wik)
      
    })
    
    w.reg <- do.call(rbind, w.reg)  
    return(w.reg)
  
}

