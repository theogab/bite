## Utility functions to define initial conditions for the MCMC algorithm


### default tuning for initial window size. The output can be manually modified using control.mcmc()
###
init.ws <- function(x, model, nreg){
  
  ws <- list()
  
  # initial window size for changes in species mean and variance (Multivariate Normal)
  if(model == "MN"){
    ws$m.sp <- 2*x
    ws$v.sp <- 10*x
  }
  
  # initial window size for changes in Brownian Motion parameters
  if(model == "BM"){
    ws$sig.bm <- rep(2, nreg)
    ws$the.bm <- sd(x)
  }
  
  # initial window size for changes in Ornstein-Uhlenbeck parameters
  if(model == "OU"){
    ws$alp.ou <- 0.5
    ws$sig.ou <- 2
    ws$the.ou <- rep(sd(x), nreg+1) # 2 in the previous version?
  }
  
  # initial window size for changes in White Noise parameters
  if(model == "WN"){
    ws$sig.wn <- rep(0.5, nreg)
    ws$the.wn <- sd(x) # 2 in the previous version?
  }
  
  return(ws)
  
}


### default tuning for initial conditions. The output can be manually modified using control.mcmc()
###
init.pv <- function(x, model, nreg){
  
  pv <- list()
  
  # initial parameter value for Brownian Motion parameters
  if(model == "BM"){
    pv$sig.bm <- runif(nreg, 0.5, 3)
    pv$the.bm <- mean(x)
  }
  
  # initial parameter value for Ornstein-Uhlenbeck parameters
  if(model == "OU"){
    pv$alp.ou <- runif(1, 0.1, 1)
    pv$sig.ou <- runif(1, 0.5, 3)
    pv$the.ou <- rep(mean(x), nreg+1)
  }
  
  # initial parameter value for White Noise parameters
  if(model == "WN"){
    pv$sig.wn <- runif(nreg, 0.5, 3)
    pv$the.wn <- mean(x)
  }
  
  return(pv)
  
}