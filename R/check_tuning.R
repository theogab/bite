## check validity of jive tuning
# input : jive object

check_tuning <- function(jive){
  
  ### Likelihood level ###
  # window sizes
  if(length(jive$lik$ws$m.sp) != length(jive$data$traits)) stop("Vector of window sizes for mean values is not the same length as the number of species")
  if(length(jive$lik$ws$v.sp) != length(jive$data$traits)) stop("Vector of window sizes for variance values is not the same length as the number of species")
  if(any(is.na(c(jive$lik$ws$m.sp, jive$lik$ws$v.sp)))) stop("Missing values are not allowed in window sizes: check obj$lik$ws")
  # initial values
  if(length(jive$lik$init$m.sp) != length(jive$data$traits))  stop("Vector of initial values for means is not the same length as the number of species")
  if(length(jive$lik$init$v.sp) != length(jive$data$traits))  stop("Vector of initial values for variances is not the same length as the number of species")
  if(any(is.na(unlist(jive$lik$init))))  stop("Missing values are not allowed in initial values: check obj$lik$init")
  
  ## Prior level ##
  for(level in c("prior.mean", "prior.var")){
    
    if(any(is.na(unlist(jive[[level]]$ws)))) stop(sprintf("Missing values are not allowed in window sizes: check obj$%s$ws", level))
    if(any(is.na(unlist(jive[[level]]$init)))) stop(sprintf("Missing values are not allowed in initial values: check obj$%s$init", level))
    if(any(unlist(jive[[level]]$ws) <= 0)) stop(sprintf("window sizes should not be negative: check obj$%s$ws", level))
    
    # White Noise #
    if (grepl("WN", jive[[level]]$name)){
      nreg <- max(do.call(cbind,jive[[level]]$map)[1,])
      # window sizes
      if(length(unlist(jive[[level]]$ws)) != (nreg + 1)) stop(sprintf("Vector of window sizes for %s is not the same length as the number of parameters: check obj$%s$ws", level))
      # initial values
      if(any(jive[[level]]$init$wn.sig <= 0)) stop(sprintf("initial values for sigma^2 should not be negative: check obj$%s$init", level))
      if(length(unlist(jive[[level]]$init)) != (nreg + 1)) stop(sprintf("Vector of initial values for %s is not the same length as the number of parameters: check obj$%s$init", level))
      for(i in 1:nreg){
        if(is.finite(jive[[level]]$hprior[[i]](-1)[[1]])){
          stop(sprintf("Hyper prior should not allow sigma^2 <= 0: check obj$%s$hprior$%s", level, names(jive[[level]]$hprior[i])))
        }  
      }
    } 
    
    # Brownian Motion
    if (grepl("BM", jive[[level]]$name)){
      nreg <- max(do.call(cbind,jive[[level]]$map)[1,])
      # window sizes
      if(length(unlist(jive[[level]]$ws)) != (nreg + 1)) stop(sprintf("Vector of window sizes for %s is not the same length as the number of parameters: check obj$%s$ws", level))
      # initial values
      if(any(jive[[level]]$init$bm.sig <= 0)) stop(sprintf("initial values for %s should not be negative: check obj$%s$init", level))
      if(length(unlist(jive[[level]]$init)) != (nreg + 1)) stop(sprintf("Vector of initial values for %s is not the same length as the number of parameters: check obj$%s$init", level))
      # hyperprior
      for(i in 1:nreg){
        if(is.finite(jive[[level]]$hprior[[i]](-1)[[1]])){
          stop(sprintf("Hyper prior should not allow sigma^2 <= 0: check obj$%s$hprior$%s", level, names(jive[[level]]$hprior[i])))
        }  
      }
    } 
    
    # Ornstein-Uhlenbeck #
    if (grepl("OU", jive[[level]]$name)){
      # number of regimes for each parameter
      nreg <- max(do.call(cbind,jive[[level]]$map)[1,])
      rsig <- ifelse(grepl("sigma",  jive[[level]]$name), nreg, 1)
      rsv <- ifelse(any(grepl("alpha", jive[[level]]$name), grepl("sigma", jive[[level]]$name)), nreg, 1)
      rthe <- ifelse(grepl("theta", jive[[level]]$name), nreg, 1) + ifelse(grepl("root", jive[[level]]$name), 1, 0)
      # window sizes
      if(length(unlist(jive[[level]]$ws)) != (rsig + rsv + rthe)) stop(sprintf("Vector of window sizes for %s is not the same length as the number of parameters: check obj$%s$ws", level))
      # initial values
      if(any(jive[[level]]$init$ou.sig <= 0)) stop(sprintf("initial values for %s should not be negative: check obj$%s$init", level))
      if(length(unlist(jive[[level]]$init)) != (rsig + rsv + rthe)) stop(sprintf("Vector of initial values for %s is not the same length as the number of parameters: check obj$%s$init", level))
      # hyperprior
      for(i in rsv + 1:rsig){
        if(is.finite(jive[[level]]$hprior[[2]](-1)[[1]])){
          stop(sprintf("Hyper prior should not allow sigma^2 <= 0: check obj$%s$hprior$%s", level, names(jive[[level]]$hprior[i])))
        } 
      }
    }
  }
}