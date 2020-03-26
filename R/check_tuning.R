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
    
    temp <- jive[[level]]
    
    if(any(is.na(temp$ws))) stop(sprintf("Missing values are not allowed in window sizes: check obj$%s$ws", level))
    if(any(is.na(temp$init))) stop(sprintf("Missing values are not allowed in initial values: check obj$%s$init", level))
    if(any(temp$ws <= 0)) stop(sprintf("window sizes should not be negative: check obj$%s$ws", level))
    
    # White Noise #
    if (grepl("WN", temp$name)){
      nreg <- ncol(temp$map$beta)
      # window sizes
      if(length(temp$ws) != (nreg + 1)) stop(sprintf("Vector of window sizes for %s is not the same length as the number of parameters: check obj$%s$ws", level))
      # initial values
      if(any(temp$init[temp$Pi[1,]==1] <= 0)) stop(sprintf("initial values for sigma^2 should not be negative: check obj$%s$init", level))
      if(length(temp$init) != (nreg + 1)) stop(sprintf("Vector of initial values for %s is not the same length as the number of parameters: check obj$%s$init", level))
      for(i in 1:nreg){
        if(is.finite(temp$hprior[[i]](-1)[[1]])){
          stop(sprintf("Hyper prior should not allow sigma^2 <= 0: check obj$%s$hprior$%s", level, names(temp$hprior[i])))
        }  
      }
    } 
    
    # Brownian Motion
    if (grepl("BM", temp$name)){
      nreg <- ncol(temp$map$beta)
      # window sizes
      if(length(unlist(temp$ws)) != (nreg + 1)) stop(sprintf("Vector of window sizes for %s is not the same length as the number of parameters: check obj$%s$ws", level))
      # initial values
      if(any(temp$init[temp$Pi[1,]==1] <= 0)) stop(sprintf("initial values for %s should not be negative: check obj$%s$init", level))
      if(length(temp$init) != (nreg + 1)) stop(sprintf("Vector of initial values for %s is not the same length as the number of parameters: check obj$%s$init", level))
      # hyperprior
      for(i in 1:nreg){
        if(is.finite(temp$hprior[[i]](-1)[[1]])){
          stop(sprintf("Hyper prior should not allow sigma^2 <= 0: check obj$%s$hprior$%s", level, names(temp$hprior[i])))
        }  
      }
    } 
    
    # Ornstein-Uhlenbeck #
    if (grepl("OU", temp$name)){
      # number of regimes for each parameter
      nreg <- ncol(temp$map$beta)
      ralp <- sum(temp$Pi[1,])
      rsig <- sum(temp$Pi[2,])
      rthe <- sum(temp$Pi[3,])
      if(nrow(temp$Pi) == 4) rthe <- rthe +1
      # window sizes
      if(length(temp$ws) != (rsig + ralp + rthe)) stop(sprintf("Vector of window sizes for %s is not the same length as the number of parameters: check obj$%s$ws", level))
      # initial values
      if(any(temp$init[temp$Pi[2,] == 1] <= 0)) stop(sprintf("initial values for %s should not be negative: check obj$%s$init", level))
      if(length(temp$init) != (rsig + ralp + rthe)) stop(sprintf("Vector of initial values for %s is not the same length as the number of parameters: check obj$%s$init", level))
      # hyperprior
      for(i in ralp + 1:rsig){
        if(is.finite(temp$hprior[[2]](-1)[[1]])){
          stop(sprintf("Hyper prior should not allow sigma^2 <= 0: check obj$%s$hprior$%s", level, names(temp$hprior[i])))
        } 
      }
    }
  }
}