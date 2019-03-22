## Default tuning for Jive analysis

default_tuning <- function(model.mean = c("BM", "OU", "WN", "OUM", "BMM", "WNM"),
                           model.var = c("BM", "OU", "WN", "OUM", "BMM", "WNM"),
                           traits, map, root.station = F){
  
  var.sp <- apply(traits, 1, var, na.rm = T)
  var.sp[is.na(var.sp)] <- var(var.sp, na.rm = T)
  mean.sp <- apply(traits, 1, mean, na.rm = T)
  ran.sp <- apply(traits, 1, range, na.rm = T)
  counts <- apply(traits, 1, function(x) sum(!is.na(x)))
  
  ### Likelihood level ###
  # update.freq
  update.freq <- 0.35
  # model
  model <- lik_multinorm
  
  # window size
  ws <- list()
  ws$m.sp <- var.sp
  ws$v.sp <- ran.sp[1,]/ran.sp[2,]
  # initial parameter value
  init  <- list()
  init$m.sp <- mean.sp
  init$v.sp <- var.sp
  
  # proposals
  prop <- list()
  prop$m.sp <- proposal("slidingWin") 
  prop$v.sp <- proposal("logSlidingWinAbs")
  
  lik <- list(model = model, ws = ws, init = init, prop = prop, update.freq = update.freq)
    
  
  ### Prior level ###
  for(level in c("prior.mean", "prior.var")){

    if (level == "prior.mean"){ ## Mean prior level ##
      x <- mean.sp
      update.freq <- 0.2
      model.evo <- model.mean
    } else { ## Var prior level ##
      x <- log(var.sp)
      update.freq <- 0.45
      model.evo <- model.var
    }
    
    
    ## check model specification
    if (!model.evo %in% c("OU", "BM", "WN", "OUM", "BMM", "WNM")){
      stop(sprintf("in %s: model %s is not supported", level, model.evo))
    }
    
    # White Noise #
    if (grepl("WN", model.evo)){
      # model
      model <- update_wn
      # map and nreg
      if(model.evo == "WNM"){
        nreg <- ncol(map)
      } else {
        nreg <- 1
        map <- as.matrix(rowSums(map))
      }
      # name
      name <- paste(model.evo," [",nreg,"]", sep="")
      # window size
      ws <- list()
      ws$wn.sig <- rep(0.5, nreg)
      ws$wm.the <- sd(x) # 2 in the previous version?
      # initial parameter values
      init <- list()
      init$wn.sig <- runif(nreg, 0.5, 3)
      init$wn.the <- mean(x)
      # proposals
      prop <- lapply(1:nreg, proposal, prop = "multiplierProposalLakner") # sigma(s)
      prop[[nreg+1]] <- proposal("slidingWin") # theta
      # hyper priors
      hprior <- lapply(1:nreg, hpfun, hpf = "Gamma", hp.pars = c(1.1,5)) # sigma(s)
      bounds <- c(ifelse(min(x) < 0, 2*min(x),min(x)/2),ifelse(max(x) < 0, max(x)/2,2*max(x)))
      hprior[[nreg+1]] <- hpfun("Uniform", bounds) # theta
      names(hprior) <- c(if(nreg==1) "wn.sig" else sprintf("wn.sig.%s", 1:nreg), "wn.the")
    } 
    
    # Brownian Motion #
    if (grepl("BM", model.evo)){
      # model
      model <- update_bm
      # map and nreg
      if(model.evo == "BMM"){
        nreg <- ncol(map)
      } else {
        nreg <- 1
        map <- as.matrix(rowSums(map))
      }
      # name
      name <- paste(model.evo," [",nreg,"]", sep="")
      # window size
      ws <- list()
      ws$bm.sig <- rep(2, nreg)
      ws$bm.the <- sd(x)
      # initial parameter values
      init <- list()
      init$bm.sig <- runif(nreg, 0.5, 3)
      init$bm.the <- mean(x)
      # proposals
      prop <- lapply(1:nreg, proposal, prop = "multiplierProposalLakner") # sigma(s)
      prop[[nreg+1]] <- proposal("slidingWin") # theta
      # hyper priors
      hprior <- lapply(1:nreg, hpfun, hpf = "Gamma", hp.pars = c(1.1,5)) # sigma(s)
      bounds <- c(ifelse(min(x) < 0, 2*min(x),min(x)/2),ifelse(max(x) < 0, max(x)/2,2*max(x)))
      hprior[[nreg+1]] <- hpfun("Uniform", bounds) # theta
      names(hprior) <- c(if(nreg==1) "bm.sig" else sprintf("bm.sig.%s", 1:nreg), "bm.the")
    }
    
    # Ornstein-Uhlenbeck #
    if (grepl("OU", model.evo)){
      # model
      model <- update_ou
      # map and nreg
      if(model.evo == "OUM"){
        nreg <- ncol(map)
      } else {
        nreg <- 1
        map <- as.matrix(rowSums(map))
      }
      # name
      name <- paste(model.evo," [",nreg,"]", sep="")
      # window size
      ws <- list()
      ws$ou.sv <- 0.5
      ws$ou.sig <- 2
      ws$ou.the <- rep(sd(x), nreg + ifelse(root.station, 0, 1)) # 2 in the previous version?
      # initial parameter values
      init <- list()
      init$ou.sv <- runif(1, 0.1, 1)
      init$ou.sig <- runif(1, 0.5, 3)
      init$ou.the <- rep(mean(x), nreg + ifelse(root.station, 0, 1))
      # proposals
      prop <- list()
      prop[[1]]	<- proposal("multiplierProposalLakner")
      prop[[2]]	<- proposal("multiplierProposalLakner")
      for(i in 3:(nreg+ ifelse(root.station, 2, 3))){
        prop[[i]]	<- proposal("slidingWin")
      }
      # hyper priors
      hprior <- list()
      hprior[[1]] <- hpfun("Gamma", c(1.1,5))
      hprior[[2]]	<- hpfun("Gamma", c(1.1,5))
      bounds <- c(ifelse(min(x) < 0, 2*min(x),min(x)/2),ifelse(max(x) < 0, max(x)/2,2*max(x)))
      for(i in 3:(nreg + ifelse(root.station, 2, 3))){
        hprior[[i]]	<- hpfun("Uniform", bounds) ## <- test loggamma??
      }
      names(hprior) <- c("ou.sv", "ou.sig", sprintf("ou.the.%s", ifelse(root.station,1,0):nreg))
    }
    
    eval(parse(text = sprintf("%s <- list(name = name, model = model, ws = ws, init = init, prop = prop, map = map, hprior = hprior, update.freq = update.freq)", level)))

  }
  
  return(list(lik = lik, prior.mean = prior.mean, prior.var = prior.var))
  
}
