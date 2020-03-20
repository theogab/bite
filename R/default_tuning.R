#' @import stats
## Default tuning for Jive analysis

default_tuning <- function(model.mean = c("BM", "OU", "WN"),
                           model.var = c("BM", "OU", "WN"),
                           phy, traits, map){
  
  var.sp <- sapply(traits, var, na.rm = T)
  var.sp[is.na(var.sp)] <- var(var.sp, na.rm = T)
  mean.sp <- sapply(traits, mean, na.rm = T)
  ran.sp <- sapply(traits, range, na.rm = T)
  counts <- sapply(traits, function(x) sum(!is.na(x)))
  
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
    if (!model.evo[1] %in% c("OU", "BM", "WN")){
      stop(sprintf("in %s: specified model is not supported", level))
    }
    
    # White Noise #
    if ("WN" %in% model.evo){
      # model
      model <- update_bm
      # map and nreg
      if("sigma" %in% model.evo){
        nreg <- ncol(map$beta)
        newmap <- map
      } else {
        nreg <- 1
        newmap <- input_to_map(phy, nreg = nreg)
      }
      # name
      name <- paste(paste(model.evo, collapse = " + ")," [",nreg,"]", sep=" ")
      rname <- list("",sprintf("_%s",colnames(newmap$beta)))
      
      # indicator matrix for parameters
      Pi <- matrix(0, 2, 1 + nreg)
      Pi[1,1:nreg] <- 1
      Pi[2,nreg + 1] <- 1
      
      # window size
      ws <- c(rep(2, nreg), sd(x))
      
      # initial parameter values
      init <- c(runif(nreg, 0.5, 3), mean(x))
      
      # proposals
      prop <- lapply(1:nreg, proposal, prop = "multiplierProposal") # sigma(s)
      prop[[nreg+1]] <- proposal("slidingWin") # root
      
      # hyper priors
      hprior <- lapply(1:nreg, hpfun, hpf = "Gamma", hp.pars = c(1.1,5)) # sigma(s)
      bounds <- c(min(x) - abs(min(x)),max(x) + abs(max(x)))
      hprior[[nreg+1]] <- hpfun("Uniform", bounds) # theta
      names(hprior) <- names(ws) <- names(init) <- names(prop) <- c(sprintf("sigma^2%s", rname[[1+Pi[1,2]]]), "root")
      
    } 
    
    # Brownian Motion #
    if ("BM" %in% model.evo){
      # model
      model <- update_bm
      # map and nreg
      if("sigma" %in% model.evo){
        nreg <- ncol(map$beta)
        newmap <- map
      } else {
        nreg <- 1
        newmap <- input_to_map(phy, nreg = nreg)
      }
      # name
      name <- paste(paste(model.evo, collapse = " + ")," [",nreg,"]", sep=" ")
      rname <- list("",sprintf("_%s",colnames(newmap$beta)))
      
      # indicator matrix for parameters
      Pi <- matrix(0, 2, 1 + nreg)
      Pi[1,1:nreg] <- 1
      Pi[2,nreg + 1] <- 1
      
      # window size
      ws <- c(rep(2, nreg), sd(x))
      
      # initial parameter values
      init <- c(runif(nreg, 0.5, 3), mean(x))
      
      # proposals
      prop <- lapply(1:nreg, proposal, prop = "multiplierProposal") # sigma(s)
      prop[[nreg+1]] <- proposal("slidingWin") # root
      
      # hyper priors
      hprior <- lapply(1:nreg, hpfun, hpf = "Gamma", hp.pars = c(1.1,5)) # sigma(s)
      bounds <- c(min(x) - abs(min(x)),max(x) + abs(max(x)))
      hprior[[nreg+1]] <- hpfun("Uniform", bounds) # theta
      names(hprior) <- names(ws) <- names(init) <- names(prop) <- c(sprintf("sigma^2%s", rname[[1+Pi[1,2]]]), "root")
    }
    
    # Ornstein-Uhlenbeck #
    if ("OU" %in% model.evo){
      # model
      model <- update_ou
      # map and nreg
      if(any(c("sigma", "theta", "sv") %in% model.evo)){
        nreg <-  ncol(map$beta)
        newmap <- map
      } else {
        nreg <- 1
        newmap <- input_to_map(phy, nreg = nreg)
      }
      # name
      name <- paste(paste(model.evo, collapse = " + ")," [",nreg,"]", sep=" ")
      rname <- list("",sprintf("_%s",colnames(newmap$beta)))
      
      # number of regimes for each parameter
      rsv <- ifelse(any(c("sv", "sigma") %in% model.evo), nreg, 1)
      rsig <- ifelse("sigma" %in% model.evo, nreg, 1)
      rthe <- ifelse("theta" %in% model.evo, nreg, 1)
      rroot <- ifelse("root" %in% model.evo, 1, 0)
      
      # indicator matrix for parameters
      Pi <- matrix(0, 3+rroot, rsv + rsig + rthe + rroot)
      Pi[1,1:rsv] <- 1
      Pi[2,rsv + 1:rsig] <- 1
      Pi[3, rsv + rsig + 1:rthe] <- 1
      if(rroot == 1) Pi[4,ncol(Pi)] <- 1
      
      # window size
      ws <- c(rep(0.5, rsv), rep(2, rsig), rep(sd(x), rroot+rthe)) # 2 in the previous version?
      
      # initial parameter values
      init<-c(runif(rsv, 0.1, 1), runif(rsig, 0.5, 3), rep(mean(x), rroot+rthe))
      
      # proposals
      prop <- list()
      i <- 1
      while(i <= (rsv + rsig + rroot + rthe)){
        if (i <= rsv) prop[[i]] <- proposal("multiplierProposal")
        else if (i <= (rsv + rsig)) prop[[i]] <- proposal("multiplierProposal")
        else prop[[i]] <- proposal("slidingWin")
        i <- i + 1
      }
      
      # hyper priors
      hprior <- list()
      i <- 1
      bounds <- c(min(x) - abs(min(x)),max(x) + abs(max(x)))
      while(i <= (rsv + rsig + rthe + rroot)){
        if (i <= rsv)  hprior[[i]] <- hpfun("Gamma", c(1.1,5))
        else if (i <= (rsv + rsig)) hprior[[i]]	<- hpfun("Gamma", c(1.1,5))
        else hprior[[i]]	<- hpfun("Uniform", bounds)
        i <- i + 1
      }
      
      names(hprior) <- names(ws) <- names(init) <- names(prop) <- c(sprintf("sv%s", rname[[1+Pi[1,2]]]), 
                                                                    sprintf("sigma^2%s", rname[[1+Pi[2,rsv + 2]]]),
                                                                    sprintf("theta%s", rname[[1+ifelse(sum(Pi[3,])>1,1,0)]]),
                                                                    "root"[rroot])
    }
    
    eval(parse(text = sprintf("%s <- list(name = name, model = model, Pi = Pi, ws = ws, init = init, prop = prop, map = newmap, hprior = hprior, update.freq = update.freq)", level)))

  }
  
  
  return(eval(parse(text = "list(lik = lik, prior.mean = prior.mean, prior.var = prior.var)")))
  
}
