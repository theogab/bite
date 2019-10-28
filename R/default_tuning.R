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
      model <- update_wn
      # map and nreg
      if("sigma" %in% model.evo){
        nreg <- max(do.call(cbind,map)[1,])
        newmap <- map
      } else {
        nreg <- 1
        newmap <- input_to_map(phy, nreg = nreg)
      }
      # name
      name <- paste(paste(model.evo, collapse = " + ")," [",nreg,"]", sep=" ")
      nr <- c(nreg,1)
      # window size
      ws <- list()
      ws$wn.sig <- rep(0.5, nreg)
      ws$wm.the <- sd(x) # 2 in the previous version?
      # initial parameter values
      init <- list()
      init$wn.sig <- runif(nreg, 0.5, 3)
      init$wn.the <- mean(x)
      # proposals
      prop <- lapply(1:nreg, proposal, prop = "multiplierProposal") # sigma(s)
      prop[[nreg+1]] <- proposal("slidingWin") # theta
      # hyper priors
      hprior <- lapply(1:nreg, hpfun, hpf = "Gamma", hp.pars = c(1.1,5)) # sigma(s)
      bounds <- c(min(x) - abs(min(x)),max(x) + abs(max(x)))
      hprior[[nreg+1]] <- hpfun("Uniform", bounds) # theta
      names(hprior) <- c(if(nreg==1) "wn.sig" else sprintf("wn.sig.%s", 1:nreg), "wn.the")
    } 
    
    # Brownian Motion #
    if ("BM" %in% model.evo){
      # model
      model <- update_bm
      # map and nreg
      if("sigma" %in% model.evo){
        nreg <- max(do.call(cbind,map)[1,])
        newmap <- map
      } else {
        nreg <- 1
        newmap <- input_to_map(phy, nreg = nreg)
      }
      # name
      name <- paste(paste(model.evo, collapse = " + ")," [",nreg,"]", sep=" ")
      nr <- c(nreg,1)
      # window size
      ws <- list()
      ws$bm.sig <- rep(2, nreg)
      ws$bm.the <- sd(x)
      # initial parameter values
      init <- list()
      init$bm.sig <- runif(nreg, 0.5, 3)
      init$bm.the <- mean(x)
      # proposals
      prop <- lapply(1:nreg, proposal, prop = "multiplierProposal") # sigma(s)
      prop[[nreg+1]] <- proposal("slidingWin") # theta
      # hyper priors
      hprior <- lapply(1:nreg, hpfun, hpf = "Gamma", hp.pars = c(1.1,5)) # sigma(s)
      bounds <- c(min(x) - abs(min(x)),max(x) + abs(max(x)))
      hprior[[nreg+1]] <- hpfun("Uniform", bounds) # theta
      names(hprior) <- c(if(nreg==1) "bm.sig" else sprintf("bm.sig.%s", 1:nreg), "bm.the")
    }
    
    # Ornstein-Uhlenbeck #
    if ("OU" %in% model.evo){
      # model
      model <- update_ou
      # map and nreg
      if(any(c("sigma", "theta", "alpha") %in% model.evo)){
        nreg <-  max(do.call(cbind,map)[1,])
        newmap <- map
      } else {
        nreg <- 1
        newmap <- input_to_map(phy, nreg = nreg)
      }
      # name
      name <- paste(paste(model.evo, collapse = " + ")," [",nreg,"]", sep=" ")
      # number of regimes for each parameter
      rsv <- ifelse(any(c("alpha", "sigma") %in% model.evo), nreg, 1)
      rsig <- ifelse("sigma" %in% model.evo, nreg, 1)
      rthe <- ifelse("theta" %in% model.evo, nreg, 1)
      rroot <- ifelse("root" %in% model.evo, 1, 0)
      nr <- c(rsv, rsig, rroot, rthe)
      # window size
      ws <- list()
      ws$ou.sv <- rep(0.5, rsv)
      ws$ou.sig <- rep(2, rsig)
      ws$ou.the <- rep(sd(x), rroot+rthe) # 2 in the previous version?
      # initial parameter values
      init <- list()
      init$ou.sv <- runif(rsv, 0.1, 1)
      init$ou.sig <- runif(rsig, 0.5, 3)
      init$ou.the <- rep(mean(x), rroot+rthe)
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
      while(i <= (rsv + rsig + rroot + rthe)){
        if (i <= rsv)  hprior[[i]] <- hpfun("Gamma", c(1.1,5))
        else if (i <= (rsv + rsig)) hprior[[i]]	<- hpfun("Gamma", c(1.1,5))
        else hprior[[i]]	<- hpfun("Uniform", bounds)
        i <- i + 1
      }
      names(hprior) <- c(sprintf("ou.sv.%s", 1:rsv), sprintf("ou.sig.%s", 1:rsig), sprintf("ou.the.%s", if("root" %in% model.evo) 0:(rthe) else 1:rthe))
    }
    
    eval(parse(text = sprintf("%s <- list(name = name, model = model, ws = ws, init = init, prop = prop, map = newmap, hprior = hprior, nr = nr, update.freq = update.freq)", level)))

  }
  
  return(list(lik = lik, prior.mean = prior.mean, prior.var = prior.var))
  
}
