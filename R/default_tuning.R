#' @import stats
## Default tuning for Jive analysis

default_tuning <- function(model.priors, phy, traits, map, init.lik, lik.f = NULL){
  
  ### Likelihood level ###
  # update.freq
  update.freq <- 0.3
  # number of species to update at the same time
  n.u <- 3
  
  if(is.null(lik.f)){
    # model
    model <- lik_norm
  } else {
    # model
    model <- lik.f
  }
    
  init <- list()
  ws <- list()
  prop <- list()
  
  for(i in 1:length(model.priors)){
    init[[i]]  <- init.lik[i,]
    ws[[i]] <- rep(sd(init.lik[i,]), ncol(init.lik))
    prop[[i]] <- proposal("slidingWin")
  }
  
  names(ws) <- names(init) <- names(prop) <- names(model.priors)
  
  lik <- list(model = model, ws = ws, init = init, prop = prop, update.freq = update.freq, n.u = n.u)
    
  
  ### Prior level ###
  priors <- list()
  
  for(i in 1:length(model.priors)){

    x <- init.lik[i,]
    update.freq <- 0.7/length(model.priors)
    
    ## check model specification
    if (class(model.priors[[i]]) == "function"){
      name <- "User defined function"
      warning(sprintf("Alternative evolutionary model for %s: You must Add MCMC tuning information using control_jive()", names(model.priors)[i]))
      model <- lik_gen(phy, model.priors[[i]])
      Pi <- NULL
      init <- ws <- c()
      prop <- hprior <- list()
      newmap <- NULL
    } else {
     
      # White Noise #
      if ("WN" %in% model.priors[[i]]){
        # model
        model <- lik_bm
        # map and nreg
        if("sigma" %in% model.priors[[i]]){
          nreg <- ncol(map$beta)
          newmap <- map
        } else {
          nreg <- 1
          newmap <- input_to_map(phy, nreg = nreg)
        }
        # name
        name <- paste(paste(model.priors[[i]], collapse = " + ")," [",nreg,"]", sep="")
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
        names(hprior) <- names(ws) <- names(init) <- names(prop) <- c(sprintf("sigma_sq%s", rname[[1+Pi[1,2]]]), "root")
        
      } 
      
      # Brownian Motion #
      if ("BM" %in% model.priors[[i]]){
        # model
        model <- lik_bm
        # map and nreg
        if("sigma" %in% model.priors[[i]]){
          nreg <- ncol(map$beta)
          newmap <- map
        } else {
          nreg <- 1
          newmap <- input_to_map(phy, nreg = nreg)
        }
        # name
        name <- paste(paste(model.priors[[i]], collapse = " + ")," [",nreg,"]", sep="")
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
        names(hprior) <- names(ws) <- names(init) <- names(prop) <- c(sprintf("sigma_sq%s", rname[[1+Pi[1,2]]]), "root")
      }
      
      # Ornstein-Uhlenbeck #
      if ("OU" %in% model.priors[[i]]){
        # model
        model <- lik_ou
        # map and nreg
        if(any(c("sigma", "theta", "aplha") %in% model.priors[[i]])){
          nreg <-  ncol(map$beta)
          newmap <- map
        } else {
          nreg <- 1
          newmap <- input_to_map(phy, nreg = nreg)
        }
        # name
        name <- paste(paste(model.priors[[i]], collapse = " + ")," [",nreg,"]", sep="")
        rname <- list("",sprintf("_%s",colnames(newmap$beta)))
        
        # number of regimes for each parameter
        ralp <- ifelse(any(c("alpha", "sigma") %in% model.priors[[i]]), nreg, 1)
        rsig <- ifelse("sigma" %in% model.priors[[i]], nreg, 1)
        rthe <- ifelse("theta" %in% model.priors[[i]], nreg, 1)
        rroot <- ifelse("root" %in% model.priors[[i]], 1, 0)
        
        # indicator matrix for parameters
        Pi <- matrix(0, 3+rroot, ralp + rsig + rthe + rroot)
        Pi[1,1:ralp] <- 1
        Pi[2,ralp + 1:rsig] <- 1
        Pi[3, ralp + rsig + 1:rthe] <- 1
        if(rroot == 1) Pi[4,ncol(Pi)] <- 1
        
        # window size
        ws <- c(rep(0.5, ralp), rep(2, rsig), rep(sd(x), rroot+rthe)) # 2 in the previous version?
        
        # initial parameter values
        init<-c(runif(ralp, 0.1, 1), runif(rsig, 0.5, 3), rep(mean(x), rroot+rthe))
        
        # proposals
        prop <- list()
        j <- 1
        while(j <= (ralp + rsig + rroot + rthe)){
          if (j <= ralp) prop[[j]] <- proposal("multiplierProposal")
          else if (j <= (ralp + rsig)) prop[[i]] <- proposal("multiplierProposal")
          else prop[[j]] <- proposal("slidingWin")
          j <- j + 1
        }
        
        # hyper priors
        hprior <- list()
        j <- 1
        bounds <- c(min(x) - abs(min(x)),max(x) + abs(max(x)))
        while(j <= (ralp + rsig + rthe + rroot)){
          if (j <= ralp)  hprior[[j]] <- hpfun("Gamma", c(1.1,5))
          else if (j <= (ralp + rsig)) hprior[[j]]	<- hpfun("Gamma", c(1.1,5))
          else hprior[[j]]	<- hpfun("Uniform", bounds)
          j <- j + 1
        }
        
        names(hprior) <- names(ws) <- names(init) <- names(prop) <- c(sprintf("alpha%s", rname[[1+Pi[1,2]]]), 
                                                                      sprintf("sigma_sq%s", rname[[1+Pi[2,ralp + 2]]]),
                                                                      sprintf("theta%s", rname[[1+ifelse(sum(Pi[3,])>1,1,0)]]),
                                                                      "root"[rroot])
      }
      
    }
    
   
    priors[[i]] <- list(name = name, model = model, Pi = Pi, ws = ws, init = init, prop = prop, map = newmap, hprior = hprior, update.freq = update.freq)

  }
  
  names(priors) <- names(model.priors)
  
  return(list(lik = lik, priors = priors))
}
