make_control_jive <- function(level = c("lik", "prior.mean", "prior.var"), model = c("OU", "BM", "WN", "OUM", "BMM", "WNM"),
                          traits, nreg = 1, initial.ws = NULL, initial.pv = NULL, proposals = NULL, hyperprior = NULL){
  
  
  var.sp <- apply(traits, 1, sd, na.rm = T)
  mean.sp <- apply(traits, 1, mean, na.rm = T)
  
  ### Likelihood level ###
  if (level == "lik"){
    # model
    model <- lik_multinorm
    
    # window size
    ws <- list()
    if (is.null(initial.ws)){
      ws$m.sp <- 2*mean.sp
      ws$v.sp <- 10*var.sp
    } else {
      ws$m.sp <- initial.ws[,1]
      ws$v.sp <- initial.ws[,2]
    }
    #checking
    if(!length(ws$m.sp) == dim(traits)[1] | !length(ws$v.sp) == dim(traits)[1] | any(is.na(c(ws$m.sp, ws$v.sp)))){
      stop("initial.ws is not valid")
    }
    
    # initial parameter value
    pv  <- list()
    if (is.null(inital.pv)){
      pv$m.sp <- mean.sp
      pv$v.sp <- var.sp
    } else {
      pv$m.sp <- initial.pv[,1]
      pv$v.sp <- initial.pv[,2]
    }
    #checking
    if(!length(pv$m.sp) == dim(traits)[1] | !length(pv$v.sp) == dim(traits)[1] | any(is.na(c(pv$m.sp, pv$v.sp)))){
      stop("initial.pv is not valid")
    }
    
    # proposals
    prop <- list()
    if (is.null(prop)){
      prop$m.sp <- make_proposal("slidingWin") 
      prop$v.sp <- make_proposal("logSlidingWinAbs")
    } else {
      prop$m.sp <- make_proposal(prop[1])
      prop$v.sp	<- make_proposal(prop[2])
    }
    
    return(list(model = model, ws = ws, pv = pv, prop = prop))
    
  }
  
  
  ### Prior level ###
  if (grepl("prior", level)){
    
    ## Mean prior level ##    
    if (level == "prior.mean"){
      x <- mean.sp
    } else { ## Var prior level ##
      x <- var.sp
    }
    
    # White Noise #
    if (grepl("WN", model)){
      # model
      model <- lik_wn
      # window size
      if(is.null(initial.ws)){
        ws <- list()
        ws$sig.wn <- rep(0.5, nreg)
        ws$the.wn <- sd(x) # 2 in the previous version?
      } else {
        ws <- list()
        ws$sig.wn <- initial.ws[1:nreg]
        ws$the.wn <- initial.ws[nreg+1]
      }
      #checking
      if(any(is.na(c(ws$sig.wn, ws$the.wn))) | any(c(ws$sig.wn, ws$the.wn) <= 0)){
        stop("initial.ws is not valid")
      }
      
      # initial parameter values
      if(is.null(initial.pv)){
        pv <- list()
        pv$sig.wn <- runif(nreg, 0.5, 3)
        pv$the.wn <- mean(x)
      }
      #checking
      if(any(is.na(c(pv$sig.wn, pv$the.wn))) | any(c(pv$sig.wn, pv$the.wn) <= 0)){
        stop("initial.pv is not valid")
      }
      
      # proposals
      if (is.null(proposals)){
        prop <- list()
        prop$sig.wn	<- make_proposal("multiplierProposalLakner")
        prop$the.wn	<- make_proposal("slidingWin")
      } else {
        prop$sig.wn <- make_proposal(prop[1])
        prop$the.wn	<- make_proposal(prop[2])
      }
      
      # hyper priors
      if (is.null(hyperprior)){
        hprior <- list()
        hprior$sig.wn		<- hpfun("Gamma", c(1.1,5))
        hprior$the.wn		<- hpfun("Uniform", c(-20,10)) 
      }
    } 
    # checking
    if(is.finite(hprior$sig.wn(-1))|is.finite(hprior$sig.wn(-1))){
      stop("Hyper prior should not allow sigma <= 0")
    }
    
    # Brownian Motion #
    if (grepl("BM", model)){
      # model
      model <- lik_bm
      # window size
      if(is.null(initial.ws)){
        ws <- list()
        ws$sig.bm <- rep(2, nreg)
        ws$the.bm <- sd(x)
      }
      # initial parameter values
      if(is.null(initial.pv)){
        pv <- list()
        pv$sig.bm <- runif(nreg, 0.5, 3)
        pv$the.bm <- mean(x)
      }
      # proposals
      if (is.null(proposals)){
        prop <- list()
        prop$sig.bm	<- make_proposal("multiplierProposalLakner")
        prop$the.bm	<- make_proposal("slidingWin")
      } else {
        prop$sig.bm <- make_proposal(prop[1])
        prop$the.bm	<- make_proposal(prop[2])
      }
      # hyper priors
      if (is.null(hyperprior)){
        hprior <- list()
        hprior$sig.bm	<- hpfun("Gamma", c(1.1,5))
        hprior$the.bm	<- hpfun("Uniform", c(-20,10)) ## <- test loggamma??
      }
    }
    
    # Ornstein-Uhlenbeck #
    if (grepl("OU", model)){
      # model
      model <- lik_ou
      # window size
      if(is.null(initial.ws)){
        ws <- list()
        ws$alp.ou <- 0.5
        ws$sig.ou <- 2
        ws$the.ou <- rep(sd(x), nreg+1) # 2 in the previous version?
      }
      # initial parameter values
      if(is.null(initial.pv)){
        pv <- list()
        pv$alp.ou <- runif(1, 0.1, 1)
        pv$sig.ou <- runif(1, 0.5, 3)
        pv$the.ou <- rep(mean(x), nreg+1)
      }
      # proposals
      prop <- list()
      if (is.null(proposals)){
        prop$alp.ou	<- make_proposal("multiplierProposalLakner")
        prop$sig.ou	<- make_proposal("multiplierProposalLakner")
        prop$the.ou	<- make_proposal("slidingWin")
      } else {
        prop$alp.ou <- make_proposal(prop[1])
        prop$sig.ou	<- make_proposal(prop[2])
        prop$the.ou	<- make_proposal(prop[3])
      }
      # hyper priors
      if (is.null(hyperprior)){
        hprior <- list()
        hprior$alp.ou <- hpfun("Gamma", c(1.1,5))
        hprior$sig.ou	<- hpfun("Gamma", c(1.1,5))
        hprior$the.ou	<- hpfun("Uniform", c(-20,10)) ## <- test loggamma??
      }
    }
    
    return(list(model = model, ws = ws, pv = pv, prop = prop, hprior = hprior))
    
  }
}



