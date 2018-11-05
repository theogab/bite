#' @title Create a list that can be used to tune jive mcmc algorithm 
#' @description This function creates an object to parse in the control argument of the \code{\link{make_jive}} function. The output will be different regarding which level of the jive model the user wants to tune ($lik, $prior.mean, $prior.var). This function allows tuning of : initial window size for proposals, starting parameter value, proposal methods and Hyperprior specification 
#' @details If arguments initial.ws, initial.pv, proposals or hyperprior are left blank, the default tuning that we found appropriate for most datasets is applied
#' 
#' If level == "lik"
#' initial.ws and initial.pv must be entered as a matrix with 2 columns (respectively mean and variance) and a number of rows equal to the number of species. proposal must be a vector of size 2 (respectively mean and variance)
#' 
#' If level == "prior.mean" or "prior.var"
#' initial.ws and initial.pv must be entered as a vector of variable size depending on the chosen evolutionary model. for OU and OUM, the window size and parameter values must be entered in the following order c(alpha, sigma, theta0, theta1, ..., thetaN). for BM, BMM, WN and WNM, the window size and parameter values must be entered in the following order c(sigma1, ..., sigmaN, theta0), proposal must be a vector of size three for OU and OUM c(alpha, sigma, thetas) and of size two for BM, BMM, WN and WNM c(sigmas, theta)
#' 
#' proposals
#' Has to be one the following : "slidingWin" for Sliding window proposal unconstrained at maximum, "multiplierProposal", for multiplier proposal
#' 
#' Hyperprior
#' list of hyperpriror functions (see \code{\link{hpfun}}). User must provide a list of size 2 for BM, BMM, WN and WNM (sigmas, theta0) and of size 3 for OU and OUM (alpha, sigma, thetas)
#' 
#' @param level character taken in c("lik", "prior.mean", "prior.var") to specify on which level of the jive model, the control will operate (see details)
#' @param model.evo character taken in c("OU", "BM", "WN", "OUM", "BMM", "WNM") specifying the evolutionary model. ignored if level == "lik"
#' @param traits matrix of traits value (see details)
#' @param map matrix mapping regimes on every edge of phy
#' @param window.size initial window size for proposals during the mcmc algorithm. matrix or vector depending on the value of level and nreg (see details)
#' @param initial.values starting parameter values of the mcmc algorithm. matrix or vector depending on the value of level and nreg (see details)
#' @param proposals vector of characters taken in c("slidingWin", "slidingWinAbs", "logSlidingWinAbs","multiplierProposal", "multiplierProposalLakner","logNormal", "absNormal") to control proposal methods during mcmc algorithm (see details)
#' @param hyperprior list of hyperprior functions that can be generated with \code{\link{hpfun}}function. Ignored if level == "lik" (see details)
#' @param root.station boolean indicating whether the theta_0 should be dropped from the OU or OUM models 
#' @export
#' @author Théo Gaboriau
#' @return A list to parse into control argument of \code{\link{make_jive}} function. The list is containing the following objects:
#' $model : a function to calculate the likelihood ($lik) or the priors ($prior.mean, $prior_var)
#' $ws : a list containing the window size of proposals for each estimated parameter
#' $init : a list containing the starting value of the mcmc chain for each parameter
#' $prop : a list containing the proposal functions for the update of each parameter
#' $map (only if level %in% c(prior.mean, prior.var) : a matrix containing the mapping of regimes onto branch for a specific model
#' $hprior (only if level %in% c(prior.mean, prior.var) : a list containing the hyperprior function for each parameter 
#' 
#' @examples
#' 

control_jive <- function(level = c("lik", "prior.mean", "prior.var"), model.evo = c("BM", "OU", "WN", "OUM", "BMM", "WNM"),
                         traits, map, window.size = NULL, initial.values = NULL, proposals = NULL, hyperprior = NULL, root.station = F){
  
  var.sp <- apply(traits, 1, sd, na.rm = T)
  mean.sp <- apply(traits, 1, mean, na.rm = T)
  
  ### Likelihood level ###
  if (level == "lik"){
    # model
    model <- lik_multinorm
    
    # window size
    ws <- list()
    if (is.null(window.size)){
      ws$m.sp <- 2*mean.sp
      ws$v.sp <- 10*var.sp
    } else {
      ws$m.sp <- window.size[,1]
      ws$v.sp <- window.size[,2]
    }
    #checking
    if(!length(ws$m.sp) == dim(traits)[1] | !length(ws$v.sp) == dim(traits)[1] | any(is.na(c(ws$m.sp, ws$v.sp)))){
      stop("window.size is not valid")
    }
    
    # initial parameter value
    init  <- list()
    if (is.null(initial.values)){
      init$m.sp <- mean.sp
      init$v.sp <- var.sp
    } else {
      init$m.sp <- initial.values[,1]
      init$v.sp <- initial.values[,2]
    }
    #checking
    if(!length(init$m.sp) == dim(traits)[1] | !length(init$v.sp) == dim(traits)[1] | any(is.na(c(init$m.sp, init$v.sp)))){
      stop("initial.values is not valid")
    }
    
    # proposals
    prop <- list()
    if (is.null(proposals)){
      prop$m.sp <- proposal("slidingWin") 
      prop$v.sp <- proposal("logSlidingWinAbs")
    } else {
      prop$m.sp <- proposal(proposals[1])
      prop$v.sp	<- proposal(proposals[2])
    }
    
    return(list(model = model, ws = ws, init = init, prop = prop))
    
  }
  
  
  ### Prior level ###
  if (grepl("prior", level)){
    
    ## check model specification
    if (!model.evo %in% c("OU", "BM", "WN", "OUM", "BMM", "WNM")){
      stop(paste("model",model,"is not supported", sep = " "))
    }
    
    
    ## Mean prior level ##    
    if (level == "prior.mean"){
      x <- mean.sp
    } else { ## Var prior level ##
      x <- var.sp
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
      # window size
      if(is.null(window.size)){
        ws <- list()
        ws$wn.sig <- rep(0.5, nreg)
        ws$wm.the <- sd(x) # 2 in the previous version?
      } else {
        ws <- list()
        ws$wn.sig <- window.size[1:nreg]
        ws$wn.the <- window.size[nreg+1]
      }
      # initial parameter values
      if(is.null(initial.values)){
        init <- list()
        init$wn.sig <- runif(nreg, 0.5, 3)
        init$wn.the <- mean(x)
      } else {
        init <- list()
        init$wn.sig <- initial.values[1:nreg]
        init$wn.the <- initial.values[nreg]
      }
      # proposals
      if (is.null(proposals)){
        prop <- lapply(1:nreg, proposal, prop = "multiplierProposalLakner") # sigma(s)
        prop[[nreg+1]] <- proposal("slidingWin") # theta
      } else {
        prop <- lapply(1:nreg, proposal, prop = proposals[1]) # sigma(s)
        prop[[nreg+1]]	<- proposal(proposals[2]) # theta
      }
      # hyper priors
      if (is.null(hyperprior)){
        hprior <- lapply(1:nreg, hpfun, hpf = "Gamma", hp.pars = c(1.1,5)) # sigma(s)
        hprior[[nreg+1]] <- hpfun("Uniform", c(-20,10)) # theta
      } else {
        hprior <- lapply(1:nreg, function(x) hyperprior[[1]]) # sigma(s)
        hprior[[nreg+1]] <- hyperprior[[2]] # theta
      }
      # checking
      if(model.evo %in% "WNM" & ncol(map) != (length(do.call(c, init)) - 1)){
        stop("in White Noise with Multiple regimes, number of parameters doesn't correspond to the number of regimes in map")
      }
      if(any(is.na(c(ws$wn.sig, ws$wn.the))) | any(c(ws$wn.sig) <= 0)){
        stop("window.size is not valid")
      }
      if(any(is.na(c(init$wn.sig, init$wn.the))) | any(c(init$wn.sig) <= 0)){
        stop("initial.values is not valid")
      }
      if(is.finite(hprior[[1]](-1))){
        stop("Hyper prior should not allow sigma <= 0")
      }
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
      # window size
      if(is.null(window.size)){
        ws <- list()
        ws$bm.sig <- rep(2, nreg)
        ws$bm.the <- sd(x)
      } else {
        ws <- list()
        ws$bm.sig <- window.size[1:nreg]
        ws$bm.the <- window.size[nreg+1]
      }
      # initial parameter values
      if(is.null(initial.values)){
        init <- list()
        init$bm.sig <- runif(nreg, 0.5, 3)
        init$bm.the <- mean(x)
      } else {
        init <- list()
        init$bm.sig <- initial.values[1:nreg]
        init$bm.the <- initial.values[nreg]
      }
      # proposals
      if (is.null(proposals)){
        prop <- lapply(1:nreg, proposal, prop = "multiplierProposalLakner") # sigma(s)
        prop[[nreg+1]] <- proposal("slidingWin") # theta
      } else {
        prop <- lapply(1:nreg, proposal, prop = proposals[1]) # sigma(s)
        prop[[nreg+1]]	<- proposal(proposals[2]) # theta
      }
      # hyper priors
      if (is.null(hyperprior)){
        hprior <- lapply(1:nreg, hpfun, hpf = "Gamma", hp.pars = c(1.1,5)) # sigma(s)
        hprior[[nreg+1]] <- hpfun("Uniform", c(-20,10)) # theta
      } else {
        hprior <- lapply(1:nreg, function(x) hyperprior[[1]]) # sigma(s)
        hprior[[nreg+1]] <- hyperprior[[2]] # theta
      }
      # checking
      if(model.evo %in% "BMM" & ncol(map) != (length(do.call(c, init)) - 1)){
        stop("in Brownian Motion with Multiple regimes, number of parameters doesn't correspond to the number of regimes in map")
      }
      if(any(is.na(c(ws$bm.sig, ws$bm.the))) | any(c(ws$bm.sig) <= 0)){
        stop("window.size is not valid")
      }
      if(any(is.na(c(init$bm.sig, init$bm.the))) | any(c(init$bm.sig) <= 0)){
        stop("initial.values is not valid")
      }
      if(is.finite(hprior[[1]](-1))){
        stop("Hyper prior should not allow sigma <= 0")
      }
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
      # window size
      if(is.null(window.size)){
        ws <- list()
        ws$ou.alp <- 0.5
        ws$ou.sig <- 2
        ws$ou.the <- rep(sd(x), nreg + ifelse(root.station, 0, 1)) # 2 in the previous version?
      } else {
        ws <- list()
        ws$ou.alp <- window.size[1]
        ws$ou.sig <- window.size[2]
        ws$ou.the <- window.size[3:(nreg + ifelse(root.station, 2, 3))]
      }
      # initial parameter values
      if(is.null(initial.values)){
        init <- list()
        init$ou.alp <- runif(1, 0.1, 1)
        init$ou.sig <- runif(1, 0.5, 3)
        init$ou.the <- rep(mean(x), nreg + ifelse(root.station, 0, 1))
      } else {
        init <- list()
        init$ou.alp <- initial.values[1]
        init$ou.sig <- initial.values[2]
        init$ou.the <- initial.values[3:(nreg+ ifelse(root.station, 2, 3))]
      }
      # proposals
      prop <- list()
      if (is.null(proposals)){
        prop[[1]]	<- proposal("multiplierProposalLakner")
        prop[[2]]	<- proposal("multiplierProposalLakner")
        for(i in 3:(nreg+ ifelse(root.station, 2, 3))){
          prop[[i]]	<- proposal("slidingWin")
        }
      } else {
        prop[[1]] <- proposal(proposals[1])
        prop[[2]]	<- proposal(proposals[2])
        for(i in 3:(nreg+ ifelse(root.station, 2, 3))){
          prop[[i]]	<-  proposal(proposals[3])
        }
      }
      # hyper priors
      if (is.null(hyperprior)){
        hprior <- list()
        hprior[[1]] <- hpfun("Gamma", c(1.1,5))
        hprior[[2]]	<- hpfun("Gamma", c(1.1,5))
        for(i in 3:(nreg+ ifelse(root.station, 2, 3))){
          hprior[[i]]	<- hpfun("Uniform", c(-20,10)) ## <- test loggamma??
        }
      } else {
        hprior[[1]] <- hyperprior[[1]]
        hprior[[2]] <- hyperprior[[2]]
        for(i in 3:(nreg+ ifelse(root.station, 2, 3))){
          hprior[[i]]	<- hyperprior[[3]]
        }
      }
      # checking
      if(model.evo %in% "OUM" & ncol(map) != ifelse(root.station, length(do.call(c, init)) - 2, length(do.call(c, init)) - 3)){
        stop("in Ornstein Uhlenbeck with Multiple regimes, number of parameters doesn't correspond to the number of regimes in map")
      }
      if(any(is.na(c(ws$ou.alp, ws$ou.sig, ws$ou.the))) | any(c(ws$ou.alp, ws$ou.sig) <= 0)){
        stop("window.size is not valid")
      }
      if(any(is.na(c(init$ou.alp, init$ou.sig, init$ou.the))) | any(c(init$ou.alp, init$ou.sig) <= 0)){
        stop("initial.values is not valid")
      }
      if(is.finite(hprior[[1]](-1))){
        stop("Hyper prior should not allow sigma <= 0 or alpha <= 0")
      }
    }
    
    return(list(model = model, ws = ws, init = init, prop = prop, map = map, hprior = hprior))
    
  }
  
}



