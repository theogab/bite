#' @title Control tuning parameters of the jive algorithm
#' @description This function modifies a jive object to tune the jive mcmc algorithm. The output will be different regarding which level of the jive model the user wants to tune ($lik, $prior.mean, $prior.var). This function allows tuning of : initial window size for proposals, starting parameter value, proposal methods, Hyperpriors and update frequencies  
#' @details 
#' If level == "lik" changes will be applied to the likelihood level of the algorithm
#' window.size and initial.values must be entered as a matrix with 2 columns (respectively mean and variance) and a number of rows equal to the number of species. proposal must be a vector of size 2 (respectively mean and variance)
#' 
#' If level == "prior.mean" or "prior.var" changes will be applied to the prior level of the algorithm
#' window.size and initial.values must be entered as a vector of variable size depending on the chosen evolutionary model. for OU and OUM, the window size and parameter values must be entered in the following order c(alpha, sigma, theta0, theta1, ..., thetaN). for BM, BMM, WN and WNM, the window size and parameter values must be entered in the following order c(sigma1, ..., sigmaN, theta0), proposal must be a vector of size three for OU and OUM c(alpha, sigma, thetas) and of size two for BM, BMM, WN and WNM c(sigmas, theta)
#' 
#' Note that if you want to change the tuning at the three levels of the algorithm, you will have to use the control_jive function three times
#' 
#' proposals
#' Has to be one the following : "slidingWin" for Sliding window proposal unconstrained at maximum, "multiplierProposal", for multiplier proposal
#' 
#' Hyperprior
#' list of hyperpriror functions (see \code{\link{hpfun}}). User must provide a list of size 2 for BM, BMM, WN and WNM (sigmas, theta0) and of size 3 for OU and OUM (alpha, sigma, thetas)
#' 
#' @param jive a jive object obtained from \code{\link{make_jive}}
#' @param level character taken in c("lik", "prior.mean", "prior.var") to specify on which level of the jive model, the control will operate (see details)
#' @param window.size initial window size for proposals during the mcmc algorithm. matrix or vector depending on the value of level and nreg (see details)
#' @param initial.values starting parameter values of the mcmc algorithm. matrix or vector depending on the value of level and nreg (see details)
#' @param proposals vector of characters taken in c("slidingWin", "slidingWinAbs", "logSlidingWinAbs","multiplierProposal", "multiplierProposalLakner","logNormal", "absNormal") to control proposal methods during mcmc algorithm (see details)
#' @param hyperprior list of hyperprior functions that can be generated with \code{\link{hpfun}}function. Ignored if level == "lik" (see details)
#' @param update.freq numeric giving the frequency at which parameters should be updated.
#' @export
#' @author Theo Gaboriau
#' @return A jive object to parse into mcmc_jive function (see \code{\link{make_jive}})
#' 
#' @examples
#' 
#' data(Anolis_traits)
#' data(Anolis_tree)
#'  
#' ## Create a jive object
#' my.jive <- make_jive(Anolis_tree, Anolis_traits,  model.var="OU", model.mean="BM")
#' 
#' ## change starting values for the species mean and variances
#' my.jive$lik$init #default values
#' new.init <- cbind(rep(40,16), rep(20, 16)) 
#' my.jive <- control_jive(my.jive, level = lik, initial.values = new.init)
#' my.jive$lik$init #mean initial values changed
#' 
#' 
#' ## change hyperpriors for prior.mean
#' plot_hp(my.jive) #default values
#' new.hprior <- list(hpfun("Gamma", hp.pars = c(2,6)), hpfun("Uniform", c(20,80)))
#' my.jive <- control_jive(my.jive, level = "prior.mean", hprior = new.hprior)
#' plot_hp(my.jive) #mean initial values changed
#' 

control_jive <- function(jive, level = c("lik", "prior.mean", "prior.var"), window.size = NULL,
                         initial.values = NULL, proposals = NULL, hyperprior = NULL, update.freq = NULL){
  
  ### Likelihood level ###
  if (level == "lik"){
    
    # window size
    if (!is.null(window.size)){
      jive$lik$ws$m.sp <- window.size[,1]
      jive$lik$ws$v.sp <- window.size[,2]
    }
    
    # initial parameter value
    if (!is.null(initial.values)){
      jive$lik$init$m.sp <- initial.values[,1]
      jive$lik$init$v.sp <- initial.values[,2]
    }
    
    # proposals
    if (!is.null(proposals)){
      jive$lik$prop$m.sp <- proposal(proposals[1])
      jive$lik$prop$v.sp	<- proposal(proposals[2])
    }
    
    # update frequency
    if (!is.null(update.freq)){
      jive$lik$update.freq <- update.freq
    }
  }
  
  ### Prior level ###
  if (level == "prior.mean"){
    nreg <- ncol(jive$prior.mean$map)
    # Ornstein-Uhlenbeck #
    if (grepl("OU", jive$prior.mean$name)){
      
      # window size
      if(!is.null(window.size)){
        jive$prior.mean$ws[[1]] <- window.size[1]
        jive$prior.mean$ws[[2]] <- window.size[2]
        jive$prior.mean$ws[[3]] <- window.size[3:(nreg + ifelse(jive$data$root.station, 2, 3))]
      }
      
      # initial parameter values
      if(!is.null(initial.values)){
        jive$prior.mean$init[[1]] <- initial.values[1]
        jive$prior.mean$init[[2]] <- initial.values[2]
        jive$prior.mean$init[[3]] <- initial.values[3:(nreg+ ifelse(jive$data$root.station, 2, 3))]
      }
      
      # proposals
      if (!is.null(proposals)){
        jive$prior.mean$prop[[1]] <- proposal(proposals[1])
        jive$prior.mean$prop[[2]]	<- proposal(proposals[2])
        for(i in 3:(nreg+ ifelse(jive$data$root.station, 2, 3))){
          jive$prior.mean$prop[[i]]	<-  proposal(proposals[3])
        }
      }
      
      # hyper priors
      if (!is.null(hyperprior)){
        jive$prior.mean$hprior[[1]] <- hyperprior[[1]]
        jive$prior.mean$hprior[[2]] <- hyperprior[[2]]
        for(i in 3:(nreg+ ifelse(jive$data$root.station, 2, 3))){
          jive$prior.mean$hprior[[i]]	<- hyperprior[[3]]
        }
      }
    } else {
      
      # window size
      if(!is.null(window.size)){
        jive$prior.mean$ws[[1]] <- window.size[1:nreg]
        jive$prior.mean$ws[[2]] <- window.size[nreg+1]
      }
      
      # initial parameter values
      if(!is.null(initial.values)){
        jive$prior.mean$init[[1]] <- initial.values[1:nreg]
        jive$prior.mean$init[[2]] <- initial.values[nreg]
      }
      
      # proposals
      if (!is.null(proposals)){
        jive$prior.mean$prop <- lapply(1:nreg, proposal, prop = proposals[1]) # sigma(s)
        jive$prior.mean$prop[[nreg+1]]	<- proposal(proposals[2]) # theta
      }
      
      # hyper priors
      if (!is.null(hyperprior)){
        jive$prior.mean$hprior <- lapply(1:nreg, function(x) hyperprior[[1]]) # sigma(s)
        jive$prior.mean$hprior[[nreg+1]] <- hyperprior[[2]] # theta
      }
      
      
    }
    
    # update frequency
    if (!is.null(update.freq)){
      jive$prior.mean$update.freq <- update.freq
    }
  }
  if (level == "prior.var"){
    nreg <- ncol(jive$prior.var$map)
    # Ornstein-Uhlenbeck #
    if (grepl("OU", jive$prior.var$name)){
      
      # window size
      if(!is.null(window.size)){
        jive$prior.var$ws[[1]] <- window.size[1]
        jive$prior.var$ws[[2]] <- window.size[2]
        jive$prior.var$ws[[3]] <- window.size[3:(nreg + ifelse(jive$data$root.station, 2, 3))]
      }
      
      # initial parameter values
      if(!is.null(initial.values)){
        jive$prior.var$init[[1]] <- initial.values[1]
        jive$prior.var$init[[2]] <- initial.values[2]
        jive$prior.var$init[[3]] <- initial.values[3:(nreg+ ifelse(jive$data$root.station, 2, 3))]
      }
      
      # proposals
      if (!is.null(proposals)){
        jive$prior.var$prop[[1]] <- proposal(proposals[1])
        jive$prior.var$prop[[2]]	<- proposal(proposals[2])
        for(i in 3:(nreg+ ifelse(jive$data$root.station, 2, 3))){
          jive$prior.var$prop[[i]]	<-  proposal(proposals[3])
        }
      }
      
      # hyper priors
      if (!is.null(hyperprior)){
        jive$prior.var$hprior[[1]] <- hyperprior[[1]]
        jive$prior.var$hprior[[2]] <- hyperprior[[2]]
        for(i in 3:(nreg+ ifelse(jive$data$root.station, 2, 3))){
          jive$prior.var$hprior[[i]]	<- hyperprior[[3]]
        }
      }
    } else {
      
      # window size
      if(!is.null(window.size)){
        jive$prior.var$ws[[1]] <- window.size[1:nreg]
        jive$prior.var$ws[[2]] <- window.size[nreg+1]
      }
      
      # initial parameter values
      if(!is.null(initial.values)){
        jive$prior.var$init[[1]] <- initial.values[1:nreg]
        jive$prior.var$init[[2]] <- initial.values[nreg]
      }
      
      # proposals
      if (!is.null(proposals)){
        jive$prior.var$prop <- lapply(1:nreg, proposal, prop = proposals[1]) # sigma(s)
        jive$prior.var$prop[[nreg+1]]	<- proposal(proposals[2]) # theta
      }
      
      # hyper priors
      if (!is.null(hyperprior)){
        jive$prior.var$hprior <- lapply(1:nreg, function(x) hyperprior[[1]]) # sigma(s)
        jive$prior.var$hprior[[nreg+1]] <- hyperprior[[2]] # theta
      }
      
      
    }
    
    # update frequency
    if (!is.null(update.freq)){
      jive$prior.var$update.freq <- update.freq
    }
  }
  
  check_tuning(jive)
  return(jive)
  
}



