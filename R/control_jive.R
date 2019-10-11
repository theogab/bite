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
    # Ornstein-Uhlenbeck #
    if (grepl("OU", jive$prior.mean$name)){
      
      # window size
      if(!is.null(window.size)){
        jive$prior.mean$ws[[1]] <- window.size[1:jive$prior.mean$nr[1]]
        jive$prior.mean$ws[[2]] <- window.size[jive$prior.mean$nr[1] + 1:jive$prior.mean$nr[2]]
        jive$prior.mean$ws[[3]] <- window.size[jive$prior.mean$nr[2] + 1:jive$prior.mean$nr[3]]
      }
      
      # initial parameter values
      if(!is.null(initial.values)){
        jive$prior.mean$init[[1]] <- initial.values[1:jive$prior.mean$nr[1]]
        jive$prior.mean$init[[2]] <- initial.values[jive$prior.mean$nr[1] + 1:jive$prior.mean$nr[2]]
        jive$prior.mean$init[[3]] <- initial.values[jive$prior.mean$nr[2] + 1:jive$prior.mean$nr[3]]
      }
      
      # proposals
      if (!is.null(proposals)){
        i <- 1
        while(i <= sum(jive$prior.mean$nr)){
          if (i <= jive$prior.mean$nr[1]) jive$prior.mean$prop[[i]] <- proposal(proposals[1])
          else if (i <= sum(jive$prior.mean$nr[1:2])) jive$prior.mean$prop[[i]] <- proposal(proposals[2])
          else jive$prior.mean$prop[[i]] <- proposal(proposals[3])
          i <- i + 1
        }
      }
      # hyper priors
      if (!is.null(hyperprior)){
        i <- 1
        while(i <= sum(jive$prior.mean$nr)){
          if (i <= jive$prior.mean$nr[1]) jive$prior.mean$hprior[[i]] <- hyperprior[[1]]
          else if (i <= sum(jive$prior.mean$nr[1:2])) jive$prior.mean$hprior[[i]] <- hyperprior[[2]]
          else jive$prior.mean$hprior[[i]] <- hyperprior[[3]]
          i <- i + 1
        }
      }
    } else {
      
      # window size
      if(!is.null(window.size)){
        jive$prior.mean$ws[[1]] <- window.size[1:jive$prior.mean$nr[1]]
        jive$prior.mean$ws[[2]] <- window.size[jive$prior.mean$nr[1]+1]
      }
      
      # initial parameter values
      if(!is.null(initial.values)){
        jive$prior.mean$init[[1]] <- initial.values[1:jive$prior.mean$nr[1]]
        jive$prior.mean$init[[2]] <- initial.values[jive$prior.mean$nr[1]]
      }
      
      # proposals
      if (!is.null(proposals)){
        jive$prior.mean$prop <- lapply(1:jive$prior.mean$nr[1], proposal, prop = proposals[1]) # sigma(s)
        jive$prior.mean$prop[[jive$prior.mean$nr[1]+1]]	<- proposal(proposals[2]) # theta
      }
      
      # hyper priors
      if (!is.null(hyperprior)){
        jive$prior.mean$hprior <- lapply(1:jive$prior.mean$nr[1], function(x) hyperprior[[1]]) # sigma(s)
        jive$prior.mean$hprior[[jive$prior.mean$nr[1]+1]] <- hyperprior[[2]] # theta
      }
      
      
    }
    
    # update frequency
    if (!is.null(update.freq)){
      jive$prior.mean$update.freq <- update.freq
    }
  }
  if (level == "prior.var"){
    # Ornstein-Uhlenbeck #
    if (grepl("OU", jive$prior.var$name)){
      
      # window size
      if(!is.null(window.size)){
        jive$prior.var$ws[[1]] <- window.size[1:jive$prior.var$nr[1]]
        jive$prior.var$ws[[2]] <- window.size[jive$prior.var$nr[1] + 1:jive$prior.var$nr[2]]
        jive$prior.var$ws[[3]] <- window.size[jive$prior.var$nr[2] + 1:jive$prior.var$nr[3]]
      }
      
      # initial parameter values
      if(!is.null(initial.values)){
        jive$prior.var$init[[1]] <- initial.values[1:jive$prior.var$nr[1]]
        jive$prior.var$init[[2]] <- initial.values[jive$prior.var$nr[1] + 1:jive$prior.var$nr[2]]
        jive$prior.var$init[[3]] <- initial.values[jive$prior.var$nr[2] + 1:jive$prior.var$nr[3]]
      }
      
      # proposals
      if (!is.null(proposals)){
        i <- 1
        while(i <= sum(jive$prior.var$nr)){
          if (i <= jive$prior.var$nr[1]) jive$prior.var$prop[[i]] <- proposal(proposals[1])
          else if (i <= sum(jive$prior.var$nr[1:2])) jive$prior.var$prop[[i]] <- proposal(proposals[2])
          else jive$prior.var$prop[[i]] <- proposal(proposals[3])
          i <- i + 1
        }
      }
      # hyper priors
      if (!is.null(hyperprior)){
        i <- 1
        while(i <= sum(jive$prior.var$nr)){
          if (i <= jive$prior.var$nr[1]) jive$prior.var$hprior[[i]] <- hyperprior[[1]]
          else if (i <= sum(jive$prior.var$nr[1:2])) jive$prior.var$hprior[[i]] <- hyperprior[[2]]
          else jive$prior.var$hprior[[i]] <- hyperprior[[3]]
          i <- i + 1
        }
      }
    } else {
      
      # window size
      if(!is.null(window.size)){
        jive$prior.var$ws[[1]] <- window.size[1:jive$prior.var$nr[1]]
        jive$prior.var$ws[[2]] <- window.size[jive$prior.var$nr[1]+1]
      }
      
      # initial parameter values
      if(!is.null(initial.values)){
        jive$prior.var$init[[1]] <- initial.values[1:jive$prior.var$nr[1]]
        jive$prior.var$init[[2]] <- initial.values[jive$prior.var$nr[1]]
      }
      
      # proposals
      if (!is.null(proposals)){
        jive$prior.var$prop <- lapply(1:jive$prior.var$nr[1], proposal, prop = proposals[1]) # sigma(s)
        jive$prior.var$prop[[jive$prior.var$nr[1]+1]]	<- proposal(proposals[2]) # theta
      }
      
      # hyper priors
      if (!is.null(hyperprior)){
        jive$prior.var$hprior <- lapply(1:jive$prior.var$nr[1], function(x) hyperprior[[1]]) # sigma(s)
        jive$prior.var$hprior[[jive$prior.var$nr[1]+1]] <- hyperprior[[2]] # theta
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



