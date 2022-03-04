#' @title Control tuning parameters of the jive algorithm
#' @description This function modifies a jive object to tune the jive mcmc algorithm. The output will be different regarding which level of the jive model the user wants to tune ($lik, $priors). This function allows tuning of : initial window size for proposals, starting parameter value, proposal methods, Hyperpriors and update frequencies  
#' @details 
#' If level == "lik" changes will be applied to the likelihood level of the algorithm. intvar is giving the variable on whic the changes will be operated
#' window.size and initial.values must be entered as a vector of length equal to the number of species. proposal must be a character
#' 
#' If level == "prior" changes will be applied to the prior level of the algorithm. intvar is giving the variable on which te change will be operated.
#' window.size and initial.values must be entered as a vector of size equal to the number of parameters or equal to the length of pars. 
#' 
#' Note that if you want to change the tuning at the three levels of the algorithm, you will have to use the control_jive function three times
#' 
#' proposals
#' Has to be one the following : "slidingWin" for Sliding window proposal unconstrained at maximum, "multiplierProposal", for multiplier proposal
#' 
#' Hyperprior
#' list of hyperpriror functions (see \code{\link{hpfun}}). User must provide a list of size equal to the number of parameters or equal to the length of pars
#' 
#' @param jive a jive object obtained from \code{\link{make_jive}}
#' @param level character taken in c("lik", "prior") to specify on which level of the jive model, the control will operate (see details)
#' @param intvar character taken in names(jive$priors) giving the variable to be edited (see details)
#' @param pars vector of character taken in names(jive$priors[[intvar]]$init) giving the names of the hyper parameter to be edited
#' @param window.size initial window size for proposals during the mcmc algorithm. matrix or vector depending on the value of level and nreg (see details)
#' @param initial.values starting parameter values of the mcmc algorithm. matrix or vector depending on the value of level and nreg (see details)
#' @param proposals vector of characters taken in c("slidingWin", "slidingWinAbs", "logSlidingWinAbs","multiplierProposal", "multiplierProposalLakner","logNormal", "absNormal") to control proposal methods during mcmc algorithm (see details)
#' @param hyperprior list of hyperprior functions that can be generated with \code{\link{hpfun}}function. Ignored if level == "lik" (see details)
#' @param update.freq numeric giving the frequency at which parameters should be updated.
#' @export
#' @author Theo Gaboriau
#' @return A JIVE (of class "JIVE" and "list") object to parse into mcmc_bite function (see \code{\link{make_jive}})
#' 
#' @examples
#' 
#' data(Anolis_traits)
#' data(Anolis_tree)
#'  
#' ## Create a jive object
#' my.jive <- make_jive(Anolis_tree, Anolis_traits[,-3],
#' model.priors = list(mean = "BM", logvar= c("OU", "root")))
#' 
#' ## change starting values for the species means
#' my.jive$lik$init #default values
#' new.init <- rep(40,16)
#' my.jive <- control_jive(my.jive, level = "lik", intvar = "mean", initial.values = new.init)
#' my.jive$lik$init #mean initial values changed
#' 
#'  ## change hyperpriors for prior.mean
#'  plot_hp(my.jive) #default values
#'  new.hprior <- list(hpfun("Gamma", hp.pars = c(2,6)), hpfun("Uniform", c(20,80)))
#'  my.jive <- control_jive(my.jive, level = "prior", intvar = "mean", hyperprior = new.hprior)
#'  plot_hp(my.jive) #mean initial values changed
#' @encoding UTF-8

control_jive <- function(jive, level = c("lik", "prior"), intvar = NULL,  pars = NULL, window.size = NULL,
                         initial.values = NULL, proposals = NULL, hyperprior = NULL, update.freq = NULL){
  
  ### Likelihood level ###
  if (level == "lik"){
    
    
    # window size
    if (!is.null(window.size)){
      lab <- names(jive$lik$ws[[intvar]])
      names(window.size) <- lab 
      jive$lik$ws[[intvar]] <- window.size
      
    }
    
    # initial parameter value
    if (!is.null(initial.values)){
      lab <- names(jive$lik$init[[intvar]])
      names(initial.values) <- lab 
      jive$lik$init[[intvar]] <- initial.values
      
    }
    
    # proposals
    if (!is.null(proposals)){
      jive$lik$prop[[intvar]] <- proposal(proposals)
    }
    
    # update frequency
    if (!is.null(update.freq)){
      jive$lik$update.freq <- update.freq
    }
  } else {### Prior level ###
    
    if(is.null(intvar)) stop("intvar must be specified")
    if(!intvar %in% names(jive$priors)) stop(sprintf("variable %s not found in jive object", intvar))
    if(is.null(pars)) pars <- names(jive$priors[[intvar]]$ws) 
    
    # window size
    if(!is.null(window.size)){
      if(length(window.size) != length(pars)) stop("window.size must be of the same size as pars")
      jive$priors[[intvar]]$ws[pars] <- window.size
    }
    
    # initial parameter values
    if(!is.null(initial.values)){
      if(length(initial.values) != length(pars)) stop("initial.values must be of the same size as pars")
      jive$priors[[intvar]]$init[pars] <- initial.values
    }
    
    # proposal functions
    if(!is.null(proposals)){
      if(length(proposals) != length(pars)) stop("proposals must be of the same size as pars")
      for(x in pars){
        jive$priors[[intvar]]$prop[[x]] <- proposal(proposals[[x]])
      }
    }
      
    # hyper priors
    if (!is.null(hyperprior)){
      if(is.null(names(hyperprior))) names(hyperprior) <- pars
      for(x in pars){
        jive$priors[[intvar]]$hprior[[x]] <- hyperprior[[x]]
      }
    }
    
    # update frequency
    if (!is.null(update.freq)){
      jive$priors[[intvar]]$update.freq <- update.freq
    }
  }
  
  ## Calculate new priors
  if(!is.null(jive$data$tree)){
    mat.priors <- list()
    
    for(i in 1:length(jive$priors)){
      # Calculate expectation and var/covar matrices #
      mat.priors[[i]] <- try(jive$priors[[i]]$model(x = jive$lik$init[[i]], n = jive$data$n, pars = jive$priors[[i]]$init,
                                                    Pi = jive$priors[[i]]$Pi, par.n = 1:ncol(jive$priors[[i]]$Pi), 
                                                    data = list(), map = jive$priors[[i]]$map), silent = TRUE)
      
      if(any(grepl("Error", mat.priors[[i]]))){
        warning(sprintf("Initial values for %s prior return an error: %s\nConsider changing initial values",
                        names(jive$priors)[i], mat.priors[[i]][1]))
      } else {
        jive$priors[[i]]$data <- mat.priors[[i]]$data
        jive$priors[[i]]$value <- mat.priors[[i]]$loglik
      }
    }
  }
  
  #### Prepare headers of log file ####
  
  jive$header <- c("iter", "posterior", "log.lik", sprintf("prior.%s", names(jive$priors)),
                   unlist(lapply(1:length(jive$priors), function(i){
                     c(paste(names(jive$priors)[i], names(jive$priors[[i]]$init), sep ="."),
                       paste(names(jive$priors)[i], names(jive$data$traits), sep="_"))
                   })), "acc", "temperature")		
  
  
  check_tuning(jive)
  return(jive)
  
}



