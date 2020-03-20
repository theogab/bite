#' @title Hyper-prior function
#' @description This function creates a hyper-prior density function. 
#' Currently supported density function are Uniform, Gamma and Normal. 
#' The resulting function is used during MCMC \code{\link{mcmc_bite}}
#' to estimate parameters of priors.
#' 
#' @details There are three currently implemented density function: 
#' Uniform, Gamma and Normal. Each of these densities requires two input parameters and hp.pars 
#' must be a vector of two values and cannot be left empty.
#' 
#' @param hpf name of a density function. Supported density functions are: Uniform, Gamma and Normal (abbreviations are not supported)
#' @param hp.pars a vector of density function parameters
#' @param ... additional parameters that can be passed to a density function
#' @import stats
#' @export
#' @author Anna Kostikova and Daniele Silvestro
#' @return A hyper-prior density function (function)
#' @examples
#' my.hp <- hpfun(hpf="Uniform", hp.pars=c(1,2))
#' 
#' @encoding UTF-8


hpfun <-function(hpf="Uniform", hp.pars = c(1,2), ...){
	# Function that makes function of hyper prior
	#
	# Args:
	# 	hpf:  		name of a density function
	#	  hp.pars:	parameters of a density function
	#
	# Returns:
	#	Hyper-prior function (function).
  
  if(!hpf %in% c("Uniform", "Gamma", "Normal","Lognormal")) stop(sprintf("%s distribution is not supported", hpf))
  
  force(hp.pars)
  
  #uniform
	if (hpf == "Uniform"){
		my.f <- function(x, ...){
			hp <- sum(dunif(x, min=hp.pars[1], max=hp.pars[2], log=TRUE))
			return(list(hp, list(hpf, hp.pars)))
		}
	}
		
  #gamma
	if (hpf == "Gamma"){
		my.f <- function(x, ...){
			hp <- sum(dgamma(x, shape=hp.pars[1], scale=hp.pars[2], log=TRUE))
			return(list(hp, list(hpf, hp.pars)))
		}
	}
		
  #normal
	if (hpf == "Normal"){
		my.f <- function(x, ...){
			hp <- sum(dnorm(x, mean=hp.pars[1], sd=hp.pars[2], log=TRUE))
			return(list(hp, list(hpf, hp.pars)))
		}
	}	
		
  #log normal
  if (hpf == "Lognormal"){
    my.f <- function(x, ...){
      hp <- sum(dlnorm(x, meanlog=hp.pars[1], sdlog=hp.pars[2], log=TRUE))
      return(list(hp, list(hpf, hp.pars)))
    }
  }	
		
	return(my.f)
  
}
