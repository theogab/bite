#' @title Hyper-prior function
#' @description This function creates a hyper-prior density function. 
#' Currently supported density function are Uniform, Gamma and Normal. 
#' The resulting function is used during MCMC \code{\link{jiveMCMC}}
#' to estimate parameters of priors.
#' 
#' @details There are three currently implemented density function: 
#' Uniform, Gamma and Normal. Each of these densities requires two input parameters and hp.pars 
#' must be a vector of two values and cannot be left empty.
#' 
#' 
#' 
#' 
#' @param hpf name of a density function. Supported density functions are: Uniform, Gamma and Normal
#' @param hp.pars a vector of density function parameters
#' @param ... additional parameters that can be passed to a density function
#' @return Hyper-prior density function (function)
#' @export
#' @author Anna Kostikova and Daniele Silvestro
#' @return A hyper-prior density function (function)
#' @examples
#' my.hp <- make.hpfun(hpf="Uniform", hp.pars=c(1,2))


hpfun <-function(hpf="Uniform", hp.pars, ...){
	# Function that makes function of hyper prior
	#
	# Args:
	# 	hpf:  		name of a density function
	#	  hp.pars:	parameters of a density function
	#
	# Returns:
	#	Hyper-prior function (function).

  #uniform
	if (hpf == "Uniform"){
		my.f <- function(x, ...){
			hp <- sum(dunif(x, min=hp.pars[1], max=hp.pars[2], log=TRUE))
			return(hp)
		}
	}
		
  #gamma
	if (hpf == "Gamma"){
		my.f <- function(x, ...){
			hp <- sum(dgamma(x, shape=hp.pars[1], scale=hp.pars[2], log=TRUE))
			return(hp)
		}
	}
		
  #normal
	if (hpf == "Normal"){
		my.f <- function(x, ...){
			hp <- sum(dnorm(x, mean=hp.pars[1], sd=hp.pars[2], log=TRUE))
			return(hp)
		}
	}	
		
  #log gamma
	if (hpf == "Loggamma"){
		my.f <- function(x, ...){
			hp <- sum(dgamma(exp(x), shape=hp.pars[1], scale=hp.pars[2], log=TRUE))
			return(hp)
		}
	}
		
	return(my.f)
  
}
