#' @title plot Hyper-prior function
#' @description This function plots a hyper-prior density function. 
#' Currently supported density function are Uniform, Gamma, Normal, Loggamma and Lognormal. 
#' The resulting function is used during MCMC \code{\link{mcmc_jive}}
#' to estimate parameters of priors.
#' 
#' @details There are three currently implemented density function: 
#' Uniform, Gamma and Normal. Each of these densities requires two input parameters and hp.pars 
#' must be a vector of two values and cannot be left empty.
#' 
#' @param hpf name of a density function. Supported density functions are: Uniform, Gamma and Normal
#' @param hp.pars a vector of density function parameters
#' @param ... additional parameters that can be passed to a density function and  \code{\link{par}}
#' @export
#' @author Theo Gaboriau
#' @return plot
#' @examples
#' my.hp <- plot_hp(hpf="Uniform", hp.pars=c(1,2))


plot_hp <- function(hpf = "Uniform", hp.pars = c(1,2), ...){

  ss <- 1e6
  
  #uniform
  if (hpf == "Uniform"){
    my.f <- function(x, ...){
      hp <- runif(x, min=hp.pars[1], max=hp.pars[2])
      return(hp)
    }
    n.pars = paste0(hpf, " density: min = ", hp.pars[1], ", max = ", hp.pars[2])
  }
  
  #gamma
  if (hpf == "Gamma"){
    my.f <- function(x, ...){
      hp <- rgamma(x, shape=hp.pars[1], scale=hp.pars[2])
      return(hp)
    }
    n.pars =paste0(hpf, " density: scale = ", hp.pars[1], ", shape = ", hp.pars[2])
  }
  
  #normal
  if (hpf == "Normal"){
    my.f <- function(x, ...){
      hp <- rnorm(x, mean=hp.pars[1], sd=hp.pars[2])
      return(hp)
    }
    n.pars = paste0(hpf, " density: mean = ", hp.pars[1], ", sd = ", hp.pars[2])
  }	
  
  #log gamma
  if (hpf == "Loggamma"){
    my.f <- function(x, ...){
      hp <- log(rgamma(x, shape=hp.pars[1], scale=hp.pars[2]))
      return(hp)
    }
    n.pars =paste0(hpf, " density: scale = ", hp.pars[1], ", shape = ", hp.pars[2])
  }
  
  #log normal
  if (hpf == "Lognormal"){
    my.f <- function(x, ...){
      hp <- log(rnorm(x, mean=hp.pars[1], sd=hp.pars[2]))
      return(hp)
    }
    n.pars = paste0(hpf, " density: mean = ", hp.pars[1], ", sd = ", hp.pars[2])
  }	
  
  plot(density(my.f(ss)), main = n.pars, xlab = "x")
      
}

    