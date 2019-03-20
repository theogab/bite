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
#' @param ... additional parameters that can be passed to a density function and  \code{\link{par}}
#' @export
#' @author Theo Gaboriau
#' @return plot
#' @examples
#' my.hp <- plot_hp(hpf="Uniform", hp.pars=c(1,2))


plot_hp <- function(hpf, range, col = NULL, border = NULL, ...){
  
  ss <- seq(range[1], range[2], length.out = 1e4)
  d <- exp(sapply(ss, hpf))
  plot(d~ss, type = "n", xlab = "x", ylab = "density", ...)
  if(is.null(border) & is.null(col)){
    col <- "#2e86ab"
  } 
  if(is.null(border)){
    polygon(ss, d, col = col, border = NA)
  } else if(is.null(col)){
    lines(ss, d, col = border)
  } else {
    polygon(ss, d, col = col, border = NA)
    lines(ss, d, col = border)
  }
      
}

    