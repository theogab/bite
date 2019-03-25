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
#'
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#'   
#' my.hp <- hpfun(hpf="Uniform", hp.pars=c(1,2))
#' plot_hp(my.hp)
#' 
#' my.jive <- make_jive(Anolis_tree, Anolis_traits, model.mean="BM", model.var="OU")
#' par(mfrow = c(2,3))
#' plot_hp(my.jive, cex.main = .8)


plot_hp <- function(hpf, col = c("#bfdbf7", "#f49e4c"), border = c("#2e86ab", "#a31621"), bty = "n", ...){
  
  if("JIVE" %in% class(hpf)){
    par(mfrow = c(length(hpf$prior.mean$hprior),length(hpf$prior.var$hprior)))
    for(i in 1:length(hpf$prior.mean$hprior)){
      plot_hyper(hpf$prior.mean$hprior[[i]], col = col[1], border = border[1], xlab = paste("M-",names(hpf$prior.mean$hprior)[i]), bty = bty, ...)
    }
    for(i in 1:length(hpf$prior.var$hprior)){
      plot_hyper(hpf$prior.var$hprior[[i]], col = col[2], border = border[2], xlab = paste("V-", names(hpf$prior.var$hprior)[i]), bty = bty, ...)
    }
  } else {
    plot_hyper(hpf, col = col[1], border = border[1], bty = bty, ...)
  }
}

    