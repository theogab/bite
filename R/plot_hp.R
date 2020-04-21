#' @title plot Hyper-prior function
#' @description This function plots a hyper-prior density function. 
#' Currently supported density function are Uniform, Gamma, Normal, Loggamma and Lognormal. 
#' The resulting function is used during MCMC \code{\link{mcmc_bite}}
#' to estimate parameters of priors.
#' 
#' @details There are three currently implemented density function: 
#' Uniform, Gamma and Normal. Each of these densities requires two input parameters and hp.pars 
#' must be a vector of two values and cannot be left empty.
#' 
#' @param hpf name of a density function. Supported density functions are: Uniform, Gamma and Normal
#' @param col color of the density area. Can be of size 2 (hpriors for the means, hpriors for the logvars) if a jive object is plotted
#' @param border color of the density curve. Can be of size 2 (hpriors for the means, hpriors for the logvars) if a jive object is plotted
#' @param bty,... additional parameters that can be passed to a density function and  \code{\link[graphics]{par}}
#' @export
#' @author Theo Gaboriau
#' @encoding UTF-8
#' @examples
#'
#'  ## Load test data
#'  data(Anolis_traits)
#'  data(Anolis_tree)
#'    
#'  my.hp <- hpfun(hpf="Uniform", hp.pars=c(1,2))
#'  plot_hp(my.hp)
#'  
#'  my.jive <- make_jive(Anolis_tree, Anolis_traits[,-3], model.priors = list(mean="BM", logvar="OU"))
#'  plot_hp(my.jive, cex.main = .8)



plot_hp <- function(hpf, col = c("#bfdbf7", "#f49e4c"), border = c("#2e86ab", "#a31621"), bty = "n", ...){
  
  if("JIVE" %in% class(hpf)){
    n <- sum(sapply(hpf$priors, function(x) length(x$hprior)))
    nrow <- floor(sqrt(n))
    ncol <- ceiling(sqrt(n))
    ro <- 1
    co <- 1
    oldpar <- par(no.readonly = T)
    on.exit(par(oldpar))
    for(p in 1:length(hpf$priors)){
      for(i in 1:length(hpf$priors[[p]]$hprior)){
        par(fig = c((co-1)/ncol,co/ncol,1-ro/nrow,1-(ro-1)/nrow), new = ifelse(ro == 1 & co == 1, FALSE, TRUE))
        plot_hyper(hpf$priors[[p]]$hprior[[i]], col = col[p], border = border[p],
                   xlab = paste(names(hpf$priors[[p]]$hprior)[i], "[", names(hpf$priors)[p], "]", sep = ""),
                   bty = bty, ...)
        if(co == ncol){
          co <- 1
          ro <- ro + 1
        } else {
          co <- co + 1
        }
      }
    }

  } else {
    plot_hyper(hpf, col = col[1], border = border[1], bty = bty, ...)
  }
}

    