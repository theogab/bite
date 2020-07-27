#' @import ape

# input: n - number of species, par.n - which parameter has been updated, pars - c(sig1, ..., sigN, the0), tree and map
# does: calculate log-likelihood; 
lik_gen <- function(phy, f.prior){
  eval(parse(text = sprintf("f.prior <- %s", paste0(deparse(f.prior), collapse = "\n"))))  
  foo <- function(x, n, pars, Pi, par.n, data, map){
    data <- NULL
    loglik <- f.prior(phy, x, pars)
    return(list(loglik = loglik, data = data))
  }
}