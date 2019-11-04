#' @title Simulate JIVE process
#' @description Generate random values of trait mean and variance simulated under a JIVE process along a phylogenetic tree
#' 
#' @details map : the list must be ordered in the same order than phy$edge. Each element represents an edge and contains a vector indicating the time spent under each regime in the branch. The name of the regimes must appear on the map
#' pars : list containing parameters depending on the chosen model.
#' parameters used in the different models:
#' White Noise model (WN):
#' -theta0: root value, abbreviated the
#' -sigma square: evolutionary rate, abbreviated sig of length n regimes n regimes if "sigma" is specified in model.mean/var
#' 
#' Brownian Motion model (BM):
#' -theta0: root value, abbreviated the
#' -sigma square: evolutionary rate, abbreviated sig of length n regimes n regimes if "sigma" is specified in model.mean/var
#' 
#' Ornstein Uhlenbeck model (OU):
#' -theta0: root value, abbreviated the0 Only used if "root" is specified in model.mean/var
#' -sigma square: evolutionary rate, abbreviated sig of length n regimes if "sigma" is specified in model.mean/var
#' -theta: optimal value, abbreviated the of length n regimes if "theta" is specified in model.mean/var
#' -alpha: strength of selection), abbreviated alp of length n regimes if "alpha" is specified in model.mean/var
#' 
#' @param phy Phylogenetic tree 
#' @param map list containing the mapping of regimes over each edge (see details). 
#' @param model.var model specification for the simulation of trait variance evolution. Supported models are c("OU", "BM", "WN")
#' @param model.mean model specification for the simulation of trait mean evolution. Supported models are c("OU", "BM", "WN")					
#' @param v.pars parameters used for the simulation of trait variance evolution (see details).
#' @param m.pars parameters used for the simulation of trait mean evolution (see details).
#'
#' @import ape
#' @export
#' @author Theo Gaboriau
#' @return returns the simulated mean and variance of the trait for each tip
#' @examples
#'
#' library(phytools)
#' phy <- pbtree(n = 50)
#' Q <- cbind(c(-.002, .002), c(.002, -.002))
#' phy <- sim.history(phy, Q = Q)
#' # MBM and VOU
#' jive_phy <- sim_jive(phy, phy$maps)
#' 
#' # MWN + sigma and VOU + theta + root + alpha
#' jive_phy <- sim_jive(phy, phy$maps, c("WN", "sigma"), c("OU", "theta", "root", "alpha"),  v.pars = list(the = c(10,5,10), 
#' sig2 = 0.1, alp = c(0.2, 0.8)), m.pars = list(the = 0, sig2 = c(0.1, 0.5)))

sim_jive <- function(phy, map = NULL, model.mean="BM", model.var="OU",
                    m.pars = list(the = 0, sig2 = 0.1), v.pars = list(the0 = 2, the = 1, sig2 = 0.1, alp = 1),
                     sampling = list(min = 0, max = 7)){
  
  ntips <- length(phy$tip.label)
  
  ## Simulations for mean
  if(model.mean %in% c("OU", "BM", "WN")) map.mean <- lapply(phy$edge.length, function(x) c('1'= x))
  else map.mean <- map
  mean <- sim_pet(phy, map.mean, model.mean, m.pars, ntips)
  
  ## Simulations for var
  if(model.var %in% c("OU", "BM", "WN")) map.var <- lapply(phy$edge.length, function(x) c('1'= x))
  else map.var <- map
  var <- sim_pet(phy, map.var, model.var, v.pars, ntips)
  
  mv <- cbind(mean, var)
  ## sample individuals in each species
  sp <- apply(mv,1,function(x){
    ind <- rnorm(sampling$max, x[1], sqrt(exp(x[2])))
    ind[sample(1:sampling$max,runif(1, 0, sampling$max-sampling$min))] <- NA
    ind
  })

  out <- cbind(mv,t(sp))
  dimnames(out) <- list(c(phy$tip.label, ntips + 1:phy$Nnode), c("Mean", "Logvariance",1:sampling$max))
  
  return(out)
  
}

