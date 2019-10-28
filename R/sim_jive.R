#' @title Simulate JIVE process
#' @description Generate random values of trait mean and variance simulated under a JIVE process along a phylogenetic tree
#' 
#' @details map : the list must be ordered in the same order than phy$edge. Each element represents an edge and contains a vector indicating the time spent under each regime in the branch. The name of the regimes must appear on the map
#' pars : depends on the chosen model.
#' BM and WN : list containing the ancestral state (theta0) and the evolutionary rate (sigma^2)
#' BMM and WNM : list containing the ancestral state (theta0) and the different evolutionary rates corresponding to the regimes (sigma^2i)
#' OU : list containing the ancestral state and optimum (theta0, theta) the evolutionary rate (sigma^2) and the strentgh of selection (alpha)
#' OUM : list containing the ancestral state and optima (theta0, thetai) the evolutionary rates (sigma^2i) and the strentghs of selection (alphai). At least one of thetas, sigmas^2 or alphas must have as many values as regimes.
#' 
#' @param phy Phylogenetic tree 
#' @param map list containing the mapping of regimes over each edge (see details). 
#' @param model.var model specification for the simulation of trait variance evolution. Supported models are c("OU", "BM", "WN", "OUM", "BMM", "WNM")
#' @param model.mean model specification for the simulation of trait mean evolution. Supported models are c("OU", "BM", "WN", "OUM", "BMM", "WNM")					
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
#' # MWNM and VOUM
#' jive_phy <- sim_jive(phy, phy$maps, "OUM", "WNM",  v.pars = list(the = c(10,5,10), 
#' sig2 = 0.1, alp = c(0.2, 0.8)), m.pars = list(the = 0, sig2 = c(0.1, 0.5)))

sim_jive <- function(phy, map = NULL, model.var="OU", model.mean="BM",
                     v.pars = list(the = c(2,1), sig2 = 0.1, alp = 1), m.pars = list(the = 0, sig2 = 0.1),
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

