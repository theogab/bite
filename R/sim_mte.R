#' @title Simulate MTE process
#' @description Generate random values of trait mean simulated under a MTE process along a phylogenetic tree
#' 
#' @details map : the list must be ordered in the same order than phy$edge. Each element represents an edge and contains a vector indicating the time spent under each regime in the branch. The name of the regimes must appear on the map
#' pars : list containing parameters depending on the chosen model. Elements of that lists must be vectors of size 1 or n, with n = number of regimes in the map.
#' Each element of pars must be named with the corresponding parameter abbreviation.
#' Parameters used in the different models:
#' 
#' White Noise model (WN):
#' \itemize{
#'  \item root: root value
#'  \item sigma_sq: evolutionary rate, n regimes if "sigma" is specified in models
#' }
#'  
#' Brownian Motion model (BM):
#' \itemize{
#'  \item root: root value
#'  \item sigma_sq: evolutionary rate, n regimes if "sigma" is specified in models
#' }
#' 
#' Ornstein Uhlenbeck model (OU):
#' \itemize{
#'  \item root: root value. Only used if "root" is specified in models
#'  \item sigma_sq: evolutionary rate, n regimes if "sigma" is specified in models
#'  \item theta: optimal value, n regimes if "theta" is specified in models
#'  \item alpha: strength of selection, n regimes if "alpha" is specified in models
#' }
#' 
#' @param phy Phylogenetic tree 
#' @param map list containing the mapping of regimes over each edge (see details). 
#' @param model model specification for the simulation of trait mean evolution. Supported models are c("OU", "BM", "WN")
#' @param pars parameters used for the simulation of trait mean evolution (see details).
#' @param sampling vector of size 2 giving the min and max number of individual per species
#' @param bounds vector of size 2 giving the bounds of the mean
#'
#' @import ape
#' @export
#' @author Theo Gaboriau
#' @return returns a numeric vector giving the simulated mean value of the trait for each species of the tree.
#' @examples
#'
#' library(phytools)
#' phy <- pbtree(n = 50)
#' Q <- cbind(c(-.002, .002), c(.002, -.002))
#' phy <- sim.history(phy, Q = Q)
#' # MBM and VOU
#' mte_phy <- sim_mte(phy, phy$maps)
#' 
#' @encoding UTF-8

sim_mte <- function(phy, map = NULL, model = "OU", pars = c(root = 2, theta = 1, sigma_sq = 0.1, alpha = 1),
                     sampling = c(1, 7), bounds = c(-Inf, Inf)){
  
  ntips <- length(phy$tip.label)
  
  ## Simulations for mean
  if(all(sapply(pars, length) == 1)) map.mean <- lapply(phy$edge.length, function(x) c('1'= x))
  else map.mean <- map
  mean <- sim_pet(phy, map.mean, model, pars, ntips, bounds)
  return(mean)
  
}

