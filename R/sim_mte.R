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
#'  \item theta0: root value, abbreviated the0 (always length 1)
#'  \item sigma square: evolutionary rate, abbreviated sig
#' }
#' 
#' Brownian Motion model (BM):
#' \itemize{
#'  \item theta0: root value, abbreviated the0 (always length 1)
#'  \item sigma square: evolutionary rate, abbreviated sig
#' }
#' 
#' Ornstein Uhlenbeck model (OU):
#' \itemize{
#'  \item theta0: root value, abbreviated the0 (always length 1)
#'  \item sigma square: evolutionary rate, abbreviated sig
#'  \item theta: optimal value, abbreviated the
#'  \item alpha: strength of selection), abbreviated alp
#' }
#' 
#' Independant Ornstein Uhlenbeck model (IOU):
#' \itemize{
#'  \item theta0: root value, abbreviated ou.the.0. Only used if "root" is specified in model.mean/var
#'  \item sigma square: evolutionary rate, abbreviated ou.sig or ou.sig.1, ou.sig.2, ..., ou.sig.n for n regimes if "sigma" is specified in model.mean/var
#'  \item optimal value, abbreviated ou.the.1 or ou.the.1, ou.the.2, ..., ou.the.n for n regimes if "theta" is specified in model.mean/var
#'  \item stationary variance (alpha/2*sigma_sq with alpha being the strength of selection), abbreviated ou.sv or ou.sv.1
#' }
#' 
#' @param phy Phylogenetic tree 
#' @param map list containing the mapping of regimes over each edge (see details). 
#' @param model model specification for the simulation of trait mean evolution. Supported models are c("OU", "BM", "WN", "IOU")
#' @param pars parameters used for the simulation of trait mean evolution (see details).
#' @param sampling vector of size 2 giving the min and max number of individual per species
#' @param bounds vector of size 2 giving the bounds of the mean
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
#' mte_phy <- sim_mte(phy, phy$maps)
#' 
#' @encoding UTF-8

sim_mte <- function(phy, map = NULL, model = "OU", pars = list(the0 = 2, the = 1, sig = 0.1, alp = 1),
                     sampling = c(1, 7), bounds = c(-Inf, Inf)){
  
  ntips <- length(phy$tip.label)
  
  ## Simulations for mean
  if(all(sapply(pars, length) == 1)) map.mean <- lapply(phy$edge.length, function(x) c('1'= x))
  else map.mean <- map
  mean <- sim_pet(phy, map.mean, model, pars, ntips, bounds)
  return(mean)
  
}

