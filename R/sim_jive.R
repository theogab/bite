#' @title Simulate JIVE process
#' @description Generate random values of trait mean and variance simulated under a JIVE process along a phylogenetic tree
#' 
#' @details map : the list must be ordered in the same order than phy$edge. Each element represents an edge and contains a vector indicating the time spent under each regime in the branch. The name of the regimes must appear on the map
#' models : trait evolution can be simulated using Ornstein-Uhlenbeck (OU), Brownian Motion (BM) or White Noise (WN) processes. Multiple regimes can be defined for both models and will apply on thetas: c("OU", "theta"), sigmas: c("OU", "sigma") or alphas: c("OU", "alpha") for OU and on sigmas only for WN: c("WN", "sigma") and BM: c("BM", "sigma"). While using the OU model, the user can also relax the stationarity of the root: c("OU", "root") and relax several assumptions at the same time c("OU", "root", "theta") 
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
#' @param map list containing the mapping of regimes over each edge (see details) 
#' @param models list giving model specification for the simulation of trait evolution. Supported models are c("OU", "BM", "WN"). The user can also specify if the assumptions of the model should be relaxed (see details)	
#' @param pars list giving parameters used for the simulation of trait evolution (see details). Name of parameters must be explecitly entered (see details)
#' @param sampling vector of size 2 giving the min and max number of individual per species
#' @param bounds list giving traits bounds
#' @param var.f  alternative function to model intraspecific variation of the form: function(n, pars) with n, number of individuals and pars, parameters of the model 
#' @import ape
#' @export
#' @author Theo Gaboriau
#' @return A list of length two containing a numeric matrix named "evo_traits" giving simulated traits representing intraspecific variation and a data.frame called "obs" containing species name in the fisrt column and individual observation of the trait in the second.
#' @examples
#'
#' library(phytools)
#' phy <- pbtree(n = 50)
#' Q <- cbind(c(-.002, .002), c(.002, -.002))
#' phy <- sim.history(phy, Q = Q)
#' # MBM and VOU
#' jive_phy <- sim_jive(phy = phy, map = phy$maps)
#' 
#' # MWN + sigma and VOU + theta + root + alpha
#' jive_phy <- sim_jive(phy = phy, map = phy$maps, 
#'    models = list(mean= c("WN", "sigma"), logvar = c("OU")),
#'    pars = list(mean = c(root = 0, sigma_sq1 = 0.1, sigma_sq2 = 0.5), 
#'               logvar = c(root = 10, theta1 = 5, theta2 = 10,
#'                          sigma_sq = 0.1, alpha1 = 0.2, alpha2 = 0.8))
#'  )
#' 
#' 
#' # With a different model of intraspecific variation:
#' unif.f <- function(n, pars){
#'  runif(n, pars[1] - exp(pars[2])/2, pars[1] + exp(pars[2])/2)
#' }
#' 
#' unif_phy <- sim_jive(phy = phy, map = phy$maps, models = list(mid=c("BM"), logrange=c("OU")),
#'  pars = list(mid = c(root = 0, sigma_sq = 0.1),
#'              logrange = c(root = 2, theta = 1, sigma_sq = 0.1, alpha = 1)),
#'  var.f = unif.f)
#' 
#' @encoding UTF-8

sim_jive <- function(phy, map = NULL, models = list(mean=c("BM"), logvar=c("OU")),
                     pars = list(mean = c(root = 0, sigma_sq = 0.1), logvar = c(root = 2, theta = 1, sigma_sq = 0.1, alpha = 1)),
                     sampling = c(1, 7), bounds = list(c(-Inf, Inf), c(-Inf, Inf)), var.f = NULL){
  
  ntips <- length(phy$tip.label)
  npars <- length(pars)
  
  ## trait evolution Simulations
  evo_traits <- sapply(1:npars, function(p){
    if(!any(c("theta", "sigma", "alpha") %in% models[[p]])) map.p <- lapply(phy$edge.length, function(x) c('1'= x))
    else map.p <- map
    sim_pet(phy = phy, map = map.p, model = models[[p]], pars = pars[[p]], ntips = ntips, bounds = bounds[[p]])
  })
  
  colnames(evo_traits) <- names(models)
  
  if(is.null(var.f)){
    var.f <- function(n, pars){
      rnorm(n, pars[1], sqrt(exp(pars[2])))
    }
  }
  
  ## sample individuals in each species
  sp <- data.frame(do.call(rbind,lapply(phy$tip.label, function(x){
    if(sampling[1] != sampling[2]) sampling <- seq(sampling[1], sampling[2])
    ind <- var.f(sample(sampling,1), evo_traits[x,])
    cbind(sp = rep(x, length(ind)), ind)
  })))

  out <- list(pars = evo_traits, obs = sp)
  
  return(out)
  
}

