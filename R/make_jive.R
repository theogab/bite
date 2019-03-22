#' @title Create a list that can be used as an input to mcmc_jive
#' @description This function creates a jive object from a matrix of intraspecific observations
#' and species phylogeny. The obtained jive object is a list that can than be used as an input to \code{\link{mcmc_jive}} function
#' Intraspecific observations should be stored as matrix, where lines are vector of observations for each species,
#' with NA for no data. Phylogenetic tree can be either a simmap object (\code{\link{make.simmap}}) or phylo object (\code{\link{as.phylo}})
#' 
#' @details This function creates a jive object needed for \code{\link{mcmc_jive}} function.  
#' Trait values must be stored as a matrix, where lines are vectors of observations for each species, with NA for no data. Rownames are species names that should match exactly tip labels of the phylogenetic tree.
#'
#' Phylogenetic tree must be provided as either simmap object or as a phylo object. If the phylogenetic tree is a phylo object but model specification indicates multiple regimes, user must provide a mapping of the regime in map.
#' 
#' map is a matrix giving the mapping of regimes on phy edges. Each row correspond to an edge in phy and each column correspond to a regime. If map is provided the map from the simmap object is ignored.   
#' 
#' variance and mean evolution can be modeled with Ornstein-Uhlenbeck (OU), Brownian Motion (BM) or White Noise (WN) processes. Multiple regimes can be defined for both models and will apply on thetas only for OU (OUM) and on sigmas only for WN (WNM) and BM (BMM)
#' Species-specific distributions are modeled as multivariate normal distributions
#' 
#' control is a list containig tuning parameters acting at different levels of the MCMC algorithm ($lik for likelihood level, $prior.mean for mean prior level and $prior.var for variance prior level). Inside each level ($lik, $prior.mean, $prior.var), the user can modify the default value of initial parameter value ($pv), initial window size ($ws), proposal methods ($prop) for $lik, $prior.mean and $prior.var and hyperpriors ($hprior) for $prior.mean and $prior.var. 
#' The \code{\link{control_jive}} function provides an easy way to modify control parameters (see examples) 
#' 
#' @param phy phylogenetic tree provided as either a simmap or a phylo object
#' @param traits matrix of traits value for every species of phy (see details)
#' @param map matrix mapping regimes on every edge of phy (see details) 
#' @param model.mean model specification for trait mean evolution. Supported models are c("OU", "BM", "WN", "OUM", "BMM", "WNM")				
#' @param model.var model specification for trait variance evolution. Supported models are c("OU", "BM", "WN", "OUM", "BMM", "WNM")
#' @param root.station boolean indicating whether the theta_0 should be dropped from the OU or OUM models 
#' @param scale boolean indicating whether the tree should be scaled to unit length for the model fitting
#' @param control list to control tuning parameters of the MCMC algorithm (see details)
#' @export
#' @import ape 
#' @author Theo Gaboriau, Anna Kostikova and Simon Joly
#' @return A list of functions and tuning parameters to parse into \code{\link{mcmc_jive}} function.
#' @examples
#' 
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' 
#' ## JIVE object to run jive with single regimes
#' my.jive <- make_jive(Anolis_tree, Anolis_traits, model.mean="BM", model.var="OU")
#'
#' ## JIVE object to run jive with multiple regimes
#' my.jive <- make_jive(Anolis_tree, Anolis_traits, map = Anolis_map, model.mean="BM", model.var="OUM")



make_jive <- function(phy, traits, map = NULL, model.mean="BM", model.var="OU", root.station = F, scale = F, control = list()){
  
  ### validity test ###
  if (!all(phy$tip.label %in% rownames(traits))) {
    stop("Species do not match in tree and traits")
  }
  
  ### dealing with the map
  if (is.null(map)) {
    if (is.null(phy$mapped.edge)){ # It is assumed there is only one regime
      map <- matrix(phy$edge.length)
    } else {
      map <- phy$mapped.edge
    } 
  }
  
  if (dim(map)[1] != length(phy$edge.length)) {
    stop("map must provide mapping for every edge of phy")
  }
  if (!all(abs(rowSums(map) - phy$edge.length) < 1e-5)) {
    stop("Mapping does not correspond to edge lengths")
  }
  
  ### scale height ###
  if(scale){
    t.len <- max(branching.times(phy))
    phy$edge.length <- phy$edge.length/t.len
    map <- map/t.len
  }
  
  jive <- list()
  
  ### Global variables ###
  jive$data$traits 					<- traits[phy$tip.label,]
  jive$data$counts 					<- apply(jive$data$traits, 1, function (x) {sum( !is.na(x) )})
  jive$data$n               <- length(phy$tip.label)
  jive$data$tree   					<- phy
  jive$data$vcv             <- vcv(phy)
  jive$data$root.station   	<- root.station
  jive$data$scale    	      <- scale
  jive$data$map             <- map
  
  
  dt <- default_tuning(model.mean = model.mean, model.var = model.var, traits = jive$data$traits, map = jive$data$map, root.station = jive$data$root.station)
  ### Likelihood parameters ###
  jive$lik <- dt$lik
  
  #### Models for means ####
  jive$prior.mean <- dt$prior.mean
  
  # rownames for map
  rownames(jive$prior.mean$map) <- apply(phy$edge, 1, paste, collapse = ",")
  
  # Calculate expectation and var/covar matrices #
  jive$prior.mean$data <- lapply(jive$prior.mean$model(n = jive$data$n, n.p = 1:length(do.call(c,jive$prior.mean$init)),
                                                       pars = do.call(c,jive$prior.mean$init), tree = jive$data$tree,
                                                       map = jive$prior.mean$map, t.vcv = jive$data$vcv, root.station = jive$data$root.station),
                                 function(x) if (x[[1]]) x[[2]])
  
  cat("Mean prior model: ",model.mean," [",ncol(jive$prior.mean$map),"]","\n",sep="")
  
  #### Models for variance ####
  jive$prior.var <- dt$prior.var

  # rownames for map
  rownames(jive$prior.var$map) <- apply(phy$edge, 1, paste, collapse = ",")
  
  # Calculate expectation and var/covar matrices #
  jive$prior.var$data <- lapply(jive$prior.var$model(n = jive$data$n, n.p = 1:length(do.call(c,jive$prior.var$init)),
                                                     pars = do.call(c,jive$prior.var$init), tree = jive$data$tree,
                                                     map = jive$prior.var$map, t.vcv = jive$data$vcv, root.station = jive$data$root.station),
                                function(x) if (x[[1]]) x[[2]])
  
  cat("Variance prior model: ",model.var," [",ncol(jive$prior.var$map),"]","\n",sep="")
  
  ## Checks
  check_tuning(jive)
  
  #### Prepare headers of log file ####
  
  jive$header <- c("iter", "posterior", "log.lik", "prior mean", "prior var", 
                   paste("mean.", rep(names(jive$prior.mean$init), lengths(jive$prior.mean$init)),
                         unlist(lapply(lengths(jive$prior.mean$init), function(x){  
                           if (x == 1) ""
                           else if (model.mean %in% c("OU", "OUM") & !root.station) c("0", seq(1, x-1))
                           else seq(1, x)
                         } )), sep =""),
                   paste("var.", rep(names(jive$prior.var$init), lengths(jive$prior.var$init)),
                         unlist(lapply(lengths(jive$prior.var$init), function(x){  
                           if (x == 1) ""
                           else if (model.var %in% c("OU", "OUM") & !root.station) c("0", seq(1, x-1))
                           else seq(1, x)
                         } )), sep =""),
                   paste(rownames(jive$data$traits), "_m", sep=""),
                   paste(rownames(jive$data$traits), "_v", sep=""),
                   "acc", "temperature")			
  
  class(jive) <- c("JIVE", "list")
  return(jive)
  
}




