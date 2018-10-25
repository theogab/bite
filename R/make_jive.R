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
#' @param model.var model specification for trait variance evolution. Supported models are c("OU", "BM", "WN", "OUM", "BMM", "WNM")
#' @param model.mean model specification for trait mean evolution. Supported models are c("OU", "BM", "WN", "OUM", "BMM", "WNM")				
#' @param root.station boolean indicating whether the theta_0 should be dropped from the model
#' @param scale boolean indicating whether the tree should be scaled to unit length for the model fitting
#' @param control list to control tuning parameters of the MCMC algorithm (see details)
#' @export
#' @author Anna Kostikova, Simon Joly and Théo Gaboriau
#' @return A jive object to parse into mcmc_jive function
#' @examples
#' library(OUwie)
#' library(phytools)
#' library(MASS)
#' 
#' data(treeOU1)
#' data(traitsOU1)
#' 
#' #VOU and MBM with one regime and without mcmc parameter tuning
#' my.jive <- make_jive(treeOU1, traitsOU1,  model.var="OU", model.mean="BM")
#' 
#' #VOU and MBM with one regime and with mcmc parameter tuning
#' # control argument can be built using control_jive
#' my.control <- list(lik = control_jive("lik", traits = traitsOU1, initial.ws = cbind(apply(traitsOU1, 1, mean, na.rm = T),apply(traitsOU1, 1, sd, na.rm = T))),
#'                    prior.mean = control_jive("prior.mean", model.evo = "BM", traits = traitsOU1, initial.pv = c(1,2)))
#' my.jive <- make_jive(treeOU1, traitsOU1,  model.var="OU", model.mean="BM", control = my.control)                   
#' 
#' data(treeOU2)
#' data(traitsOU2)
#' data(regimesOU2)
#' 
#' #VOUM and MBM with two regimes using simmap tree
#' my.jive <- make_jive(treeOU2, traitsOU2,  model.var="OUM", model.mean="BM")
#' 
#' #VOUM and MBM with two regimes using map
#' my.jive <- make_jive(treeOU2, traitsOU2, map = mapOU2, model.var="OUM", model.mean="BM")

make_jive <- function(phy, traits, map = NULL, model.var="OU", model.mean="BM", root.station = F, scale = F, control = list()){

	### validity test ###
	if (geiger::name.check(phy, traits) != "OK") {
	  stop("Species do not match in tree and traits")
	}

  ### dealing with the map
  if (is.null(map)) {
    if (is.null(phy$mapped.edge)){ # It is assumed there is only one regime
      if (any(c(model.var, model.mean) %in% c("OU","BM","WN")))
        map <- matrix(phy$edge.length)
      else{
        stop("No map provided for multi-regime process")
      }
    } else {
      map <- phy$mapped.edge
    } 
  }
  rownames(map) <- apply(treeOU1$edge, 1, paste, collapse = ",")
  
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
	jive$data$counts 					<- apply(traits, 1, function (x) {sum( !is.na(x) )})
	jive$data$tree   					<- phy
	jive$data$map             <- map
	jive$data$root.station   	<- root.station
	jive$data$scale    	      <- scale
	
	
	
	### Likelihood parameters ###
	jive$lik <- control_jive("lik", traits = traits)
	if(!is.null(control$lik)){ # evaluate control$lik provided by the user and change jive$lik when specified 
	  jive$lik <- mapply(function(a,b){
	    if(is.null(b)) a else b}, jive$lik, control$lik, SIMPLIFY = F)
	}



	#### Models for means ####
	if(model.mean %in% c("OUM", "WNM", "BMM")){
	  nreg.mean <- dim(map)[2]
	} else {
	  nreg.mean <- 1
	}
	
	jive$prior.mean <- control_jive("prior.mean", model.evo = model.mean, traits = traits, nreg = nreg.mean)
	if(!is.null(control$prior.mean)){ # evaluate control$prior.mean provided by the user and change jive$prior.mean when specified 
	  jive$prior.mean <- mapply(function(a,b){
	    if(is.null(b)) a else b}, jive$prior.mean, control$prior.mean, SIMPLIFY = F)
	}
	
	cat("Mean prior model: ",model.mean," [",nreg.mean,"]","\n",sep="")
	
	
	
	#### Models for variance ####
	if(model.var %in% c("OUM", "WNM", "BMM")){
	  nreg.var <- dim(map)[2]
	} else {
	  nreg.var <- 1
	}
	
	jive$prior.var <- control_jive("prior.var", model.evo = model.var, traits = traits, nreg = nreg.var)
	if(!is.null(control$prior.var)){ # evaluate control$prior.var provided by the user and change jive$prior.var when specified 
	  jive$prior.var <- mapply(function(a,b){
	    if(is.null(b)) a else b}, jive$prior.var, control$prior.var)
	}
	
	cat("Variance prior model: ",model.var," [",nreg.var,"]","\n",sep="")
	
	
	
	
	#### Prepare headers of log file ####

	jive$header <- c("Iter", "Posterior", "log.lik", "Prior mean", "Prior var", 
	                 paste("mean.", rep(names(jive$prior.mean$init), lengths(jive$prior.mean$init)),
	                       unlist(lapply(lengths(jive$prior.mean$init), function(x){  
	                          if (x == 1) ""
	                          else if (model.mean %in% c("OU", "OUM")) c("0", seq(1, x-1))
	                            else seq(1, x)
	                        } )), sep =""),
	                 paste("var.", rep(names(jive$prior.var$init), lengths(jive$prior.var$init)),
	                       unlist(lapply(lengths(jive$prior.var$init), function(x){  
	                         if (x == 1) ""
	                         else if (model.var %in% c("OU", "OUM")) c("0", seq(1, x-1))
	                           else seq(1, x)
	                       } )), sep =""),
						       paste(rownames(jive$data$traits), "_m", sep=""),
						       paste(rownames(jive$data$traits), "_v", sep=""),
						       "acc", "temperature")			

	
	return(jive)

}




