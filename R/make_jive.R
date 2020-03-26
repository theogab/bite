#' @title Create a list that can be used as an input to mcmc_bite
#' @description This function creates a jive object from a matrix of intraspecific observations
#' and species phylogeny. The obtained jive object is a list that can than be used as an input to \code{\link{mcmc_bite}} function
#' Intraspecific observations should be stored as matrix, where lines are vector of observations for each species,
#' with NA for no data. Phylogenetic tree can be either a simmap object (\code{\link{make.simmap}}) or phylo object (\code{\link{as.phylo}})
#' 
#' @details This function creates a jive object needed for \code{\link{mcmc_bite}} function.  
#' Trait values must be stored as a matrix, where lines are vectors of observations for each species, with NA for no data. Rownames are species names that should match exactly tip labels of the phylogenetic tree.
#'
#' Phylogenetic tree must be provided as either simmap object or as a phylo object. If the phylogenetic tree is a phylo object but model specification indicates multiple regimes, user must provide a mapping of the regime in map. If you keep the phy = NULL options the JIVE object can only be parsed to the \code{\link{xml_bite}} function.
#' 
#' map is a matrix giving the mapping of regimes on phy edges. Each row correspond to an edge in phy and each column correspond to a regime. If map is provided the map from the simmap object is ignored.   
#' 
#' variance and mean evolution can be modeled with Ornstein-Uhlenbeck (OU), Brownian Motion (BM) or White Noise (WN) processes. Multiple regimes can be defined for both models and will apply on thetas: model.mean/var = c("OU", "theta"), sigmas: model.mean/var = c("OU", "sigma") or stationary variances: model.mean/var = c("OU", "sv") for OU and on sigmas only for WN: model.mean/var = c("WN", "sigma") and BM: model.mean/var = c("BM", "sigma"). While using the OU model, the user can also relax the stationarity of the root: model.mean/var = c("OU", "root") and relax several assumptions at the same tme model.mean/var = c("OU", "root", "theta") 
#' Species-specific distributions are modeled as multivariate normal distributions
#' 
#' control is a list containig tuning parameters acting at different levels of the MCMC algorithm ($lik for likelihood level, $prior.mean for mean prior level and $prior.var for variance prior level). Inside each level ($lik, $prior.mean, $prior.var), the user can modify the default value of initial parameter value ($pv), initial window size ($ws), proposal methods ($prop) for $lik, $prior.mean and $prior.var and hyperpriors ($hprior) for $prior.mean and $prior.var. 
#' The \code{\link{control_jive}} function provides an easy way to modify control parameters (see examples) 
#' 
#' parameters used in the different models:
#' 
#' White Noise model (WN):
#' \itemize{
#'  \item root: root value
#'  \item sigma^2: evolutionary rate,n regimes if "sigma" is specified in model.mean/var
#' }
#'  
#' Brownian Motion model (BM):
#' \itemize{
#'  \item root: root value
#'  \item sigma^2: evolutionary rate,n regimes if "sigma" is specified in model.mean/var
#' }
#'
#' Ornstein Uhlenbeck model (OU):
#' \itemize{
#'  \item root: root value. Only used if "root" is specified in model.mean/var
#'  \item sigma^2: evolutionary rate, n regimes if "sigma" is specified in model.mean/var
#'  \item theta: optimal value, n regimes if "theta" is specified in model.mean/var
#'  \item alpha: strength of selection, n regimes if "alpha" is specified in model.mean/var
#' }
#' 
#' @param phy phylogenetic tree provided as either a simmap or a phylo object
#' @param traits matrix of traits value for every species of phy (see details)
#' @param map matrix mapping regimes on every edge of phy (see details) 
#' @param model.mean model specification for trait mean evolution. Supported models are "OU", "BM", "WN". The user can also specify if the assumptions of the model should be relaxed (see details)				
#' @param model.var model specification for trait variance evolution. Supported models are "OU", "BM", "WN". The user can also specify if the assumptions of the model should be relaxed (see details)
#' @param scale boolean indicating whether the tree should be scaled to unit length for the model fitting
#' @param control list to control tuning parameters of the MCMC algorithm (see details)
#' @param nreg integer giving the number of regimes for a Beast analysis. Only evaluated if phy == NULL
#' @export
#' @import ape stats
#' @author Theo Gaboriau, Anna Kostikova, Daniele Silvestro and Simon Joly
#' @return A list of functions and tuning parameters to parse into \code{\link{mcmc_bite}} function.
#' @seealso \code{\link{xml_bite}}, \code{\link{mcmc_bite}} 
#' @examples
#' 
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' 
#' ## JIVE object to run jive with single regimes
#' my.jive <- make_jive(phy = Anolis_tree, traits = Anolis_traits,
#'  model.mean="BM", model.var= c("OU", "root"))
#'
#' ## JIVE object to run jive with multiple regimes
#' my.jive <- make_jive(Anolis_tree, Anolis_traits, map = Anolis_map,
#'  model.mean="BM", model.var=c("OU", "theta", "alpha"))
#' 
#' ## JIVE object to run jive from an ancestral state reconstruction (stochastic mapping)
#' # First generate simmap object
#' library(phytools)
#' n= length(Anolis_tree$tip.label)
#' trait = rep(0,n)
#' trait[c(4,3,14,16, 6,5)] = 1
#' names(trait) =  Anolis_tree$tip.label
#' 
#' mapped_tree=make.simmap(Anolis_tree, trait, model='SYM')
#' plotSimmap(mapped_tree)
#' 
#' my.jive <- make_jive(mapped_tree, Anolis_traits, model.mean=c("OU"), model.var=c("OU", "theta"))
#'  
#' @encoding UTF-8

make_jive <- function(phy = NULL, traits, map = NULL, model.mean=c("BM"), model.var=c("OU"),
                      scale = F, control = list(), nreg = NULL){
  
  ### dealing with the tree
  no.tree <- is.null(phy)
  
  if(no.tree){
    phy <- list()
    phy$tip.label <- unique(traits[,1])
  } else {
    if(!is.null(map)){
      rownames(map) <- sprintf("%s,%s", phy$edge[,1], phy$edge[,2])
    }
    phy <- reorder(phy, "postorder")
    if(!is.null(map)){
      map <- map[sprintf("%s,%s", phy$edge[,1], phy$edge[,2]),]
    }
  }

  ### dealing with traits
  traits <- sapply(phy$tip.label, function(sp){
    if(any(traits[,1]  == sp)){
      traits[traits[,1]  == sp,2]
    } else {
      NA  
    }
  }, USE.NAMES = T, simplify = F) 
  
  ### validity test ###
  missing <- character(0)
  for(i in 1:length(traits)){
    if(all(is.na(traits[[i]]))){
      missing <- c(missing, names(traits)[[i]])
    } 
  }
  if(length(missing) > 0){
    warning(sprintf("species: %s can not be found in traits. Check matching between species names in phy and traits", paste0(missing, collapse = ", ")))
  }
  
  if(!no.tree){
    ### dealing with the map
    if (is.null(map)) {
      if (is.null(phy$maps)){ 
        if (is.null(phy$nodelabels)){
          map <- input_to_map(phy) # It is assumed there is only one regime
          reg.names <- NA
        } else {
          map <- input_to_map(phy, ndlabels = phy$nodelabels)
          reg.names <- unique(phy$nodelabels)
        }
      } else {
        map <- input_to_map(phy, simmap = phy$maps)
        reg.names <- colnames(phy$mapped.edge)
      } 
    } else {
      reg.names <- colnames(map)
      map <- input_to_map(phy, map = map)
    }
    
    ### scale height ###
    if(scale){
      t.len <- max(branching.times(phy))
      phy$edge.length <- phy$edge.length/t.len
      map$S <- map$S/t.len
    }
    
  } else {
    map <- input_to_map(phy, nreg = nreg)
    reg.names <- NA
  }
  
  
  jive <- list()
  
  ### Global variables ###
  jive$data$traits          <- traits
  jive$data$counts 					<- sapply(jive$data$traits, function (x) {sum( !is.na(x) )})
  jive$data$n               <- length(phy$tip.label)
  jive$data$map             <- map
  jive$data$reg             <- reg.names
  jive$data$tree   					<- phy
  
  if(!no.tree){
    jive$data$scale    	      <- scale
  }
  
  
  dt <- default_tuning(model.mean = model.mean, model.var = model.var, phy = jive$data$tree, traits = jive$data$traits, map = jive$data$map)
  
  ### Likelihood parameters ###
  jive$lik <- dt$lik
  
  #### Models for means ####
  jive$prior.mean <- dt$prior.mean
  
  #### Models for variance ####
  jive$prior.var <- dt$prior.var
  
  if(!no.tree){
    done <- F
    while(!done){
      # Calculate expectation and var/covar matrices #
      mat.mean <- try(jive$prior.mean$model(x = jive$lik$init$m.sp, n = jive$data$n, pars = jive$prior.mean$init,
                                                        Pi = jive$prior.mean$Pi, par.n = 1:ncol(jive$prior.mean$Pi), 
                                                        data = list(), map = jive$prior.mean$map), silent = T)

      # Calculate expectation and var/covar matrices #
      mat.var <- try(jive$prior.var$model(x = log(jive$lik$init$v.sp), n = jive$data$n, pars = jive$prior.var$init,
                                                        Pi = jive$prior.var$Pi, par.n = 1:ncol(jive$prior.var$Pi), 
                                                        data = list(), map = jive$prior.var$map), silent = T)
      
      if(all(!grepl("Error", mat.mean)) & all(!grepl("Error", mat.var))){
        
        jive$prior.mean$data <- mat.mean$data
        jive$prior.mean$value <- mat.mean$loglik
        jive$prior.var$data <- mat.var$data
        jive$prior.var$value <- mat.var$loglik
        done <- T
        
      } else {
        # new initial conditions
        dt <- default_tuning(model.mean = model.mean, model.var = model.var, phy = jive$data$tree, traits = jive$data$traits, map = jive$data$map)
        
        jive$prior.mean <- dt$prior.mean
        jive$prior.var <- dt$prior.var
      }
    }
  }
  
  cat("Mean prior model: ",model.mean[1]," [",ncol(dt$prior.mean$map$beta),"]","\n",sep="")
  cat("Variance prior model: ",model.var[1]," [",ncol(dt$prior.var$map$beta),"]","\n",sep="")
  
  ## Checks
  check_tuning(jive)
  
  #### Prepare headers of log file ####
  
  jive$header <- c("iter", "posterior", "log.lik", "prior.mean", "prior.var", 
                   paste("mean", names(jive$prior.mean$init), sep ="."),
                   paste("var", names(jive$prior.var$init), sep ="."),
                   paste(names(jive$data$traits), "_m", sep=""),
                   paste(names(jive$data$traits), "_v", sep=""),
                   "acc", "temperature")			
  
  class(jive) <- c("JIVE", "list")
  return(jive)
  
}




