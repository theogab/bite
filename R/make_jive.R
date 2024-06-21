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
#' trait evolution can be modeled with Ornstein-Uhlenbeck (OU), Brownian Motion (BM) or White Noise (WN) processes. Multiple regimes can be defined for both models and will apply on thetas: c("OU", "theta"), sigmas: c("OU", "sigma") or alphas: c("OU", "alpha") for OU and on sigmas only for WN: c("WN", "sigma") and BM: c("BM", "sigma"). While using the OU model, the user can also relax the stationarity of the root: c("OU", "root") and relax several assumptions at the same time c("OU", "root", "theta") 
#' Species-specific distributions are modeled as multivariate normal distributions. User defined functions of trait evolution can be used in model.priors. The function should be of the form: function(tree, x, pars) and return a loglikelihood value with "tree" being the phylogenetic tree, x being a vector of trait value of size equal to the number of species and ordered as tree$tip.label and pars should be a vector of model parameters (see examples)
#' 
#' parameters used in the different pre-defined models:
#' 
#' White Noise model (WN):
#' \itemize{
#'  \item root: root value
#'  \item sigma_sq: evolutionary rate, n regimes if "sigma" is specified in model.priors
#' }
#'  
#' Brownian Motion model (BM):
#' \itemize{
#'  \item root: root value
#'  \item sigma_sq: evolutionary rate, n regimes if "sigma" is specified in model.priors
#' }
#'
#' Ornstein Uhlenbeck model (OU):
#' \itemize{
#'  \item root: root value. Only used if "root" is specified in model.priors
#'  \item sigma_sq: evolutionary rate, n regimes if "sigma" is specified in model.priors
#'  \item theta: optimal value, n regimes if "theta" is specified in model.priors
#'  \item alpha: strength of selection, n regimes if "alpha" is specified in model.priors
#' }
#' 
#' @param phy phylogenetic tree provided as either a simmap or a phylo object
#' @param traits matrix of traits value for every species of phy (see details)
#' @param map matrix mapping regimes on every edge of phy (see details) 
#' @param model.priors list giving model specification for trait evolution preferably given along with variable names. Supported models are "OU", "BM", "WN". The user can also specify if the assumptions of the model should be relaxed and can also enter a function (see details)				
#' @param scale boolean indicating whether the tree should be scaled to unit length for the model fitting
#' @param nreg integer giving the number of regimes for a Beast analysis. Only evaluated if phy == NULL
#' @param lik.f alternative likelihood function of the form function(pars.lik, traits, counts) to model intraspecific variation (see details)
#' @param init matrix giving initial values for parameters with the variables in rows and the species in columns (see examples)
#' @export
#' @import ape stats
#' @author Theo Gaboriau, Anna Kostikova, Daniele Silvestro and Simon Joly
#' @return A list of functions and tuning parameters (of class "JIVE" and "list") representing the plan of the hierarchical model to parse into \code{\link{mcmc_bite}}.
#' @seealso \code{\link{xml_bite}}, \code{\link{mcmc_bite}} 
#' @examples
#' 
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' 
#' ## JIVE object to run jive with single regimes
#' my.jive <- make_jive(phy = Anolis_tree, traits = Anolis_traits[,-3],
#'  model.priors = list(mean = "BM", logvar= c("OU", "root")))
#'
#' ## JIVE object to run jive with multiple regimes
#' my.jive <- make_jive(Anolis_tree, Anolis_traits[,-3], map = Anolis_map,
#'  model.priors =list(mean = "BM", logvar = c("OU", "theta", "alpha")))
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
#' my.jive <- make_jive(mapped_tree, Anolis_traits[,-3]
#' , model.priors = list(mean = "OU" , logvar = c("OU", "theta")))
#'  
#'  ## Jive object using another model of trait evolution (EB from mvMORPH)
#'  library(mvMORPH)
#'  early_burst <- function(tree, x, pars){
#'   suppressMessages(mvEB(tree, x, method = "inverse", optimization = "fixed", 
#'    echo = FALSE)$llik(pars, root.mle = FALSE))
#'  }
#'  
#'  my.jive <- make_jive(phy = Anolis_tree, traits = Anolis_traits[,-3]
#' , model.priors = list(mean = early_burst , logvar = c("OU", "root")))
#'  initial.values <- c(0.1, 1, 50)
#'  window.size <- c(0.1, 0.2, 1)
#'  proposals <- list("slidingWin", "slidingWin", "slidingWin")
#'  hyperprior <- list(hpfun("Gamma", hp.pars = c(1.1, 5)), hpfun("Gamma", hp.pars = c(3, 5)),
#'                      hpfun("Uniform", hp.pars = c(30, 80)))
#'  names(initial.values) <- names(window.size) <- c("sigma_sq", "beta", "root")
#'  names(proposals) <- names(hyperprior) <- c("sigma_sq", "beta", "root")
#'  my.jive <- control_jive(jive = my.jive, level = "prior", intvar = "mean",
#'   pars = names(initial.values), window.size = window.size,
#'   initial.values = initial.values, proposals = proposals, hyperprior = hyperprior)
#'  
#'  ## Jive object using another model of intraspecific variation (uniform model)
#'  lik_unif <- function(pars.lik, traits, counts){
#'    if(!"mid" %in% names(pars.lik)) stop("'mid' parameter cannot be found in model.priors")
#'    if(!"logrange" %in% names(pars.lik)){
#'     stop("'logrange' parameter cannot be found in model.priors")
#'    }
#'
#'    min.sp <- pars.lik$mid - 1/2*exp(pars.lik$logrange)
#'    max.sp <- pars.lik$mid + 1/2*exp(pars.lik$logrange)
#'    
#'    log.lik.U <- sapply(1:length(traits), function(i){
#'    sum(dunif(traits[[i]], min.sp[i], max.sp[i], log = TRUE))
#'    })
#'    
#'    if (is.na(sum(log.lik.U))) {
#'      return(-Inf)
#'    } else {
#'      return(log.lik.U)
#'    }
#'  }
#'  
#'  init_unif <- sapply(Anolis_tree$tip.label, function(sp){
#'   logrange <- log(diff(range(Anolis_traits[Anolis_traits[,1] == sp, 3])) + 2)
#'   mid <- mean(range(Anolis_traits[Anolis_traits[,1] == sp, 3]))
#'   c(mid = mid, logrange = logrange)
#'  })
#'  
#'  my.jive <- make_jive(phy = Anolis_tree, traits = Anolis_traits[,-2],  
#'  model.priors = list(mid = "BM" , logrange = c("OU", "root")),
#'  lik.f = lik_unif, init = init_unif)
#'  
#' @encoding UTF-8

make_jive <- function(phy = NULL, traits, map = NULL, model.priors = list(mean = "BM", logvar = "OU"), scale = FALSE,
                      nreg = NULL, lik.f = NULL, init = NULL){
  
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
  }, USE.NAMES = TRUE, simplify = FALSE) 
  
  ### validity test ###
  missing <- character(0)
  for(i in 1:length(traits)){
    if(all(is.na(traits[[i]]))){
      missing <- c(missing, names(traits)[[i]])
    } 
  }
  if(length(missing) > 0){
      warning(sprintf("species: %s can not be found in traits. Check matching between species names in phy and traits\nIgnore if you have no data for %s", paste0(missing, collapse = ", "), paste0(missing, collapse = ", ")))
  }
  
  if(is.null(names(model.priors))) names(model.priors) <- sprintf("Prior.mod%s", 1:length(model.priors))
  
  if(!is.null(lik.f)){
    if(is.null(init)){
      warning("Alternative likelihood function: You must add initial values for the likelihood function in init \nRandom initial values set")
      init <- sapply(traits, function(x) rnorm(length(model.priors)))
      rownames(init) <- names(model.priors)
    }
    
  } else {
    if(is.null(init)){
      var.sp <- sapply(traits, var, na.rm = TRUE)
      var.sp[is.na(var.sp)|var.sp == 0] <- var(var.sp, na.rm = TRUE)
      mean.sp <- sapply(traits, mean, na.rm = TRUE)
      mean.sp[is.na(mean.sp)] <- mean(mean.sp, na.rm = TRUE)
      init <- rbind(mean = mean.sp, logvar = log(var.sp))
    }
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
  
  dt <- default_tuning(model.priors = model.priors, phy = jive$data$tree, traits = jive$data$traits, map = jive$data$map, init.lik = init, lik.f = lik.f)
  
  ### Likelihood parameters ###
  jive$lik <- dt$lik
  
  #### priors ####
  jive$priors <- dt$priors
  
  if(!no.tree){
      mat.priors <- list()
      
      for(i in 1:length(model.priors)){
        # Calculate expectation and var/covar matrices #
        mat.priors[[i]] <- try(jive$priors[[i]]$model(x = jive$lik$init[[i]], n = jive$data$n, pars = jive$priors[[i]]$init,
                                            Pi = jive$priors[[i]]$Pi, par.n = 1:ncol(jive$priors[[i]]$Pi), 
                                            data = list(), map = jive$priors[[i]]$map), silent = TRUE)
        
        if(any(grepl("Error", mat.priors[[i]]))){
          warning(sprintf("Initial values for %s prior return an error: %s\nConsider changing initial values using control_jive()",
                  names(jive$priors)[i], mat.priors[[i]][1]))
        } else {
          jive$priors[[i]]$data <- mat.priors[[i]]$data
          jive$priors[[i]]$value <- mat.priors[[i]]$loglik
        }
    }
  }
  
  for(i in 1:length(model.priors)){
    message(names(jive$priors)[i], " prior model: ",jive$priors[[i]]$name)
  }
  
  
  ## Checks
  check_tuning(jive)
  
  #### Prepare headers of log file ####
  
  jive$header <- c("iter", "posterior", "log.lik", sprintf("prior.%s", names(model.priors)),
                   unlist(lapply(1:length(model.priors), function(i){
                     c(paste(names(model.priors)[i], names(jive$priors[[i]]$init), sep ="."),
                     paste(names(model.priors)[i], names(jive$data$traits), sep="_"))
                   })), "acc", "temperature")			
  
  class(jive) <- c("JIVE", "list")
  return(jive)
  
}




