#' @title Make jive object
#' @description This function makes a jive object from a matrix of intraspecific observations
#' and species phylogeny. The obtained jive object can than be used as an input to \code{\link{jiveMCMC}} function
#' Intraspecific observations should be stored as matrix, where lines are vector of observations for each species,
#' with NA for no data. Phylogenetic tree can be either a simmap object (\code{\link{make.simmap}}) or phylo object (\code{\link{as.phylo}})
#' 
#' @details This function creates a jive object needed for \code{\link{jiveMCMC}} function.  
#' Trait values must be stored as a matrix, where lines are vectors of observations for each species, with NA for no data.
#' Rownames are species names. Phylogenetic tree must be provided as either simmap object (for models with multiple regimes)
#' or as a phylo object (for BM or OU1 models). Rownames and tip labels of a phylogenetic tree should match exactly. 
#' There are three models implemeted for estimation of species variances evolution - BM, OU1 and OUM. Evolution of 
#' species means is only implemented with BM model. Species-specific distribution are models as multivariate normal distribution
#' 
#' 
#' @param simmap an object of class "jive" (see details)
#' @param traits name of the output file that will store the log of MCMC chain
#' @param model.var sampling frequency of the MCMC chain (how often chain will be saved into output file
#' @param model.mean printing frequency of the MCMC chain (how often chain will be printed in the R console)					
#' @param model.lik number of classes for thermodynamic integration (see details)
#' @param root.station boolean indicating whether the theta_0 should be dropped from the model (see details)
#' @param scaleHeight boolean indicating whether the tree should be scaled to unit length for the model fitting (see details)
#' @export
#' @author Anna Kostikova and Simon Joly
#' @return An object of class jive
#' @examples
#' library(OUwie)
#' library(phytools)
#' library(MASS)
#' ## number of species we want to simulate
#' n <- 50
#' 
#' ## generate tree with a pure birth model and scale it to the height of 1
#' tree  <- pbtree(b = 1, n = n, scale = 1, nsim = 1, ape = TRUE)
#' treeb <- tree
#' 
#' ## set parameters for OU1 model of species-specific variances
#' sig.sq <- 0.9
#' alpha  <- 0.1
#' theta0 <- 1
#' theta  <- 5
#' 
#' ## set parameters for BM model of specific-specific means
#' sig.sq.bm <- 0.5
#' mu0       <- 350
#' 
#' ## set mean number of observations per species
#' mean.obs <- 20
#' 
#' ## get selective regimes (all 1s because of OU1 model)
#' y <- data.frame(tree$tip.label, rep(1, n))
#' 
#' ## add node labels
#' tree$node.label <- rep("1", n-1)
#' 
#' ## simulate species-specific variances 
#' sigma.val <- abs(OUwie.sim(tree, y, simmap.tree=FALSE, 
#' scaleHeight=TRUE, alpha=rep(alpha,2),
#' sigma.sq=rep(sig.sq,2), theta0=theta0, theta=theta)$X)
#' 
#' ## simulate species-specific means
#' mean.val <- mvrnorm(mu=rep(mu0, length(tree$tip)), Sigma=(sig.sq.bm * vcv(tree)))
#' 
#' ## draw a random number of intraspecific observations for each species
#' spec.obs <- rpois(n, mean.obs)
#' 
#' ## generate a data matrix where rows are species and columns are individual observations	
#' traits <- matrix(rnorm(max(spec.obs) * n, mean=mean.val, sd=sqrt(sigma.val)), 
#' nrow=n, ncol=max(spec.obs))
#' traits <- cbind(as.matrix(max(spec.obs) - spec.obs), traits)
#' 
#' ## function to replace empty cells with NA
#' foo <- function(x){
#'		to <- x[1]
#'		x[1:(to + 1)] <- NA
#'		return(x[-1])
#' }
#' 
#' ## apply to data matrix	
#' traits <- as.matrix(t(apply(traits, 1, foo)))
#' 
#' ## add species names to rownames
#' rownames(traits) <- tree$tip.label
#' my.jive <- jiveMake(treeb, traits,  model.var="OU1", model.mean="BM", model.lik="Multinorm")

#TODO: allow non simmap phylo

make_jive <- function(phy, traits, model.var="OU", model.mean="BM", root.station = T, scale = F, control = list()){

	### validity tests ###
	if (name.check(phy, traits) != "OK") {
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
	jive$lik <- make_control_jive("lik", traits = traits)
	if(!is.null(control$lik)){ # evaluate control$lik provided by the user and change jive$lik when specified 
	  control$lik$model <- NULL # user cannot change the likelihood function
	  jive$lik <- mapply(function(a,b){
	    if(is.null(b)) a else b}, jive$lik, control$lik)
	}



	#### Models for means ####
	if(model.mean %in% c("OUM", "WNM", "BMM")){
	  nreg.mean <- dim(map)[2]
	} else {
	  nreg.mean <- 1
	}
	
	jive$prior.mean <- make_control_jive("prior.mean", model = model.mean, traits = traits, nreg = nreg.mean)
	if(!is.null(control$prior.mean)){ # evaluate control$prior.mean provided by the user and change jive$prior.mean when specified 
	  control$prior.mean$model <- NULL # user cannot change the likelihood function
	  jive$prior.mean <- mapply(function(a,b){
	    if(is.null(b)) a else b}, jive$prior.mean, control$prior.mean)
	}
	
	cat("Mean prior model: ",model.mean," [",nreg,"]","\n",sep="")
	
	
	
	#### Models for variance ####
	if(model.var %in% c("OUM", "WNM", "BMM")){
	  nreg.var <- dim(map)[2]
	} else {
	  nreg.var <- 1
	}
	
	jive$prior.var <- make_control_jive("prior.var", model = model.var, traits = traits, nreg = nreg.var)
	if(!is.null(control$prior.var)){ # evaluate control$prior.var provided by the user and change jive$prior.var when specified 
	  control$prior.var$model <- NULL # user cannot change the likelihood function
	  jive$prior.var <- mapply(function(a,b){
	    if(is.null(b)) a else b}, jive$prior.var, control$prior.var)
	}
	
	cat("Variance prior model: ",model.var," [",nreg,"]","\n",sep="")
	
	
	
	
	#### Prepare headers of log file ####

	jive$header <- c("Iter", "Posterior", "log.lik", "Prior mean", "Prior var", 
	                 paste("mean.", rep(names(jive$prior_mean$pv), lengths(jive$prior_mean$pv)),
	                       unlist(lapply(lengths(jive$prior_mean$pv), function(x){  
	                          if (x == 1) ""
	                          else paste(".regime", seq(1, x), sep ="")
	                        } )), sep =""),
	                 paste("var.", rep(names(jive$prior_var$pv), lengths(jive$prior_var$pv)),
	                       unlist(lapply(lengths(jive$prior_var$pv), function(x){  
	                         if (x == 1) ""
	                         else paste(".regime", seq(1, x), sep ="")
	                       } )), sep =""),
						       paste(rownames(jive$data$traits), "_m", sep=""),
						       paste(rownames(jive$data$traits), "_v", sep=""),
						       "acc", "temperature")			

	
	return(jive)

}




