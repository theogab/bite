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

jiveMake <- function(phy, traits, model.var="OU1", model.mean="BM", root.station = T, scaleHeight = F, control = list()){

	jive <- list()
	
	if (name.check(phy, traits) != "OK") {
	
		stop("Species do not match in tree and traits")
	
	} else {


		### Global variables ###
		jive$data$traits 					<- traits[phy$tip.label,]
		jive$data$counts 					<- apply(traits, 1, function (x) {sum( !is.na(x) )})
		jive$data$tree   					<- phy
		jive$data$vcv    					<- vcv(simmap)  # useless so far
		jive$data$scaleHeight    	<- scaleHeight


		### Likelihood parameters ###
		jive$lik$model 					<- likMultinorm
		jive$lik$mspws 					<- initWinSizeMVN(jive$data$traits)$msp
		jive$lik$sspws 					<- initWinSizeMVN(jive$data$traits)$ssp # normal scale
		jive$lik$mspinit				<- initParamMVN(jive$data$traits)$mspA
		jive$lik$sspinit				<- initParamMVN(jive$data$traits)$sspA # log scale
		jive$lik$prop$msp				<- make.proposal("slidingWin") 
		jive$lik$prop$ssp				<- make.proposal("logSlidingWinAbs") ######### <- HERE


		#### Models for means ####
					
		if (model.mean %in% c("WN", "BM1", "OU1")) { # this is needed to overcome phytools limitation about making simmap obj from a trait with a single regime# td$phy$node.label <- rep("1", n - 1)
			jive$prior_mean$data$map <- matrix(rep(1, ((dim(traits)[1]) * 2) - 2))
		} else {
			if (is.null(map)){
				jive$prior_mean$data$map <- relSim(simmap)$mapped.edge
			} else {
				jive$prior_mean$data$map <- map
			}
		}
				
		if (model.mean %in% c("OUM", "WNM", "BMM")) {
			jive$prior_mean$data$regimes   	<- getStates(simmap,type ="tips")	
			jive$prior_mean$data$nreg   	<- dim(jive$prior_mean$data$map)[2]
		} else {
			jive$prior_mean$data$regimes   	<- "oneregime"
			jive$prior_mean$data$nreg   	<- 1
		}

		cat("Mean prior model: ",model.mean," [",jive$prior_mean$data$nreg,"]","\n",sep="")
		
		if (model.mean == "WN" ) {
			jive$prior_mean$model 			<- likWN
			jive$prior_mean$modelname 		<- model.mean
			jive$prior_mean$init  			<- initParamMWN(jive$data$traits, 1) # vector of two values: the variance and the mean of the species means
			jive$prior_mean$ws	  			<- initWinSizeMWN(jive$data$traits, 1) # vector of two values: both are the standard deviation of the species standard deviations; why?
			jive$prior_mean$hprior$r		<- make.hpfun("Gamma", c(1.1,5)) # sigma
			jive$prior_mean$hprior$m		<- make.hpfun("Uniform", c(-20,10)) # anc.mean ######### <- HERE Loggamma
			jive$prior_mean$prop$r			<- make.proposal("multiplierProposalLakner") 
			jive$prior_mean$prop$m			<- make.proposal("slidingWin") ######### <- HERE logsliding window
		}

		if (model.mean == "WNM" ) {
			jive$prior_mean$model 			<- likWNM
			jive$prior_mean$modelname 		<- model.mean
			jive$prior_mean$init  			<- initParamMWN(jive$data$traits, jive$prior_mean$data$nreg) # check
			jive$prior_mean$ws	  			<- initWinSizeMWN(jive$data$traits, jive$prior_mean$data$nreg) # check
			jive$prior_mean$hprior$r		<- make.hpfun("Gamma", c(1.1,5)) # sigma
			for (i in 1:jive$prior_mean$data$nreg){									 # different means
				ti = paste("t",i,sep="")
				jive$prior_mean$hprior[[ti]]<- make.hpfun("Uniform", c(-20,10))
			}
			jive$prior_mean$prop$r			<- make.proposal("multiplierProposalLakner") 
			for (i in 1:jive$prior_mean$data$nreg){									 # theta
				ti = paste("t",i,sep="")
				jive$prior_mean$prop[[ti]]	<- make.proposal("slidingWin") ######### <- HERE
			}
		}
		
		if (model.mean %in% c("BM1","BMM")) {
			if (model.mean == "BM1") {
				jive$prior_mean$modelname 	<- model.mean			
			} else {
				jive$prior_mean$modelname 	<- "BMS" #This is how it is called in OUwie
			}
			jive$prior_mean$root.station	<- root.station			
			jive$prior_mean$init	  		<- initParamMBM(jive$data$traits, jive$prior_mean$data$nreg, jive$prior_mean$modelname, root.station)  # check
			jive$prior_mean$ws	  			<- initWinSizeMBM(jive$data$traits, jive$prior_mean$data$nreg, jive$prior_mean$modelname, root.station)  # check
			# Initialize hyperpriors for the sigmas
			jive$prior_mean$hprior$r		<- make.hpfun("Gamma", c(1.1,5))
			if (jive$prior_mean$data$nreg > 1) {
				for (i in 2:jive$prior_mean$data$nreg){
					ri = paste("r",i,sep="")
					jive$prior_mean$hprior[[ri]]	<- make.hpfun("Gamma", c(1.1,5))
				}			
			}
			# Initialize hyperpriors for the thetas
			if (model.mean == "BM1") {
				jive$prior_mean$hprior$m	<- make.hpfun("Uniform", c(-20,10)) # theta ######### <- HERE Loggamma
			} else if (model.mean == "BMM") {
				if (!root.station) {
					jive$prior_mean$hprior$m	<- make.hpfun("Uniform", c(-20,10)) # theta ######### <- HERE Loggamma
				} else { # root.station == TRUE; one theta per regime
					for (i in 1:jive$prior_mean$data$nreg){									 # theta
						ti = paste("t",i,sep="")
						jive$prior_mean$hprior[[ti]]	<- make.hpfun("Uniform", c(-20,10))  ######### <- HERE Loggamma
					}
				}
			}
			#print(jive$prior_mean$hprior)
			# Initialize hyperpriors for the sigmas
			jive$prior_mean$prop$r		<- make.proposal("multiplierProposalLakner")
			if (jive$prior_mean$data$nreg > 1) {
				for (i in 2:jive$prior_mean$data$nreg){
					ri = paste("r",i,sep="")
					jive$prior_mean$prop[[ri]]	<- make.proposal("multiplierProposalLakner")
				}			
			}
			# Initialize proposals for the thetas
			if (model.mean == "BM1") {
				jive$prior_mean$prop$m	<- make.proposal("slidingWin") 
			} else if (model.mean == "BMM") {
				if (!root.station) {
					jive$prior_mean$prop$m	<- make.proposal("slidingWin")  
				} else { # root.station == TRUE; one theta per regime
					#jive$prior_mean$prop$m	<- make.proposal("slidingWin") 
					for (i in 1:jive$prior_mean$data$nreg){								
						ti = paste("t",i,sep="")
						jive$prior_mean$prop[[ti]]	<- make.proposal("slidingWin")  
					}
				}
			}
			#print(jive$prior_mean$prop)

			# jive$prior_mean$prop$r			<- make.proposal("multiplierProposalLakner") 
			# if (!root.station) jive$prior_mean$prop$m	<- make.proposal("slidingWin") ######### <- HERE
			# for (i in 1:jive$prior_mean$data$nreg){									 # theta
			# 	ti = paste("t",i,sep="")
			# 	jive$prior_mean$prop[[ti]]	<- make.proposal("slidingWin") ######### <- HERE
			# }		
		}

		if (model.mean == "OU1" || model.mean == "OUM") {
			jive$prior_mean$modelname 		<- model.mean
			jive$prior_mean$root.station	<- root.station
			jive$prior_mean$init  			<- initParamMOU(jive$data$traits, jive$prior_mean$data$nreg, root.station)  # check
			jive$prior_mean$ws	  			<- initWinSizeMOU(jive$data$traits, jive$prior_mean$data$nreg, root.station)  # check
			jive$prior_mean$hprior$a		<- make.hpfun("Gamma", c(1.1,5)) # alpha
			jive$prior_mean$hprior$r		<- make.hpfun("Gamma", c(1.1,5)) # sigma
			if (!root.station) jive$prior_mean$hprior$m			<- make.hpfun("Uniform", c(-20,10)) # anc.mean ######### <- HERE Loggamma
			for (i in 1:jive$prior_mean$data$nreg){									 # theta
				ti = paste("t",i,sep="")
				jive$prior_mean$hprior[[ti]]	<- make.hpfun("Uniform", c(-20,10))  ######### <- HERE Loggamma
			}
			
			jive$prior_mean$prop$a			<- make.proposal("multiplierProposalLakner") 
			jive$prior_mean$prop$r			<- make.proposal("multiplierProposalLakner") 
			if (!root.station) jive$prior_mean$prop$m	<- make.proposal("slidingWin") ######### <- HERE
			for (i in 1:jive$prior_mean$data$nreg){									 # theta
				ti = paste("t",i,sep="")
				jive$prior_mean$prop[[ti]]	<- make.proposal("slidingWin") ######### <- HERE
			}		
		}


		### Models for variances

		if (model.var %in% c("WN", "BM1", "OU1")) { # this is needed to overcome phytools limitation about making simmap obj from a trait with a single regime# td$phy$node.label <- rep("1", n - 1)
			jive$prior_var$data$map <- matrix(rep(1, ((dim(traits)[1]) * 2) - 2))

		} else {
			if (is.null(map)){
				jive$prior_var$data$map <- relSim(simmap)$mapped.edge
			} else {
				jive$prior_var$data$map <- map
				
			}
		}

		if (model.var %in% c("OUM", "WNM", "BMM")) {
			jive$prior_var$data$regimes   	<- getStates(simmap,type ="tips")			
			jive$prior_var$data$nreg   		<- dim(jive$prior_var$data$map)[2]
		} else {
			jive$prior_var$data$regimes   	<- "oneregime"
			jive$prior_var$data$nreg   		<- dim(jive$prior_var$data$map)[2]
		}

		cat("Variance prior model: ",model.var," [",jive$prior_var$data$nreg,"]",sep="")

		if (model.var == "WN" ) {
			jive$prior_var$model 			<- likWN
			jive$prior_var$modelname 		<- model.var
			jive$prior_var$init  			<- initParamVWN(jive$data$traits, 1) # check
			jive$prior_var$ws	  			<- initWinSizeVWN(jive$data$traits, 1) # check
			jive$prior_var$hprior$r			<- make.hpfun("Gamma", c(1.1,5)) # sigma
			jive$prior_var$hprior$m			<- make.hpfun("Uniform", c(-20,10)) # anc.mean ######### <- HERE Loggamma
			jive$prior_var$prop$r			<- make.proposal("multiplierProposalLakner") 
			jive$prior_var$prop$m			<- make.proposal("slidingWin") ######### <- HERE logsliding window
			jive$prior_var$header			<- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
												"mbm_sig.sq", "mbm_anc.st",  "vbm_sig.sq", "vbm_anc.st", 
												paste(rownames(jive$data$traits), "_m", sep=""),
												paste(rownames(jive$data$traits), "_v", sep=""),
												"acc", "temperature")
		}

		if (model.var == "WNM" ) {
			jive$prior_var$model 			<- likWNM
			jive$prior_var$modelname 		<- model.var
			jive$prior_var$init  			<- initParamVWN(jive$data$traits, jive$prior_var$data$nreg) # check
			jive$prior_var$ws	  			<- initWinSizeVWN(jive$data$traits, jive$prior_var$data$nreg) # check
			jive$prior_var$hprior$r			<- make.hpfun("Gamma", c(1.1,5)) # sigma
#			jive$prior_var$hprior$m			<- make.hpfun("Uniform", c(-20,10)) # anc.mean ######### <- HERE Loggamma
			for (i in 1:jive$prior_var$data$nreg){									 # different means
				ti = paste("t",i,sep="")
				jive$prior_var$hprior[[ti]]	<- make.hpfun("Uniform", c(-20,10))
			}
			jive$prior_var$prop$r			<- make.proposal("multiplierProposalLakner") 
			# jive$prior_var$prop$m			<- make.proposal("slidingWin") ######### <- HERE logsliding window
			for (i in 1:jive$prior_var$data$nreg){									 # theta
				ti = paste("t",i,sep="")
				jive$prior_var$prop[[ti]]	<- make.proposal("slidingWin") ######### <- HERE
			}
			jive$prior_var$header			<- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
												"mbm_sig.sq", "mbm_anc.st",  "vbm_sig.sq", 
												paste("vwn_theta", seq(1:jive$prior_var$data$nreg),sep=""),
												paste(rownames(jive$data$traits), "_m", sep=""),
												paste(rownames(jive$data$traits), "_v", sep=""),
												"acc", "temperature")
		}
		
		
		if (model.var %in% c("BM1","BMM")) {
			if (model.var == "BM1") {
				jive$prior_var$modelname 	<- model.var			
			} else {
				jive$prior_var$modelname 	<- "BMS" #This is how it is called in OUwie
			}
			jive$prior_var$root.station		<- root.station			
			jive$prior_var$init	  			<- initParamVBM(jive$data$traits, jive$prior_var$data$nreg, jive$prior_var$modelname, root.station)  # check
			jive$prior_var$ws	  			<- initWinSizeMBM(jive$data$traits, jive$prior_var$data$nreg, jive$prior_var$modelname, root.station)  # check
			# Initialize hyperpriors for the sigmas
			jive$prior_var$hprior$r		<- make.hpfun("Gamma", c(1.1,5))
			if (jive$prior_var$data$nreg > 1) {
				for (i in 2:jive$prior_var$data$nreg){
					ri = paste("r",i,sep="")
					jive$prior_var$hprior[[ri]]	<- make.hpfun("Gamma", c(1.1,5))
				}			
			}
			# Initialize hyperpriors for the thetas
			if (model.var == "BM1") {
				jive$prior_var$hprior$m	<- make.hpfun("Uniform", c(-20,10)) # theta ######### <- HERE Loggamma
			} else if (model.var == "BMM") {
				if (!root.station) {
					jive$prior_var$hprior$m	<- make.hpfun("Uniform", c(-20,10)) # theta ######### <- HERE Loggamma
				} else { # root.station == TRUE; one theta per regime
					for (i in 1:jive$prior_var$data$nreg){									 # theta
						ti = paste("t",i,sep="")
						jive$prior_var$hprior[[ti]]	<- make.hpfun("Uniform", c(-20,10))  ######### <- HERE Loggamma
					}
				}
			}
			#print(jive$prior_mean$hprior)
			# Initialize hyperpriors for the sigmas
			jive$prior_var$prop$r		<- make.proposal("multiplierProposalLakner")
			if (jive$prior_var$data$nreg > 1) {
				for (i in 2:jive$prior_var$data$nreg){
					ri = paste("r",i,sep="")
					jive$prior_var$prop[[ri]]	<- make.proposal("multiplierProposalLakner")
				}			
			}
			# Initialize proposals for the thetas
			if (model.var == "BM1") {
				jive$prior_var$prop$m	<- make.proposal("slidingWin") 
			} else if (model.var == "BMM") {
				if (!root.station) {
					jive$prior_var$prop$m	<- make.proposal("slidingWin")  
				} else { # root.station == TRUE; one theta per regime
					#jive$prior_mean$prop$m	<- make.proposal("slidingWin") 
					for (i in 1:jive$prior_var$data$nreg){								
						ti = paste("t",i,sep="")
						jive$prior_var$prop[[ti]]	<- make.proposal("slidingWin")  
					}
				}
			}


			# jive$prior_var$hprior$r			<- make.hpfun("Gamma", c(1.1,5)) # sigma
			# if (!root.station) jive$prior_var$hprior$m		<- make.hpfun("Uniform", c(-20,10)) # anc.mean ######### <- HERE Loggamma
			# for (i in 1:jive$prior_var$data$nreg){									 # theta
			# 	ti = paste("t",i,sep="")
			# 	jive$prior_var$hprior[[ti]]	<- make.hpfun("Uniform", c(-20,10))  ######### <- HERE Loggamma
			# }
			# jive$prior_var$prop$r			<- make.proposal("multiplierProposalLakner") 
			# if (!root.station) jive$prior_var$prop$m	<- make.proposal("slidingWin") ######### <- HERE
			# for (i in 1:jive$prior_var$data$nreg){									 # theta
			# 	ti = paste("t",i,sep="")
			# 	jive$prior_var$prop[[ti]]	<- make.proposal("slidingWin") ######### <- HERE
			# }		
		}

		
		if (model.var == "OU1" || model.var == "OUM") {
			jive$prior_var$modelname 		<- model.var
			jive$prior_var$root.station		<- root.station
			jive$prior_var$init  			<- initParamVOU(jive$data$traits, jive$prior_var$data$nreg, root.station)  # check
			jive$prior_var$ws	  			<- initWinSizeVOU(jive$data$traits, jive$prior_var$data$nreg, root.station)  # check
			jive$prior_var$hprior$a			<- make.hpfun("Gamma", c(1.1,5)) # alpha
			jive$prior_var$hprior$r			<- make.hpfun("Gamma", c(1.1,5)) # sigma
			if (!root.station) jive$prior_var$hprior$m			<- make.hpfun("Uniform", c(-20,10)) # anc.mean ######### <- HERE Loggamma
			for (i in 1:jive$prior_var$data$nreg){									 # theta
				ti = paste("t",i,sep="")
				jive$prior_var$hprior[[ti]]	<- make.hpfun("Uniform", c(-20,10))  ######### <- HERE Loggamma
			}
			
			jive$prior_var$prop$a			<- make.proposal("multiplierProposalLakner") 
			jive$prior_var$prop$r			<- make.proposal("multiplierProposalLakner") 
			if (!root.station) jive$prior_var$prop$m	<- make.proposal("slidingWin") ######### <- HERE
			for (i in 1:jive$prior_var$data$nreg){									 # theta
				ti = paste("t",i,sep="")
				jive$prior_var$prop[[ti]]	<- make.proposal("slidingWin") ######### <- HERE
			}		
		}
		
		### Prepare headers of log file

		if (model.mean == "WN" && model.var == "WN") {
			jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
								"mwn_sig.sq", "mwn_anc.st", "vwn_sig.sq", "vwn_anc.st", 
								paste(rownames(jive$data$traits), "_m", sep=""),
								paste(rownames(jive$data$traits), "_v", sep=""),
								"acc", "temperature")				
		}

		if (model.mean == "WN" && model.var == "WNM") {
			jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
								"mwn_sig.sq", "mwn_anc.st", "vwn_sig.sq",  
								paste("vwn_theta", seq(1:jive$prior_var$data$nreg),sep=""),
								paste(rownames(jive$data$traits), "_m", sep=""),
								paste(rownames(jive$data$traits), "_v", sep=""),
								"acc", "temperature")				
		}

		if (model.mean == "WNM" && model.var == "WN") {
			jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
								"mwn_sig.sq", 
								paste("mwn_theta", seq(1:jive$prior_mean$data$nreg),sep=""), 
								"vwn_sig.sq", "vwn_anc.st", 
								paste(rownames(jive$data$traits), "_m", sep=""),
								paste(rownames(jive$data$traits), "_v", sep=""),
								"acc", "temperature")				
		}

		if (model.mean == "WNM" && model.var == "WNM") {
			jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
								"mwn_sig.sq", 
								paste("mwn_theta", seq(1:jive$prior_mean$data$nreg),sep=""), 
								"vwn_sig.sq",  
								paste("vwn_theta", seq(1:jive$prior_var$data$nreg),sep=""),
								paste(rownames(jive$data$traits), "_m", sep=""),
								paste(rownames(jive$data$traits), "_v", sep=""),
								"acc", "temperature")				
		}

		if ((model.mean == "WN") && (model.var == "BM1")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq", "mwn_anc.st", 
									"vbm_sig.sq", "vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq", "mwn_anc.st", 
									"vbm_sig.sq", "vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "WN") && (model.var == "BMM")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq", "mwn_anc.st", 
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									paste("vbm_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq", "mwn_anc.st", 
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									"vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "WN") && (model.var %in% c("OU1", "OUM"))) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq", "mwn_anc.st", "vou_alpha", "vou_sig.sq", 
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq", "mwn_anc.st", "vou_alpha", "vou_sig.sq", "vou_anc.st",
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "WNM") && (model.var == "BM1")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq",
									paste("mwn_theta", seq(1:jive$prior_mean$data$nreg),sep=""),									 
									"vbm_sig.sq", "vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq",
									"vbm_sig.sq", "vbm_anc.st",
									paste("vbm_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "WNM") && (model.var == "BMM")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq",
									paste("mwn_theta", seq(1:jive$prior_mean$data$nreg),sep=""),									 
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									paste("vbm_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq",
									paste("mwn_theta", seq(1:jive$prior_mean$data$nreg),sep=""),									 
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									"vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "WNM") && (model.var %in% c("OU1", "OUM"))) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq",
									paste("mwn_theta", seq(1:jive$prior_mean$data$nreg),sep=""),									 
									"vou_alpha", "vou_sig.sq", 
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mwn_sig.sq",
									paste("mwn_theta", seq(1:jive$prior_mean$data$nreg),sep=""),									 
									"vou_alpha", "vou_sig.sq", "vou_anc.st",
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "BM1") && (model.var == "WN")) {
			jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
								"mbm_sig.sq", "mbm_anc.st",
								"vwn_sig.sq", "vwn_anc.st", 
								paste(rownames(jive$data$traits), "_m", sep=""),
								paste(rownames(jive$data$traits), "_v", sep=""),
								"acc", "temperature")				
		}
		if ((model.mean == "BMM") && (model.var == "WN")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									paste("mbm_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vwn_sig.sq", "vwn_anc.st", 
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									"mbm_anc.st",
									paste("mbm_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vwn_sig.sq", "vwn_anc.st", 
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "BM1") && (model.var == "WNM")) {
			jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
								"mbm_sig.sq", "mbm_anc.st",
								"vwn_sig.sq",  
								paste("vwn_theta", seq(1:jive$prior_var$data$nreg),sep=""),
								paste(rownames(jive$data$traits), "_m", sep=""),
								paste(rownames(jive$data$traits), "_v", sep=""),
								"acc", "temperature")				
		}

		if ((model.mean == "BMM") && (model.var == "WNM")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									paste("mbm_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vwn_sig.sq",  
									paste("vwn_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									"mbm_anc.st",
									"vwn_sig.sq",  
									paste("vwn_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "BM1") && (model.var == "BM1")) {
			jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
								"mbm_sig.sq", "mbm_anc.st",
								"vbm_sig.sq", "vbm_anc.st",
								paste(rownames(jive$data$traits), "_m", sep=""),
								paste(rownames(jive$data$traits), "_v", sep=""),
								"acc", "temperature")				
		}

		if ((model.mean == "BMM") && (model.var == "BM1")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									paste("mbm_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vbm_sig.sq", "vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									"mbm_anc.st",
									"vbm_sig.sq", "vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "BM1") && (model.var == "BMM")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mbm_sig.sq", "mbm_anc.st",
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									paste("vbm_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mbm_sig.sq", "mbm_anc.st",
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									"vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "BMM") && (model.var == "BMM")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									paste("mbm_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									paste("vbm_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									"mbm_anc.st",
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									"vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "BM1") && (model.var %in% c("OU1", "OUM"))) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mbm_sig.sq", "mbm_anc.st",
									"vou_alpha", "vou_sig.sq",
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mbm_sig.sq", "mbm_anc.st",
									"vou_alpha", "vou_sig.sq", "vou_anc.st",
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean == "BMM") && (model.var %in% c("OU1", "OUM"))) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									paste("mbm_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vou_alpha", "vou_sig.sq",
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									paste("mbm_sig.sq", seq(1:jive$prior_mean$data$nreg),sep=""),
									"mbm_anc.st",
									"vou_alpha", "vou_sig.sq", "vou_anc.st",
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean %in% c("OU1", "OUM")) && (model.var == "BM1")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vbm_sig.sq", "vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq", "mou_anc.st",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vbm_sig.sq", "vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean %in% c("OU1", "OUM")) && (model.var == "BMM")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									paste("vbm_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq", "mou_anc.st",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									paste("vbm_sig.sq", seq(1:jive$prior_var$data$nreg),sep=""),
									"vbm_anc.st",
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean %in% c("OU1", "OUM")) && (model.var %in% c("OU1", "OUM"))) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vou_alpha", "vou_sig.sq", 
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq", "mou_anc.st",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vou_alpha", "vou_sig.sq", "vou_anc.st",
									paste("vou_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean %in% c("OU1", "OUM")) && (model.var == "WNM")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vwn_sig.sq",  
									paste("vwn_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq", "mou_anc.st",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vwn_sig.sq",  
									paste("vwn_theta", seq(1:jive$prior_var$data$nreg),sep=""),
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}

		if ((model.mean %in% c("OU1", "OUM")) && (model.var == "WN")) {
			if (root.station==TRUE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vwn_sig.sq", "vwn_anc.st", 
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
			if (root.station==FALSE) {
				jive$header <- c("real.iter", "postA", "log.lik", "Prior_mean", "Prior_var",  "sumHpriorA", 
									"mou_alpha", "mou_sig.sq", "mou_anc.st",
									paste("mou_theta", seq(1:jive$prior_mean$data$nreg),sep=""),
									"vwn_sig.sq", "vwn_anc.st", 
									paste(rownames(jive$data$traits), "_m", sep=""),
									paste(rownames(jive$data$traits), "_v", sep=""),
									"acc", "temperature")				
			}
		}


		
	}	
	return(jive)

}




