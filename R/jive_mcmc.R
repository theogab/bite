#' @title Jive MCMC
#' @description Implements Markov chain Monte Carlo sampling for trait evolutionary models with intraspecific data
#' 
#' @details This function runs MCMC sampling on jive object \code{\link{jiveMake}}. The jive object contains 
#' both the dataset and set of model to be used in MCMC. This function implements both a conventional MCMC
#' and an MCMC with thermodynamic integration. The latter option is turned off by default and can be changed by
#' setting ncat to values > 1. The recommended ncat for TI is 10. When setting ncat > 1, make sure to specify burning.
#' As a rule of thumb set burning to 1/10 fraction of ngen. 
#' 
#' @param jive an object of class "jive" (see details)
#' @param log.file name of the output file that will store the log of MCMC chain
#' @param sampling.freq sampling frequency of the MCMC chain (how often chain will be saved into output file
#' @param print.freq printing frequency of the MCMC chain (how often chain will be printed in the R console)					
#' @param ncat number of classes for thermodynamic integration (see details)
#' @param beta.param beta value to define classes for thermodynamic integration (see details)
#' @param ngen number of generation in MCMC chain
#' @param burnin a burning phase of MCMC chain (has to be specified for thermodynamic integration)
#' @param update.freq update frequencies for likelihood and prior level parameters
#' @export
#' @author Anna Kostikova, Daniele Silvestro, and Simon Joly
#' @return none
#' @examples
#' ## Load test data
#' data(traitsOU1)
#' data(treeOU1)
#' ## Run a simple MCMC chain
#' my.jive <- jiveMake(treeOU1, traitsOU1,  model.var="OU1", model.mean="BM", model.lik="Multinorm")
#' jiveMCMC(my.jive, log.file="my.jive_MCMC.log", sampling.freq=10, print.freq=100, ngen=5000) 
#'
#' ## Run an MCMC chain with thermodynamic integration
#' jiveMCMC(my.jive, log.file="my.jive_MCMC.log", ncat=10, sampling.freq=10, print.freq=100, ngen=5000, burnin=500) 


# MCMC MAIN ALGORITHM
jiveMCMC <- function(jive, log.file="jive_mcmc.log", sampling.freq=1000, print.freq=1000, 
				ncat=1, beta.param=0.3, ngen=5000000, burnin=0, discr.time = 5000, learn.time = 20000, update.freq=c(0.35,0.2,0.45))
{
		 
	# counter for the total chain length 
	real.iter <- 1 

	# here we define the chain length for each category
	it <- (ngen - burnin)/ncat 
	if (ncat > 1) {
		# here we get the heating parameter for a chain - scaling classes (MIGRATE)
		beta.class <- defClasses(ncat, beta.param) 
	} else {
		beta.class <- 1
	}

	# start looping over scaling classes 
	for (i.beta in 1:length(beta.class)){
		# apply burning only to a first class
		if (i.beta > 1) { 
			burnin <- 0
		}
		#extract temperature
		temperature <- beta.class[i.beta]
		acc <- 0 #acceptance count
		proposals <- c(0,0,0) # 1st number: update means and variance; 2nd: update mean priors, 3rd: update variance prior
		proposals.accepted <- c(0,0,0)
		for (iteration in 1:(it + burnin)) {

			proposal.type <- c(0,0,0) # 1st number: update means and variance; 2nd: update mean priors, 3rd: update variance prior

			hasting.ratio <- 0
			if (real.iter == 1){

			# initialize update frequencies
			update.freq <- initUpdateFreq(update.freq)
			
			# initialize MCMC parameters
			mspA  <- jive$lik$mspinit # initialize means for species
			sspA  <- jive$lik$sspinit # initialize sigma.sq for species, log scaled ** NOT TRUE **
			mbmA  <- jive$prior_mean$init
			bmouA <- jive$prior_var$init # could be aither more realistic values such as means and sds of true data 
			}
		
			msp  <- mspA
			ssp  <- sspA
			mbm  <- mbmA
			bmou <- bmouA # Variance parameters, in that order: alpha, sigma, anc.state (if appropriate), theta1, theta2, ...
			
			r    <- runif(1)
			
			# Update on all parameters level
		
			# update msp, ssp 
			if (r < update.freq[1]) {
				proposal.type[1] <- 1
				 # 5 is just for now number of species updated
				 ind <- sample(1:length(msp), 5, replace=F)
				
				 if (runif(1) < 0.5) {
					  # updating random 5 values from the vector of means 
					  msp[ind] <- jive$lik$prop$msp(i=mspA[ind], d=jive$lik$mspws[ind])$v
					  
					  
				 } else {
					  # updating random 5 values from the vector of sigmas
					  u = runif(5)
					  tmp <- jive$lik$prop$ssp(i=sspA[ind], d=jive$lik$sspws[ind], u)
					  #print(tmp)
					  ssp[ind] <- tmp$v # log scaled in/out, transformed in prop
					  hasting.ratio <- sum(tmp$lnHastingsRatio)
					  #print(hasting.ratio)

				 }
			}
			# update parameters of the mean prior
			else if (r < update.freq[2]) {
				proposal.type[2] <- 1
			
				ind <- sample(1:length(jive$prior_mean$ws), 1)
				u = runif(1)
				tmp <- jive$prior_mean$prop[[ind]](i=mbmA[ind], d=jive$prior_mean$ws[ind], u)
				mbm[ind] <- tmp$v
				hasting.ratio <- tmp$lnHastingsRatio

				 
			} else {# update BMOU parameters 
				proposal.type[3] <- 1
				
				#TODO: update thetas simultaneously to get better posterior distribution
				ind <- sample(1:length(jive$prior_var$ws), 1)
				u = runif(1)
				tmp <- jive$prior_var$prop[[ind]](i=bmouA[ind], d=jive$prior_var$ws[ind], u)
				bmou[ind] <- tmp$v
				hasting.ratio <- tmp$lnHastingsRatio


			}
			
			# Hyperpriors on all parameters level
			# mean of MVN can be negative due to PCA - mean hprior

			hprior_mean <- mapply(do.call, jive$prior_mean$hprior, lapply(mbm, list))
			hprior_var <- mapply(do.call, jive$prior_var$hprior, lapply(bmou, list))
			hprior <- c(hprior_mean, hprior_var)

			if (real.iter > 1) {
				Lik <- LikA
				Prior_mean <- Prior_meanA #priorMVNA
				Prior_var <- Prior_varA #priorBMOUA
			}
			
			# do this for first step always (because we need to have all probabiliiles)     
			if (r < update.freq[1] || real.iter == 1) {
				 

				### Likelihood ###
				Lik <- jive$lik$model(msp, ssp, jive$data$traits, jive$data$counts) # traits, counts
			
				
				### Prior probabilities for mean ###

				if (jive$prior_mean$modelname %in% c("BM1","BMS")) {
					ouwie.data <- data.frame(species=jive$data$tree$tip.label, regime=jive$prior_mean$data$regimes, traits=msp)
					Prior_mean 	<- OUwie.fixed(jive$data$tree, ouwie.data, model=jive$prior_mean$modelname, simmap.tree=TRUE, root.age=NULL, scaleHeight=jive$data$scaleHeight, root.station=jive$prior_mean$root.station, alpha=0.000000001, sigma.sq=mbm[1], theta=mbm[-1], clade=NULL, mserr="none", quiet=TRUE)$loglik							 	
				} else if (jive$prior_mean$modelname %in% c("WN","WNM")) {
					Prior_mean 	<- jive$prior_mean$model(mbm, msp, jive$data$tree, jive$prior_mean$data$regimes, scaleHeight=jive$data$scaleHeight)						 	
				} else { # Model = OU1 or OUM
					ouwie.data <- data.frame(species=jive$data$tree$tip.label, regime=jive$prior_mean$data$regimes, traits=msp)
					Prior_mean 	<- OUwie.fixed(jive$data$tree, ouwie.data, model=jive$prior_mean$modelname, simmap.tree=TRUE, root.age=NULL, scaleHeight=jive$data$scaleHeight, root.station=jive$prior_mean$root.station, alpha=mbm[1], sigma.sq=mbm[2], theta=mbm[-c(1,2)], clade=NULL, mserr="none", quiet=TRUE)$loglik							 	
				}
			
				### Prior probabilities for variance ###

				if (jive$prior_var$modelname %in% c("BM1","BMS")) {
					ouwie.data <- data.frame(species=jive$data$tree$tip.label, regime=jive$prior_var$data$regimes, traits=log(ssp))
					Prior_var 	<- OUwie.fixed(jive$data$tree, ouwie.data, model=jive$prior_var$modelname, simmap.tree=TRUE, root.age=NULL, scaleHeight=jive$data$scaleHeight, root.station=jive$prior_var$root.station, alpha=0.000000001, sigma.sq=bmou[1], theta=bmou[-1], clade=NULL, mserr="none", quiet=TRUE)$loglik							 	
				} else if (jive$prior_var$modelname %in% c("WN","WNM")) {
					Prior_var 	<- jive$prior_var$model(bmou, log(ssp), jive$data$tree, jive$prior_var$data$regimes, scaleHeight=jive$data$scaleHeight)						 	
				} else {
					ouwie.data <- data.frame(species=jive$data$tree$tip.label, regime=jive$prior_var$data$regimes, traits=log(ssp))
					Prior_var 	<- OUwie.fixed(jive$data$tree, ouwie.data, model=jive$prior_var$modelname, simmap.tree=TRUE, root.age=NULL, scaleHeight=jive$data$scaleHeight, root.station=jive$prior_var$root.station, alpha=bmou[1], sigma.sq=bmou[2], theta=bmou[-c(1,2)], clade=NULL, mserr="none", quiet=TRUE)$loglik							 	
				} # Model = OU1 or OUM
				#Prior_var 	<- jive$prior_var$model(bmou, log(ssp), jive$data$tree, jive$data$map)
			
				 										 
			} else if (r<update.freq[2]) {

				### Prior probabilities for mean ###

				if (jive$prior_mean$modelname %in% c("BM1","BMS")) {
					ouwie.data <- data.frame(species=jive$data$tree$tip.label, regime=jive$prior_mean$data$regimes, traits=msp)
					Prior_mean 	<- OUwie.fixed(jive$data$tree, ouwie.data, model=jive$prior_mean$modelname, simmap.tree=TRUE, root.age=NULL, scaleHeight=jive$data$scaleHeight, root.station=jive$prior_mean$root.station, alpha=0.000000001, sigma.sq=mbm[1], theta=mbm[-1], clade=NULL, mserr="none", quiet=TRUE)$loglik							 	
				} else if (jive$prior_mean$modelname %in% c("WN","WNM")) {
					Prior_mean 	<- jive$prior_mean$model(mbm, msp, jive$data$tree, jive$prior_mean$data$regimes, scaleHeight=jive$data$scaleHeight)						 	
				} else { # Model = OU1 or OUM
					ouwie.data <- data.frame(species=jive$data$tree$tip.label, regime=jive$prior_mean$data$regimes, traits=msp)
					Prior_mean 	<- OUwie.fixed(jive$data$tree, ouwie.data, model=jive$prior_mean$modelname, simmap.tree=TRUE, root.age=NULL, scaleHeight=jive$data$scaleHeight, root.station=jive$prior_mean$root.station, alpha=mbm[1], sigma.sq=mbm[2], theta=mbm[-c(1,2)], clade=NULL, mserr="none", quiet=TRUE)$loglik							 	
				}

			} else {

				### Prior probabilities for variance ###

				if (jive$prior_var$modelname %in% c("BM1","BMS")) {
				 	ouwie.data <- data.frame(species=jive$data$tree$tip.label, regime=jive$prior_var$data$regimes, traits=log(ssp))
					Prior_var 	<- OUwie.fixed(jive$data$tree, ouwie.data, model=jive$prior_var$modelname, simmap.tree=TRUE, root.age=NULL, scaleHeight=jive$data$scaleHeight, root.station=jive$prior_var$root.station, alpha=0.000000001, sigma.sq=bmou[1], theta=bmou[-1], clade=NULL, mserr="none", quiet=TRUE)$loglik							 	
				} else if (jive$prior_var$modelname %in% c("WN","WNM")) {
					Prior_var 	<- jive$prior_var$model(bmou, log(ssp), jive$data$tree, jive$prior_var$data$regimes, scaleHeight=jive$data$scaleHeight)						 	
				} else { # Model = OU1 or OUM
					ouwie.data <- data.frame(species=jive$data$tree$tip.label, regime=jive$prior_var$data$regimes, traits=log(ssp))
					Prior_var 	<- OUwie.fixed(jive$data$tree, ouwie.data, model=jive$prior_var$modelname, simmap.tree=TRUE, root.age=NULL, scaleHeight=jive$data$scaleHeight, root.station=jive$prior_var$root.station, alpha=bmou[1], sigma.sq=bmou[2], theta=bmou[-c(1,2)], clade=NULL, mserr="none", quiet=TRUE)$loglik							 	
				}
			}

			# Posterior calculation

			# just for 1 real.iter we need to copy all calculated likelihoods and priors
			if (real.iter == 1) {
				 LikA <- Lik
				 Prior_meanA <- Prior_mean
				 Prior_varA <- Prior_var
				 hpriorA <- hprior
				 postA <- (sum(Lik) + Prior_mean + Prior_var * temperature + sum(hprior))
			}
			post <- (sum(Lik) + Prior_mean + Prior_var * temperature + sum(hprior))

			# acceptance probability
			tryCatch(
			{
			if (post - postA + hasting.ratio >= log(runif(1))){
				 acc = acc + 1
				 proposals.accepted <- proposals.accepted + proposal.type #acceptance rate
				 LikA = Lik
				 Prior_meanA = Prior_mean
				 Prior_varA = Prior_var
				 hpriorA = hprior
				 postA = post
				 mspA = msp
				 sspA = ssp
				 mbmA  = mbm
				 bmouA = bmou
				}
			}
			,error = function(e) NULL
			)

			# Record number of proposal of each type
			proposals <- proposals + proposal.type
			
			# Prepare headers for logfile and screen
			if (real.iter == 1){
				cat("generation",'\t',"posterior",'\n')
				cat(sprintf("%s\t", jive$header), "\n", append=FALSE, file=log.file)
				} 
			
			# log to file with frequency sampling.freq
			if (real.iter %% sampling.freq == 0 & real.iter >= burnin) {
				 cat(sprintf("%s\t", c(real.iter, postA, sum(LikA),
				 Prior_meanA, Prior_varA, sum(hpriorA), mbmA, 
				 bmouA, mspA, sspA, (acc/iteration), temperature)),
				 "\n", append=TRUE, file=log.file) 
				 
			}
			
			# Print to screen
			if (real.iter %% print.freq == 0 & real.iter >= burnin) {
				#cat(real.iter,'\t',postA,'\t',mbmA,'\n')
				cat(real.iter,'\t',postA,'\n') 
			}
			  
			real.iter = real.iter + 1
		} # end for (iteration in 1:(it + burnin))
		# Calculate acceptance rate
		cat("\nUpdate cumulative frequencies\n")
		names(update.freq) <- c("means&variances","prior.means","prior.variances")
		print(update.freq)
		acceptance.results <- proposals.accepted / proposals
		names(acceptance.results) <- c("means&variances","prior.means","prior.variances")
		cat("\nAcceptance ratios\n")
		print(acceptance.results)
	} # end of for (i.beta in 1:length(beta.class))

}



