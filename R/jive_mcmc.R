#' @title Jive MCMC algorithm
#' @description Implements Markov chain Monte Carlo sampling for trait evolutionary models with intraspecific data
#' 
#' @details This function runs MCMC sampling on jive object \code{\link{make_jive}}. The jive object contains 
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
#'
#'  ## Load test data
#' data(traitsOU1)
#' data(treeOU1)
#' data(regimesOU1)
#' 
#' ## Run a simple MCMC chain
#' my.jive <- make_jive(treeOU1, traitsOU1,  model.var="OU", model.mean="BM")
#' mcmc_jive(my.jive, log.file="my.jive_MCMC.log", sampling.freq=10, print.freq=100, ngen=5000) 
#'
#' ## Run an MCMC chain with thermodynamic integration
#' mcmc_jive(my.jive, log.file="my.jive_MCMC.log", ncat=10, sampling.freq=10, print.freq=100, ngen=5000, burnin=500) 


# MCMC MAIN ALGORITHM
mcmc_jive <- function(jive, log.file = "my_jive_mcmc.log", sampling.freq = 1000, print.freq = 1000, 
				ncat = 1, beta.param = 0.3, ngen = 5000000, burnin = 0, update.freq = c(0.35,0.2,0.45))
{
	
  ## checking
  if (length(update.freq) != 3 && !is.null(update.freq)) {
    stop("update.freq must contain 3 elements" )
  }
  
  
  ## define the chain length for each category
  it <- ngen/ncat
  
  
  ## get the heating parameter for a chain - scaling classes
  if (ncat > 1) {
    beta.class <- heat_par(ncat, beta.param) 
  } else {
    beta.class <- 1
  }
  
  
  ## initial conditions
  cat("setting initial conditions\n")
  # likelihood level
  m.sp0 <- jive$lik$init$m.sp
  v.sp0 <- jive$lik$init$v.sp
  lik0 <- jive$lik$model(m.sp0, v.sp0, jive$data$traits, jive$data$counts)
  
  # prior mean level
  pars.m0 <- do.call(c, jive$prior.mean$init)
  prior.mean0 <- jive$prior.mean$model(pars.m0, m.sp0, jive$data$tree, jive$data$map)
  hprior.mean0 <- mapply(do.call, jive$prior.mean$hprior, lapply(pars.m0, list))
  
  # prior var level
  pars.v0 <- do.call(c, jive$prior.var$init)
  prior.var0 <- jive$prior.var$model(pars.v0, log(v.sp0), jive$data$tree, jive$data$map)
  hprior.var0 <- mapply(do.call, jive$prior.var$hprior, lapply(pars.v0, list))
  
  # mcmc parameters
  cat("generation\tposterior\n")
  cat(sprintf("%s\t", jive$header), "\n", append = FALSE, file = log.file)
  it.beta <- 1
  bet <- beta.class[it.beta]
  if(ncat > 1) cat("beta = ", bet, "\n")
  update.freq <- cumsum(update.freq/sum(update.freq))
  proposals <- c(0,0,0) # 1st number: update means and variance; 2nd: update mean priors, 3rd: update variance prior
  proposals.accepted <- c(0,0,0)
  
  # posterior level
  post0 <- (sum(lik0) + (prior.mean0 * bet) + (prior.var0 * bet) + sum(c(hprior.mean0, hprior.var0)))
  
  
  ## Start iterations
	for (i in 1:(it*ncat)) {
		
	  # test whether we update lik (r[1]), prior_mean (r[2]) or prior_var (r[3])
	  r <- runif(1) <= update.freq
	  r[2] <- r[2] & !r[1]
	  r[3] <- r[3] & !r[2] & !r[1]

		proposals[r] <- proposals[r] + 1
		
		if (r[1]) # update m.sp, v.sp, calculate new lik
		{ 
			ind <- sample(1:length(m.sp0), 5, replace = F) # 5 is just for now number of species updated
			u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
			m.sp1 <- m.sp0
			v.sp1 <- v.sp0
			lik1 <- lik0
			
			if (runif(1) < 0.5) # change the vector of mean
			{
			  tmp <- jive$lik$prop$m.sp(i = m.sp0[ind], d = jive$lik$ws$m.sp[ind], u)
				m.sp0[ind] <- tmp$v
			} else # change the vector of var
			{
			  tmp <- jive$lik$prop$v.sp(i = v.sp0[ind], d = jive$lik$ws$v.sp[ind], u)
			  v.sp0[ind] <- tmp$v
			}
			
			lik0 <- jive$lik$model(m.sp0, v.sp0, jive$data$traits, jive$data$counts)
		  hasting.ratio <- sum(tmp$lnHastingsRatio)
			
		}

	  if (r[2]) # update parameters of the mean prior, calculate new prior_mean and new hprior_mean
	  {
	    par.n <- sample(1:length(jive$prior.mean$prop), 1) # one parameter at a time for now (update simultaneously sigmas/thetas for multiple regimes?)
			u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
			pars.m1 <- pars.m0 # ancient parameter values are kept
			prior.mean1 <- prior.mean0 # ancient prior_mean value is kept
			hprior.mean1 <- hprior.mean0 # ancient hprior_mean value is kept
			tmp <- jive$prior.mean$prop[[par.n]](i = pars.m0[par.n], d = do.call(c,jive$prior.mean$ws)[par.n], u) #update with proposal function
			pars.m0[par.n] <- tmp$v 
			prior.mean0 <- jive$prior.mean$model(pars.m0, m.sp0, jive$data$tree, jive$data$map) #calculate new values
			hprior.mean0 <- mapply(do.call, jive$prior.mean$hprior, lapply(pars.m0, list))
			hasting.ratio <- tmp$lnHastingsRatio
		} 

		if (r[3]) # update parameters of the var prior, calculate new prior_var and new hprior_var
		{
		  par.n <- sample(1:length(jive$prior.var$prop), 1) # one parameter at a time for now (update simultaneously sigmas/thetas for multiple regimes?)
		  u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
		  pars.v1 <- pars.v0 # ancient parameter values are kept
		  prior.var1 <- prior.var0 # ancient prior_var value is kept
		  hprior.var1 <- hprior.var0 # ancient hprior_var value is kept
		  tmp <- jive$prior.var$prop[[par.n]](i = pars.v0[par.n], d = do.call(c, jive$prior.var$ws)[par.n], u) #update with proposal function
		  pars.v0[par.n] <- tmp$v 
		  prior.var0 <- jive$prior.var$model(pars.v0, log(v.sp0), jive$data$tree, jive$data$map) #calculate new values
		  hprior.var0 <- mapply(do.call, jive$prior.var$hprior, lapply(pars.v0, list))
		  hasting.ratio <- tmp$lnHastingsRatio
		}
		
    # Posterior calculation
		post1 <- post0
		post0 <- (sum(lik0) + (prior.mean0 * bet) + (prior.var0 * bet) + sum(c(hprior.mean0, hprior.var0)))
		
		# acceptance probability (log scale)
		if(any(is.infinite(c(lik0, prior.mean0, prior.var0)))){
		  pr <- -Inf
		} else {
		  pr <- post0 - post1 + hasting.ratio
		}
		
		if (pr >= log(runif(1))) # count acceptance
		{
			proposals.accepted[r] <- proposals.accepted[r] + 1
		} else # cancel changes
		{
		  post0 <- post1
		  
		  if (r[1]){
		    m.sp0 <- m.sp1
		    v.sp0 <- v.sp1
		    lik0 <- lik1
		  }

		  if (r[2]){
		    pars.m0 <- pars.m1
		    prior.mean0 <- prior.mean1
		    hprior.mean0 <- hprior.mean1
		  }

		  if (r[3]){
		    pars.v0 <- pars.v1
		    prior.var0 <- prior.var1
		    hprior.var0 <- hprior.var1
		  }
		}

		# log to file with frequency sampling.freq
		if (i %% sampling.freq == 0 & i >= burnin) {
			 cat(sprintf("%s\t", c(i, post0, sum(lik0),
			 prior.mean0, prior.var0, pars.m0, pars.v0,
			 m.sp0, v.sp0, (sum(proposals.accepted)/i), bet)),
			 "\n", append=TRUE, file=log.file) 
		}
		
		# Print to screen
		if (i %% print.freq == 0) {
			cat(i,'\t',post0,'\n') 
		}
		
		# change beta value if the length of the category is reach 
		if(i%%it == 0 & i < ngen){
		  it.beta = it.beta+1
		  bet <- beta.class[it.beta]
		  cat("beta = ", bet, "\n")
		}
		
	} # end of for 
  
	# Calculate acceptance rate
  names(proposals) <- c("means & variances","prior.means","prior.variances")
  cat("\nEffective proposal frequency\n")
  print(proposals/ngen)
	acceptance.results <- proposals.accepted / proposals
	names(acceptance.results) <- c("means & variances","prior.means","prior.variances")
	cat("\nAcceptance ratios\n")
	print(acceptance.results)
	
}



