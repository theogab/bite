#' @title Process jive MCMC
#' @description Process the jiveMCMC output log file
#' 
#' @details This function processes the output log file of the \code{\link{jiveMCMC}} function. 
#' It summarizes posterior sample for each variable into summary statistics (
#' e.g. mean, mode, median) and calculates HPD invervals
#' 
#' 
#' 
#' 
#' @param log.file log file recorded by \code{\link{jiveMCMC}} function
#' @param n.spec number of species
#' @param stat which statistics to use to summarize MCMC. By default, set to mode 
#' for prior level parameters and mean for likelihood level parameters. Can also be \code{\link{mean}}, \code{\link{median}}.
#' @param burning how much of burning to disregard		
#' @param probHPD set HPD intervals 
#' @param verbose how much of statistics to return
#' @param ... additional parameters that can be passed to HPDinterval function
#' @author Anna Kostikova
#' @return A list of averaged statistics from MCMC chain for each parameter (list)
#' @export
#' @examples
#' ## Load test data
#' data(traitsOU1)
#' data(treeOU1)
#' my.jive <- jiveMake(treeOU1, traitsOU1,  model.var="OU1", model.mean="BM", model.lik="Multinorm")
#' jiveMCMC(my.jive, log.file="OU1_log.log", sampling.freq=10, print.freq=100, ngen=5000)  
#' my.summary <- jiveProc(log.file="OU1_log.log", n = 50, verbose=FALSE)


jiveProc<-function(log.file = "jive_mcmc_OU1.log", n.spec, stat=jiveMode, burning = 0, probHPD = 0.95, verbose=TRUE, ...){
		# Processing of the MCMC log file
		# User should provide number of species in an analysed tree
		#
		# Args:
		# 	log.file:	log file recorded by jiveMCMC function
		#	n.spec:  		number of species
		#	stat:		which statistics to use to summarize MCMC
		#	burning:	how much of burning to disregard
		#	probHPD:	set HPD intervals 	
		#	verbose:	how much of statistics to return	
		#
		# Returns:
		#	List (averaged statistics from MCMC chain for each parameter).
		
		n <- n.spec
		
		# calculate mode
		jiveMode<-function(x){
			y=hdrcde::hdr(x, all.modes=F)$mode
			return(y)#print(y)
		}
		
		# supress R obsessive desire to plot scientific notation
		options(scipen=999)
		ll 	= list() # empty list to store results
		l	= read.csv(log.file, stringsAsFactors=FALSE, header=T, sep="\t") # read in csv
		g	= dim(l)[1] # number of generations
		p   = dim(l)[2] # number of parameters
		
		if (burning == 0) {
			b = c(p%/%2, p) # keep only half of chain if no burning is given
		} else {
			b = c(burning, g)
		}
		
		# e.g. start of params(7)# end of params (12)# start of spsp params (13)# end of spsp params(112)		
		s = c(7, (p - 3 - n * 2), (p - 2 - n * 2), (p - 3)) 

		
		# apply function to calculate summary
		ll$prior_pars  = apply(l[b[1]:b[2], c(s[1]:s[2])], 2, stat) # mode for prior level pars
		ll$mean_pars    = apply(l[b[1]:b[2], c(s[3]:(s[3] + n - 1))], 2, mean) # mean for lik level pars (mode is too slow)
		ll$var_pars    = apply(l[b[1]:b[2], c((s[3] + n):s[4])], 2, mean)
		names(ll$mean_pars) = sub("_mean","",names(ll$mean_pars))
		names(ll$var_pars) = sub("_var","",names(ll$var_pars))
	
		
		# calculate mcmc and hpd (coda package)
		my.mcmc        = mcmc(data=l[b[1]:b[2], 1:(p-1)])
		my.hpd         = data.frame(HPDinterval(my.mcmc), prob = probHPD, ...)
		
		# get the HPD intervals
		ll$prior_hpd.l = data.frame(my.hpd)[c(s[1]:s[2]), 1]
		ll$prior_hpd.u = data.frame(my.hpd)[c(s[1]:s[2]), 2]
		
		ll$mean_hpd.l   = data.frame(my.hpd)[c(s[3]:(s[3] + n - 1)), 1]
		ll$mean_hpd.u   = data.frame(my.hpd)[c(s[3]:(s[3] + n - 1)), 2]
		
		ll$var_hpd.l   = data.frame(my.hpd)[c((s[3] + n):s[4]), 1]
		ll$var_hpd.u   = data.frame(my.hpd)[c((s[3] + n):s[4]), 2]
		
		names(ll$prior_hpd.l) <- names(ll$prior_pars)
		names(ll$prior_hpd.u) <- names(ll$prior_pars)	

		names(ll$mean_hpd.l) <- names(ll$mean_pars)
		names(ll$mean_hpd.u) <- names(ll$mean_pars)			
		
		names(ll$var_hpd.l) <- names(ll$var_pars)
		names(ll$var_hpd.u) <- names(ll$var_pars)	
			
		if (verbose==TRUE){
			ss <- list()
			ss$prior_pars  <- ll$prior_pars
			ss$prior_hpd.l <- ll$prior_hpd.l 
			ss$prior_hpd.u <- ll$prior_hpd.u 
			return(ss)
		}
		return(ll)

}
