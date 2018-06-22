#' @title Calculate the marginal likelihood my thermodynamic integration
#' 
#' @description 
#' 
#' @details 
#' 
#' @param file the output file of a jiveMCMC run


marginal.lik <- function(file) {

	require(dplyr)
	
	# Read the file
	mcmc.table <- read.table(file, header = TRUE)

	# mean likelihood for each temperature
	likelihoods <- mcmc.table %>%
		group_by(temperature) %>%
		summarize(Mean = mean(postA, na.rm=TRUE))

	# The marginal likelihood
	lik.ti = 0
	# Calculate the intervals in temperature
	intervals <-  likelihoods$temperature[-1] - likelihoods$temperature[-length(likelihoods$temperature)]
	# Estimate the integral
	for (i in 1:length(intervals)) {
	  lik.ti = lik.ti + ( (likelihoods$Mean[i]+likelihoods$Mean[i+1])/2 * intervals[i])
	}

	return(lik.ti)

}
