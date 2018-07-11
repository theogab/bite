#' @title Calculate the marginal likelihood my thermodynamic integration
#' 
#' @description 
#' 
#' @details 
#' 
#' @param file the output file of a jiveMCMC run


marginal.lik <- function(file) {

#	require(dplyr)
	
	# Read the file
	mcmc.table <- read.table(file, header = TRUE)

	# Do some checks
	if (is.null(mcmc.table)) stop(cat("file",file,"does not exist"))
	if (!("postA" %in% colnames(mcmc.table))) stop('No \'postA\' column in the log file')
	if (!("temperature" %in% colnames(mcmc.table))) stop('No \'temperature\' column in the log file')

	# Get mean marginal likelihood for each temperature
	likelihoods <- mcmc.table %>%
		group_by(temperature) %>%
		summarize(Mean = mean(postA, na.rm=TRUE))

	if (dim(likelihoods)[1]<=1) warning('Only one temerature for the thermodynamic integration')

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



marginal.lik.2 <- function(file) {

#	require(dplyr)
	
	# Read the file
	mcmc.table <- read.table(file, header = TRUE)

	# Do some checks
	if (is.null(mcmc.table)) stop(cat("file",file,"does not exist"))
	if (!("postA" %in% colnames(mcmc.table))) stop('No \'postA\' column in the log file')
	if (!("temperature" %in% colnames(mcmc.table))) stop('No \'temperature\' column in the log file')

	# Get mean marginal likelihood for each temperature
	likelihoods <- mcmc.table %>%
		group_by(temperature) %>%
		summarize(Mean = mean(postA, na.rm=TRUE))

	if (dim(likelihoods)[1]<=1) warning('Only one temerature for the thermodynamic integration')

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
