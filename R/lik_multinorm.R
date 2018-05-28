# input: xmean = vector of means, xsd = vector of sigmas, traits = matrix of species observations, N_OBSERV vector of observation counts
# does: calculate individual log-likelihoods for each species based on normal distribution
likMultinorm <- function(xmean, xsd, x, counts){#m - mean (horizontal), s - sigma^2 (horizontal), vec - observations for a species
	
	
	log.lik.MN <- -counts/2 * log(2 * pi) - 1/2 * counts * log(xsd) - 1/2 * (apply((x - xmean)^2, 1, sum, na.rm=T)/xsd)
	
	if (is.na(sum(log.lik.MN))) {
			return(-Inf)
	} else {
		return(log.lik.MN)
	}

	
	
	
} # Gaussian density

