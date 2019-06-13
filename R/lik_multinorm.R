# input: m.sp = vector of means, v.sp = vector of sigmas, traits = matrix of species observations, counts = number of observation by species
# does: calculate individual log-likelihoods for each species based on normal distribution
lik_multinorm <- function(m.sp, v.sp, traits, counts){#m - mean (horizontal), s - sigma^2 (horizontal), vec - observations for a species
	
	log.lik.MN <- -counts/2 * log(2 * pi) - 1/2 * counts * log(v.sp) - 1/2 * (sapply(1:length(traits), function(i) sum((traits[[i]] - m.sp[i])^2))/v.sp)
	
	if (is.na(sum(log.lik.MN))) {
			return(-Inf)
	} else {
		return(log.lik.MN)
	}

} # Gaussian density
