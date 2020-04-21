# input: m.sp = vector of mean, logv.sp = vector of log(sigma^2), traits = matrix of species observations, counts = number of observation by species
# does: calculate individual log-likelihoods for each species based on normal distribution
lik_norm <- function(pars.lik, traits, counts){#m - mean (horizontal), s - sigma^2 (horizontal), vec - observations for a species
	
  if(!"mean" %in% names(pars.lik)) stop("'mean' parameter cannot be found in model.priors")
  if(!"logvar" %in% names(pars.lik)) stop("'logvar' parameter cannot be found in model.priors")
  m.sp <- pars.lik$mean
  logv.sp <- pars.lik$logvar
  
	log.lik.N <- -counts/2 * log(2 * pi) - 1/2 * counts * logv.sp - 1/2 * (sapply(1:length(traits), function(i) sum((traits[[i]] - m.sp[i])^2))/exp(logv.sp))
	
	if (is.na(sum(log.lik.N))) {
			return(-Inf)
	} else {
		return(log.lik.N)
	}

} # Gaussian density
