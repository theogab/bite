# input: xmean = vector of means, xsd = vector of sigmas, traits = matrix of species observations, N_OBSERV vector of observation counts
# does: calculate individual log-likelihoods for each species based on normal distribution
lik_multinorm <- function(xmean, xsd, x, counts){#m - mean (horizontal), s - sigma^2 (horizontal), vec - observations for a species
	
	log.lik.MN <- -counts/2 * log(2 * pi) - 1/2 * counts * log(xsd) - 1/2 * (apply((x - xmean)^2, 1, sum, na.rm=T)/xsd)
	
	if (is.na(sum(log.lik.MN))) {
			return(-Inf)
	} else {
		return(log.lik.MN)
	}

} # Gaussian density

### default tuning for initial window size. The output can be manually modified using control.mcmc()
# input: x - species trait value, nreg - number of regimes 
ws_multinorm <- function(x){
  
  ws <- list()
  ws$m.sp <- 2*x
  ws$v.sp <- 10*x
  return(ws)
  
}


init_multinorm <- function (x){
  
  pv  <- list()
  pv$m.sp  <- apply(x, 1, mean, na.rm = T) # initialize means for species
  pv$v.sp  <- apply(x, 1, sd, na.rm = T) # initialize sigma.sq for species
  return(pv)
  
}
