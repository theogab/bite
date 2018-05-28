#conditional prior on sd of species-specific normal likelihoods under BM model, x - vector of values
likBMM<-function(pars, x, tree, ...){#M - ancestral mean, S - trend, S0 - starting point of a trend, ti - total phylogenetic time, sig.sq  - sigma^2 (phylogenetic variance)
	
	Y      <- as.matrix(x)	

	sig.sq <- pars[1] # sigma
	vcv.m  <- vcv(tree)
	n      <- dim(vcv.m)[1]
	m      <- matrix(1, n, 1)
	m[, ]  <- pars[2] # ancestral mean
	DET    <- determinant(sig.sq * vcv.m, logarithm=T)
	
	log.lik.BM <- try((-n/2 * log(2 * pi) - (as.numeric(DET$modulus))/2 - 1/2 * (t(Y - m)%*%ginv(sig.sq * vcv.m)%*%(Y - m))), silent=T)
	
	if (is.na(log.lik.BM) | (class(log.lik.BM) == "try-error" )) {
			return(-Inf)
	} else {
		return(log.lik.BM)
	}

}


