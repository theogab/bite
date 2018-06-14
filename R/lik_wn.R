#conditional prior on sd of species-specific normal likelihoods under WN model, x - vector of values
# require(ape)
# require(MASS)

likWN<-function(pars, x, tree, scaleHeight, ...){#M - ancestral mean, S - trend, S0 - starting point of a trend, ti - total phylogenetic time, sig.sq  - sigma^2 (phylogenetic variance)
	
	Y      <- as.matrix(x)	
	sig.sq <- pars[1] # sigma
	tvcv   <- vcv(tree)
	n      <- dim(tvcv)[1]
	vcv.m  <- matrix(nrow=n,ncol=n,0)
	if (scaleHeight) { # Rescale tree to unit length
		vcv.m <- vcv.m/max(diag(vcv.m))
	} else { # Keep original tree length
		diag(vcv.m) <- tvcv[1]		
	}	
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


# tree<-rtree(10)
# tree<-chronopl(tree, 0.5, 12)
# likWN(c(3, 235), c(2,3,4,3,3,2,3,4,4,3), tree)
