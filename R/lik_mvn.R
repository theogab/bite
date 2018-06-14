# input: pars=c(anc.mean,sig.sq), mean.val - trait values, tree - tree
# does: return a log-likelihood of a classical BM-model (multivariate normal distribution)

likMVN <- function(pars, x, tree, scaleHeight=FALSE){#M - ancestral mean, sig.sq - sigma^2, x - species means
	
	Y      <- as.matrix(x)
	M      <- pars[1] # ancestral mean
	sig.sq <- pars[2] # sigma
	vcv.m  <- vcv(tree)
	if (scaleHeight) { # Rescale tree to unit length
		vcv.m <- vcv.m/max(diag(vcv.m))
	}
	n      <- dim(vcv.m)[1]
	m      <- matrix(1, n, 1)
	m[, ]  <- M
	DET    <- determinant(sig.sq * vcv.m, logarithm=T)

	log.lik.MVN <- try((-n/2 * log(2 * pi) - (as.numeric(DET$modulus))/2 - 1/2 * (t(Y - m)%*%ginv(sig.sq * vcv.m)%*%(Y - m))), silent=T)
	
	if (is.na(log.lik.MVN) | (class(log.lik.MVN) == "try-error" )) {
			return(-Inf)
	} else {
		return(log.lik.MVN)
	}

} 


