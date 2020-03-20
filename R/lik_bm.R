#' @import ape

# input: n - number of species, n.p - which parameter has been updated, pars - c(sig1, ..., sigN, the0), tree and map
# does: calculate log-likelihood; 
update_bm <- function(pars, Pi, map){
  
  sigma <- pars[Pi[1,]==1]
  root <- pars[Pi[2,]==1]
  
  S <- map$S ##Mapping of the start and end of each epoch nepochs x 2
  gamma <- map$gamma ##Indicator mapping of the epochs lived by each species nepochs x n
  beta <- map$beta ##Indicator mapping of the regimes on each epoch nepochs x nreg

	## Expectation
	E  <- matrix(1, ncol(gamma), 1)
	E[,] <- root # ancestral mean
	rownames(E) <- colnames(gamma)
	
	var.reg <- sigma * t((S[,2]-S[,1]) * beta)
	# Variance Covariance Matrix (SIGMA)
	V <-  t(gamma)%*%(colSums(var.reg)*gamma)
	# determinant
  det <- as.numeric(determinant(V)$modulus)
  # inverse
  inv <- solve(V)
	
	return(list(E = E, det = det, inv = inv))

}
