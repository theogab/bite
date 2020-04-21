#' @import ape

# input: n - number of species, n.p - which parameter has been updated, pars - c(sig1, ..., sigN, the0), tree and map
# does: calculate log-likelihood; 
lik_bm <-  function(x, n, pars, Pi, par.n, data, map){
  
  sigma <- pars[Pi[1,]==1]
  root <- pars[Pi[2,]==1]
  upd <- apply(Pi, 2,function(x) which(x == 1))[par.n]
  
  S <- map$S ##Mapping of the start and end of each epoch nepochs x 2
  gamma <- map$gamma ##Indicator mapping of the epochs lived by each species nepochs x n
  beta <- map$beta ##Indicator mapping of the regimes on each epoch nepochs x nreg
  
  ## Expectation
  if(any(upd == 2)){
    data$E  <- matrix(1, ncol(gamma), 1)
    data$E[,] <- root # ancestral mean
    rownames(data$E) <- colnames(gamma)
  }
  
  # Variance Matrix (SIGMA)
  if(any(upd == 1)){
    var.reg <- sigma * t((S[,2]-S[,1]) * beta)
    # Variance Covariance Matrix (SIGMA)
    V <-  crossprod(gamma,(colSums(var.reg)*gamma))
    # determinant
    data$det <- as.numeric(determinant(V)$modulus)
    # inverse
    data$inv <- solve(V)
  }

  ## Log likelihood (multinormal)
  loglik <- (-n/2 * log(2 * pi) - 1/2 * data$det - 1/2 * (crossprod(x - data$E,data$inv)%*%(x - data$E)))

  if (is.na(loglik)) {
    return(list(loglik = -Inf, data = data))
  } else {
    return(list(loglik = loglik, data = data))
  }
  
}
