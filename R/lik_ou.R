#' @import ape

# input: pars - c(alp1, alp2, ..., sig1, sig2,...,the1...theN, root), map, vcv and vector of size 4 giving the number of regime every parameter is undergoing for nr 
# does: calculate log-likelihood; see Butler and King 2004, appendix eq. A8 and A9 and Beaulieu et al. 2012
lik_ou <- function(x, n, pars, Pi, par.n, data, map){
  
  
  alpha <- pars[Pi[1,]==1]
  sigma <- pars[Pi[2,]==1]
  theta <- pars[Pi[3,]==1]
  if(sum(Pi[1,]) > sum(Pi[3,])|sum(Pi[1,]) > sum(Pi[3,])) theta <- rep(theta, max(sum(Pi[1,]), sum(Pi[2,])))
  if(nrow(Pi) == 4) theta <- c(pars[Pi[4,]==1], theta)
  upd <- apply(Pi, 2,function(x) which(x == 1))[par.n]
  
  S <- map$S ##Mapping of the start and end of each epoch nepochs x 2
  gamma <- map$gamma ##Indicator mapping of the epochs lived by each species nepochs x n
  beta <- map$beta ##Indicator mapping of the regimes on each epoch nepochs x nreg
  
  ## Weight matrix
  var.reg <- -alpha * crossprod((S[,2]-S[,1]) * beta, gamma)
  W <- exp(var.reg) * (exp(alpha*t(S[,2]*beta))-exp(alpha*t(S[,1]*beta)))%*%gamma
  
  if(nrow(Pi) == 4){
    W <- rbind(exp(colSums(var.reg)), W)
    W <- W/colSums(W)
  }
  
  ## Expectation
  if(any(upd %in% c(1,3,4))){
    data$E <- t(theta %*% W)
  }
  
  ## Variance Covariance Matrix (SIGMA)
  if(any(upd %in% c(1,2))){
    V <- exp(t(matrix(0, ncol(gamma), ncol(gamma)) + colSums(var.reg)) + colSums(var.reg)) * crossprod(gamma,colSums((sigma/(2*alpha)*t(beta))*(exp(2*alpha*t(S[,2]*beta))-exp(2*alpha*t(S[,1]*beta))))*gamma)
    data$inv <- solve(V)
    data$det <- as.numeric(determinant(V)$modulus)
  }
  
  ## Log likelihood (multinormal)
  loglik <- -n/2 * log(2 * pi) - 1/2 * data$det - 1/2 * crossprod(x - data$E, data$inv)%*%(x - data$E)
  
  if (is.na(loglik)) {
    return(list(loglik = -Inf, data = data))
  } else {
    return(list(loglik = loglik, data = data))
  }
}
