#' @import ape

# input: pars - c(sv1, sv2, ..., sig1, sig2,...,the1...theN, root), map, vcv and vector of size 4 giving the number of regime every parameter is undergoing for nr 
# does: calculate log-likelihood; see Butler and King 2004, appendix eq. A8 and A9 and Beaulieu et al. 2012
update_ou <- function(pars, Pi, map){
  
  
  alpha <- pars[Pi[2,]==1]/(2*pars[Pi[1,]==1])
  sigma <- pars[Pi[2,]==1]
  theta <- pars[Pi[3,]==1]
  if(nrow(Pi) == 4) theta <- c(pars[Pi[4,]==1], theta)
  
  S <- map$S ##Mapping of the start and end of each epoch nepochs x 2
  gamma <- map$gamma ##Indicator mapping of the epochs lived by each species nepochs x n
  beta <- map$beta ##Indicator mapping of the regimes on each epoch nepochs x nreg
  
  ## Weight matrix
  var.reg <- -alpha * t((S[,2]-S[,1]) * beta) %*% gamma
  W <- exp(var.reg) * (exp(alpha*t(S[,2]*beta))-exp(alpha*t(S[,1]*beta)))%*%gamma
  
  if(nrow(Pi) == 4){
    W <- rbind(exp(colSums(var.reg)), W)
    W <- W/colSums(W)
  }
  
  E <- t(theta %*% W)
  
  ## Variance Covariance Matrix (SIGMA)
  V <- exp(t(matrix(0, ncol(gamma), ncol(gamma)) + colSums(var.reg)) + colSums(var.reg)) * (t(gamma)%*%(colSums((sigma/(2*alpha)*t(beta))*(exp(2*alpha*t(S[,2]*beta))-exp(2*alpha*t(S[,1]*beta))))*gamma))
  inv <- solve(V)
  det <- as.numeric(determinant(V)$modulus)
  
  return(list(E = E, det = det, inv = inv))
}