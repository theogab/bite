# input: n - number of species, mat - list of e (expectation matrix), det (determinant of the varcovar) and inv (inverse of the varcovar), x - var or mean of trait by sp
# does: calculate log-likelihood; 
calc_prior <- function(n, mat, x){
 
  E <- mat$E
  det <- mat$det
  inv <- mat$inv
  
  #loglik <- try((-n/2 * log(2 * pi) - 1/2 * det - 1/2 * (t(x - E)%*%inv%*%(x - E))), silent=T)
  loglik <- (-n/2 * log(2 * pi) - 1/2 * det - 1/2 * (t(x - E)%*%inv%*%(x - E)))
  
  if (is.na(loglik) | (class(loglik) == "try-error" )) {
    return(-Inf)
  } else {
    return(loglik)
  }
}