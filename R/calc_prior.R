# input: n - number of species, mat - list of e (expectation matrix), det (determinant of the varcovar) and inv (inverse of the varcovar), x - var or mean of trait by sp
# does: calculate log-likelihood; 
calc_prior <- function(n, mat, x){
 
  e <- mat[[1]]
  det <- mat[[2]]
  inv <- mat[[3]]
  
  loglik <- try((-n/2 * log(2 * pi) - 1/2 * det - 1/2 * (t(x - e)%*%inv%*%(x - e))), silent=T)
  
  if (is.na(loglik) | (class(loglik) == "try-error" )) {
    return(-Inf)
  } else {
    return(loglik)
  }
}