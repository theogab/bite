

heat_par <- function(ncat=10, beta.param=0.3){ 
	# Defines classes for thermodynamic integration.
	# For details of the method see Xie et al 2011 Sys Bio.
	#
	# Args:
	# 	ncategories: number of classes that will be used in thermodynamic integration.
	#	beta.param:  parameter describing the shape of a beta distribution.
	# 
	# Returns:
	#	The vector of temperatures for thermodynamic integration.
	
    K <- ncat-1
    k <- 0:K
    b <- k/K
    temp<- rev(b^(1/beta.param))
    temp[length(temp)] <- 0.00001 # last category is not exactly 0 to avoid -inf likelihoods
	return(temp)
}
