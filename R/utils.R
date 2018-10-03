

defClasses <- function(ncat=10, beta.param=0.3){ 
	# Defines classes for thermodynamic integration.
	# For details of the method see Xia et al 2011 Sys Bio.
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

       
initUpdateFreq <- function(update.freq=NULL){
	# Initializes update frequencies for likelihood and two prior levels.
	#
	# Args:
	# 	update.freq: the vector (length = 3) of update frequencies (likelihood, priorMBM, priorVOU/VBM).
	#
	# Returns:
	#	The vector of update frequencies which sums to 1. 
	
	if (length(update.freq) != 3 && !is.null(update.freq)) {
		stop("Update.freq must contain 3 elements" )
	}
	
	
	# calculate update frequencies
	if (!is.null(update.freq)) {
		update.freq	<- cumsum(update.freq/sum(update.freq))
	} else {
		update.freq	<- cumsum(c(0.35,0.2,0.45))
	}
	
	return(update.freq)

}


## tranform simmap into map, input - simmap object
relSim <- function(x) {
	
	foo<-function(x) {
		x/sum(x)
	}
	
	x$mapped.edge <- t(apply(x$mapped.edge, 1, FUN=foo))
	x$mapped.edge <- x$mapped.edge[, order(colnames(x$mapped.edge))]
	
	return(x)
				
}






