
slidingWin <- function(i, d) {
	# Slidign window proporal unconstrained at maximum 
	# For details of the method see http://people.sc.fsu.edu/~pbeerli/BSC-5936/10-12-05/Lecture_13.pdf
	#
	# Args:
	# 	i:  current value
	#	d:  window size
	#
	# Returns:
	#	Proposal value (integer).
	

    ii <- i + (runif(length(i), 0, 1) - 0.5) * d #MrBayes trick
    return(ii)
}     
     

multiplierProposal <- function(i, d, u) {
	# Multiplier proporal 
	# For details of the method see http://people.sc.fsu.edu/~pbeerli/BSC-5936/10-12-05/Lecture_13.pdf
	#
	# Args:
	# 	i:  current value
	#	d:  window size
	#	u:  a random value from a uniform distribution [0,1]
	#
	# Returns:
	#	Proposal value (integer).
	
    lambda <- 2 * log(d)
    m <- exp(lambda * (u - 0.5))
    ii <- i * m
    return(ii)
}
