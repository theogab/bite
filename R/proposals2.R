

make.proposal <- function(prop, i=1, d=1, ...){

	if (prop == "slidingWin"){
		prop.f <- function(i, d, u=0) {
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
			return(list(v=ii, lnHastingsRatio=0))
		} 
	}

	if (prop == "slidingWinAbs"){
		prop.f <- function(i, d, u=0) {
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
			return(list(v=abs(ii), lnHastingsRatio=0))
		}
	}	
	
	if (prop == "logSlidingWinAbs"){
		prop.f <- function(i, d, u=0) {
			# Slidign window proporal unconstrained at maximum 
			# For details of the method see http://people.sc.fsu.edu/~pbeerli/BSC-5936/10-12-05/Lecture_13.pdf
			#
			# Args:
			# 	i:  current value
			#	d:  window size
			#
			# Returns:
			#	Proposal value (integer).
			
			i  <- log(i)
			ii <- i + (runif(length(i), 0, 1) - 0.5) * d #MrBayes trick

			return(list(v=exp(ii), lnHastingsRatio=0))
		}
	}	
	
		 
	if (prop == "multiplierProposal"){
		prop.f <- function(i, d, u) {
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
			return(list(v=ii, lnHastingsRatio=log(u)))
		}
	}
	
	
	if (prop == "multiplierProposalLakner"){
		prop.f <- function(i, d, u) {
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
			# TODO: if d = 1 then error -> needs fixing
			
			
			lambda=2*log(d)
			m=exp(lambda*(u-0.5))
			ii=i*m
			
			return(list(v=ii, lnHastingsRatio=log(m)))
			
		}
	}
	
	
	
	if (prop == "logNormal"){
		# i - current value, d - sigma (sd) of normal dist
		prop.f <- function(i, d, u=0){
			 
			ii <- rnorm(1, mean = i, sd = d)

			return(list(v=ii, lnHastingsRatio=0))
		}
		
	}
	
	if (prop == "absNormal"){
		# i - current value, d - sigma (sd) of normal dist
		prop.f <- function(i, d, u=0){
			 
			ii <- rnorm(1, mean = i, sd = d)

			return(list(v=abs(ii), lnHastingsRatio=0))
		
		}
		
	}

	

	return(prop.f)

}

# hastings ratio
# i - is a mean of normal distribution in current state log sd 
# ii - is a mean of normal distribution in proposed state l
calc.hasting.ratio <-  function(i, ii, d){
	
	prob.i  <- dlnorm(i, meanlog = ii, sdlog = d, log=TRUE)
	prob.ii <- dlnorm(ii, meanlog = i, sdlog = d, log=TRUE)
	
	# minus because of log
	hs <- prob.ii - prob.i
	
	return(hs)
}

