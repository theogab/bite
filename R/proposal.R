

proposal <- function(prop, i=1, d=1, ...){

  # checking
  if(!prop %in% c("slidingWin", "slidingWinAbs", "logSlidingWinAbs",
                  "multiplierProposal", "multiplierProposalLakner",
                  "logNormal", "absNormal") ){
    stop("Unknown proposal algorithm")
  }
  
	if (prop == "slidingWin"){
	  slidingWin <- function(i, d, u=0) {
			# Sliding window proposal unconstrained at maximum 
			# For details of the method see http://people.sc.fsu.edu/~pbeerli/BSC-5936/10-12-05/Lecture_13.pdf
			#
			# Args:
			# 	i:  current value
			#	  d:  window size
			#
			# Returns:
			#	Proposal value (integer).
			
			ii <- i + (runif(length(i), 0, 1) - 0.5) * d #MrBayes trick
			return(list(v=ii, lnHastingsRatio=0))
	  } 
	  
	  return(slidingWin)
	}

	if (prop == "slidingWinAbs"){
	  slidingWinAbs <- function(i, d, u=0) {
			# Sliding window proposal unconstrained at maximum 
			# For details of the method see http://people.sc.fsu.edu/~pbeerli/BSC-5936/10-12-05/Lecture_13.pdf
			#
			# Args:
			# 	i:  current value
			#	  d:  window size
			#
			# Returns:
			#	Proposal value (integer).
			
			ii <- i + (runif(length(i), 0, 1) - 0.5) * d #MrBayes trick

			return(list(v=abs(ii), lnHastingsRatio=0))
	  }
	  
	  return(slidingWinAbs)
	}	
	
	if (prop == "logSlidingWinAbs"){
	  logSlidingWinAbs <- function(i, d, u=0) {
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
	  
	  return(logSlidingWinAbs)
	}	
	
		 
	if (prop == "multiplierProposal"){
	  multiplierProposal <- function(i, d, u) {
			# Multiplier proposal 
			# For details of the method see http://people.sc.fsu.edu/~pbeerli/BSC-5936/10-12-05/Lecture_13.pdf
			#
			# Args:
			# 	i:  current value
			#	  d:  window size
			#	  u:  a random value from a uniform distribution [0,1]
			#
			# Returns:
			#	Proposal value (integer).
			
			
			lambda <- 2 * log(d)
			m <- exp(lambda * (u - 0.5))
			ii <- i * m

			return(list(v=ii, lnHastingsRatio=log(u)))
	  }
	  
	  return(multiplierProposal)
	}
	
	
	if (prop == "multiplierProposalLakner"){
	  multiplierProposalLakner <- function(i, d, u) {
			# Multiplier proposal 
			# For details of the method see http://people.sc.fsu.edu/~pbeerli/BSC-5936/10-12-05/Lecture_13.pdf
			#
			# Args:
			# 	i:  current value
			#	  d:  window size
			#	  u:  a random value from a uniform distribution [0,1]
			#
			# Returns:
			#	Proposal value (integer).
			# TODO: if d = 1 then error -> needs fixing
			
			
			lambda=2*log(d)
			m=exp(lambda*(u-0.5))
			ii=i*m

			return(list(v=ii, lnHastingsRatio=log(m)))
			
	  }
	  return(multiplierProposalLakner)
	  
	}
	
	
	
	if (prop == "logNormal"){
		# i - current value, d - sigma (sd) of normal dist
	  logNormal <- function(i, d, u=0){
			 
			ii <- rnorm(1, mean = i, sd = d)

			return(list(v=ii, lnHastingsRatio=0))
		}
		
	  return(logNormal)
	  
	}
	
	if (prop == "absNormal"){
		# i - current value, d - sigma (sd) of normal dist
	  absNormal <- function(i, d, u=0){
			 
			ii <- rnorm(1, mean = i, sd = d)

			return(list(v=abs(ii), lnHastingsRatio=0))
		
		}
		return(absNormal)
	}

}


