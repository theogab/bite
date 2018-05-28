#' @title Simulate trees
#' @description Simulate species trees with complex diversification scenarios
#' 
#' @details The aim of the function is to simulate species trees with complex diversification scenarios,
#' with changing diversification rates at inner nodes of the phylogeny and mass extinction events at 
#' particular random times in the past. The user can choose the total number of species, the number 
#' of diversification rate changes, the number of mass extinction events, the minimum and maximum 
#' values for lambda and mu, and the minimum survival rate for mass extinction events. Then the actual rates are choosen randomly.
#' User gets as output the tree simulated, and all the parameters used for its simulation.
#' 
#' 
#' 
#' 
#' @param nb.divers.shifts number of diversification shifts on the tree
#' @param nb.mass.ext number of mass extinction events
#' @param nb.spec total number of species in the final tree
#' @param lambda.max maximum value of speciation rate		
#' @param lambda.min minimum value of speciation rate
#' @param mu.min minimum value of extinction rate
#' @param surv.rate.min survival rate
#' @param seed seed value for reproducibility
#' @export
#' @author Sacha Laurent
#' @return An object of class phylo  and a set of parameters used for simulation.
#' @examples
#' tree <- classicalModel(1, 1, 500, 0.1, 0.01, .01, .2, seed=1506621)


classicalModel <- function(nb.divers.shifts, nb.mass.ext, nb.spec, lambda.max, lambda.min, mu.min, surv.rate.min, seed=1){
		# Simulate species trees with complex diversification scenarios
		# User should provide simulation parameters
		#
		# Args:
		# 	nb.divers.shifts number of diversification shifts on the tree
		#	nb.mass.ext number of mass extinction events
		#	nb.spec total number of species in the final tree
		#   lambda.max maximum value of speciation rate		
		# 	lambda.min minimum value of speciation rate
		# 	mu.min minimum value of extinction rate
		# 	surv.rate.min survival rate
		# 	seed seed value for reproducibility
		# 			
		# Returns:
		#	An object of class phylo  and a set of parameters used for simulation.

		set.seed(seed)
		print(seed)
		vect_nb_spec <- rep(0, nb.divers.shifts+1)
		while (sum(vect_nb_spec<2) | sum(vect_nb_spec)!=nb.spec)   vect_nb_spec <- round(MCMCpack::rdirichlet(1, rep(1, nb.divers.shifts+1))*nb.spec)
		param <- data.frame(l=rep(0, nb.divers.shifts+1), m=rep(0, nb.divers.shifts+1), spec=t(vect_nb_spec))
		for (i in 1:(nb.divers.shifts+1)){
		param[i, "l"] <- runif(1, min=lambda.min, max=lambda.max)
		param[i, "m"] <- runif(1, min=mu.min, max=param[i, "l"])
		}
		#  param <- data.frame(l=runif(nb.divers.shifts+1, min=lambda.min, max=lambda.max))
		# print(param)
		# param <- cbind(param, data.frame(m=runif(nb.divers.shifts+1, min=mu.min, max=param[1, "l"]), spec=t(vect_nb_spec)))
		print(param)
				 
		rownames(param) <- paste("s", 1:(nb.divers.shifts+1), sep="")
		mee_times <- c(0, runif(nb.mass.ext, min=0, max=min(log(param[,"spec"])/(param[,"l"]-param[,"m"]))))  
		surv_rate <- c(1, runif(nb.mass.ext, min=surv.rate.min, max=1))
		mee_times <- sort(mee_times)
		spec_names <- c(0, cumsum(vect_nb_spec)) ##definition of parameters
		br.ti <- c(0)
		print(param)
		print(mee_times)
		print(surv_rate)
		while (max(br.ti) < tail(mee_times, 1) | !max(br.ti)){
		tree <- list()
		for (i in 1:(nb.divers.shifts+1)){ #simulation of all sub trees
		  tmp <- sim.rateshift.taxa(param[i, "spec"], 1, rep(param[i, "l",], length(surv_rate)), rep(param[i, "m"], length(surv_rate)), surv_rate, mee_times, complete = TRUE, K=0, norm = TRUE)
		  tree[[i]] <- list()
		  if (Ntip(tmp[[1]][[1]]) > vect_nb_spec[i]) tree[[i]][[1]] <- drop.tip(tmp[[1]][[1]], (vect_nb_spec[i]+1):Ntip(tmp[[1]][[1]])) #deletes extinct tips
		  else tree[[i]][[1]] <- tmp[[1]][[1]]
		  tree[[i]][[1]]$tip.label <- paste("t",(spec_names[i]+1):spec_names[i+1], sep="")
		  tree[[i]][[1]]$node.label <-  paste("n", (spec_names[i]+1):(spec_names[i+1]-1), sep="")
		  tree[[i]][[1]]$node.label[1] <- paste("s", i, sep="") #we mark the root as special to keep track of the changes in lambda mu
		  tree[[i]][[2]] <- tmp[[2]] #we keep the root length to bind the trees in the next step
		}
		br.ti <- sapply(tree, function(x) branching.times(x[[1]])[1])
		}
		root_times <- sapply(tree, function(x) x[[2]])
		param <- cbind(param, root_times, mrca_times=br.ti)
		print(param)
		return(c(sortBindAllTrees(tree, param), list(mee=mee_times, seed=seed, survival=surv_rate)))
}



