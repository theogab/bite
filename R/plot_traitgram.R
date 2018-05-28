#' @title Plot traitgrams object
#' @description Plots traitgrams under different evolutionary models
#' 
#' @details This function plots the dynamics of trait evolution under Ornstein Uhlenbeck (OU) models.
#' The user needs to specify the parameters of an OU model (it can include multiple regimes OUM, OUMV..),
#' which are typically obtained from trait optimization (e.g. using OUwie package). A user-defined 
#' number of trait simulations are generated under the provided model and parameters, and assuming a 
#' trait value at the root estimated under Brownian motion. Each simulation records trait values at tips and nodes. 
#' Traitgrams (cf. Phytools or Picante functions) are inferred by estimating the trait trajectories along the 
#' branches based on the OU bridge (see Serrano-Serrano et. al in review). The function generates a plot of the 95% 
#' CI of the trait space through time based on all simulations and/or traitgrams of single simulations with branches
#' colored in function of the OU regimes.	
#' 
#' @param fitted.model a fitted object of class "OUwie" (see details)
#' @param anc.model a name of the model for ancstral state reconstruction (see details)
#' @param reps number of replicates for simulating a trait 
#' @param plot.grey.only if to plot only grey area or more					
#' @param ... additional parameters passed to XXX functions
#' @export
#' @author Martha Serrano and Daniele Silvestro
#' @return nothing
#' @examples
#' ## Load test data
#' library(ape)
#' library(phytools)
#' library(OUwie)
#' tree1		<- rtree(40)
#' tree1 		<- chronopl(tree1, 1)
#' regime1		<- rep(c(0,1),each=20)
#' names(regime1) <- tree1$tip.label
#' sim1 		<- make.simmap(tree1, regime1, model="ER", nsim=1, pi="estimated")
#' traits1 <- OUwie.sim(sim1, simmap.tree=TRUE, alpha=c(1.5,1.5), 
#' sigma.sq=c(0.1,0.1), theta0=2.5, theta=c(1.0,6.0))
#' traits1 <- data.frame(traits1[,1],regime1,traits1[,2])
#' OUwie_OUM <- OUwie(sim1, traits1, model="OUM", simmap.tree=TRUE)
#' par(mfrow=c(1,2))
#' plotTraitgram(fitted.model=OUwie_OUM, anc.model="ER", reps=5, plot.grey.only = "TRUE")
#' plotTraitgram(fitted.model=OUwie_OUM, anc.model="ER", reps=5, plot.grey.only = "FALSE")





plotTraitgram <- function (fitted.model, anc.model, reps=5, plot.grey.only = "FALSE", ...){
	
	# parse fitted.model
	tree			<- fitted.model$phy
	regime			<- (fitted.model$data)[,1] - 1
	trait			<- (fitted.model$data)[,2]
	names(regime)	<- rownames(fitted.model$data)
	names(trait)	<- rownames(fitted.model$data)

	# initial vectors
	minY		<- 0
	maxY		<- 0
	varY		<- 0
	
	q			<- 1-branching.times(tree)

	table.min	<- matrix(nrow=length(q), ncol=reps)
	table.max	<- matrix(nrow=length(q), ncol=reps)
	table.var	<- matrix(nrow=length(q), ncol=reps)
	

	alpha		<- c(fitted.model$solution[1, 1], fitted.model$solution[1, 2])
	sigma.sq	<- c(fitted.model$solution[2, 1], fitted.model$solution[2, 2])
	theta		<- c(fitted.model$theta[1, 1],	fitted.model$theta[2, 1])
	theta0		<- fitContinuous(tree, trait, model ="BM")$opt$z0 # ancestral root is difficult to estimate when multiple OU regimes, BM used
	
	#return(list(alpha,sigma.sq,theta,theta0))
	
	sim.traits	<- list()
	
	# simulate stoch map reps number of times
    sim.obj		<- make.simmap(tree, regime, model=anc.model, nsim=reps, pi="estimated")
	
	for (i in c(1:reps)){
		
		print(i)
		
		# simulating a trait with all parameters and stochastic mapped tree of regimes 
		y	<- OUwie.sim2(sim.obj[[i]], simmap.tree=TRUE, scaleHeight=F, alpha=alpha, sigma.sq=sigma.sq, theta0=theta0, theta=theta)
		x	<- as.numeric(y[[1]])
		y2	<- y[[2]][,1]
		names(x) <- c(y2, c((length(y2) + 1):length(x)))
		root <- length(sim.obj[[i]]$tip) + 1
		
		# based on the edge number (=74 branches, coordenates where it starts and where it finish)
		node.height <- matrix(NA ,nrow(sim.obj[[i]]$edge), 3)  
		
		for(j in 1:nrow(sim.obj[[i]]$edge)){
		
			if(sim.obj[[i]]$edge[j,1] == root){
		
				node.height[j, 1] <- 0.0
				node.height[j, 2] <- sim.obj[[i]]$edge.length[j]
		
			} else {
		
				node.height[j, 1] <- node.height[match(sim.obj[[i]]$edge[j, 1], sim.obj[[i]]$edge[, 2]), 2]
				node.height[j, 2] <- node.height[j, 1] + sim.obj[[i]]$edge.length[j]
		
			} 
		}

		sim.traits[[i]] <- x   # list of simulated traits for the single realization traitgram
      	  
		# trait range is estimated at each branching time 
		q <- sort(q)
		traitSimul <- matrix(ncol = 4, nrow = length(q))
 
	    # loop over each branching time
		for(s in 3:length(q)){   # # first branching time is root, starts at 2-second btime
			
			tsplit	<- as.numeric(q[s])
			edges 	<- which(node.height[, 2]>tsplit&node.height[, 1]<tsplit)
			b 		<- sim.obj[[i]]$edge[edges, ]  # nodes that belongs to the branches cut at point tsplit
			wtrait  <- vector()  # collecting traits at this point time, but weigthed by bl from/to each node
	
			# loop over all branches cut at s branching time, result estimates a values of traits, range and var at that point 
			if(length(b) == 2) { 
				len <- length(b) 
			} else {
				len <- nrow(b)
			}
		
			for (l in 1:len) {	
				
				if (length(b) == 2) {
					p <- x[b[l]]
				  
				} else {			
					p <- x[b[l, ]] 
				}
				
				if (names(p)[2] == names(q[s])){	
				
					v1 <- as.numeric(p[2])
					
				} else { 
				
					reg <- as.numeric(names(sim.obj[[i]]$maps[[edges[l]]]))	#annotated regime on edges of simmap
					v0  <- as.numeric(p[1])   #starting trait
					v2  <- as.numeric(p[2])
					re  <- sample(reg,1)+1  # transform character states into index and randomly draw a regime in case of multiple per branch 
					t0  <- as.numeric(q[which(names(q)==as.numeric(names(p[1])))])
					if ((names(p[2]))%in% sim.obj[[i]]$tip.label) {
						t2 <- 1	
					} else { 
						t2 <- as.numeric(q[which(names(q)==as.numeric(names(p[2])))])
					}
					v1 <- theta[re] + (v0 - theta[re]) * sinh(alpha[re]*(t2 - tsplit)) / sinh(alpha[re] * (t2 - t0)) + (v2 - theta[re]) * sinh(alpha[re] * (tsplit - t0)) / sinh(alpha[re] * (t2 - t0))
					
				} 
				wtrait <- append(wtrait, v1) 
			}
	 
			names(wtrait)	 <- edges
			traitSimul[s, 1] <- tsplit
			traitSimul[s, 2] <- min(wtrait)
			traitSimul[s, 3] <- max(wtrait)
			traitSimul[s, 4] <- var(wtrait)
		}
	
	n	 <- length(sim.obj[[i]]$tip.label)
	minY <- append(minY, min(x[1:n]))
	maxY <- append(maxY, max(x[1:n]))
	varY <- append(varY, var(x[1:n]))
	table.min[, i] <- traitSimul[, 2] 
	table.max[, i] <- traitSimul[, 3] 
	table.var[, i] <- traitSimul[, 4]
	}
  
	########## END OF FIRST PART
    #return(list(n, minY, maxY, varY, table.min, table.max, table.var))

    summary_min 	<- apply(table.min, 1, mean)
    summary_max 	<- apply(table.max, 1, mean)
    summary_min[1]	<- theta0
    summary_max[1] 	<- theta0
    summary_min[length(tree$tip.label)] <- mean(minY)
    summary_max[length(tree$tip.label)] <- mean(maxY)
    slices 			<- round(q, 3)
    slices 			<- c(slices, get.rooted.tree.height(tree))
    summary_min		<- cbind(slices, summary_min)
    summary_max 	<- cbind(slices, summary_max)
   
    probs <- c(0.05, 0.95)
  
    trait_q5 		<- apply(table.min, 1, function(x) quantile(x, probs, na.rm =TRUE)[1])
    trait_q5[1] 	<- theta0
	
    trait_q5[length(tree$tip.label)]  <- quantile(minY, probs)[1]
	
    trait_q95		<- apply(table.max, 1, function(x) quantile(x, probs, na.rm=TRUE)[2])
    trait_q95[1] 	<- theta0
	
    trait_q95[length(tree$tip.label)] <- quantile(maxY, probs)[2]
	
    var_q5			<- apply(table.var, 1, function(x) quantile(x, probs, na.rm=TRUE)[1])
    var_q5[1] 		<- 0  # first branching time is root...zero variance
	
    var_q5[length(tree$tip.label)]	  <- quantile(varY, probs)[1]
	
    var_q95			<- apply(table.var, 1, function(x) quantile(x, probs, na.rm=TRUE)[2])
    var_q95[1]		<- 0
	
    var_q95[length(tree$tip.label)]	  <- quantile(varY, probs)[2]

	
	## plot param pre-calculation
	
	## x-limits
	xll = 0
	xuu = vcv(tree)[1]
	mmin = min(trait)
	mmax = max(trait)
	dd = mmax - mmin
	yll = mmin - dd
	yuu = mmax + dd
	step.size = xuu/10000
	#print(list(xll,xuu,yll,yuu))
	

	if (plot.grey.only == TRUE){
	
		plot(0, 0.95, type="n", ylim=c(yll, yuu), xlim=c(xll, xuu),  ylab="Estimated Trait range", xlab="Time")
		polygon(c(slices, rev(slices)), c(trait_q95, rev(trait_q5)), col = "#C0C0C0AA", border = NA)
		lines(summary_min[, 1], summary_min[, 2], col="black", lwd=3)
		lines(summary_max[, 1], summary_max[, 2], col="black", lwd=3) 
		segments(x0=xuu, x1=xuu, y0=min(trait), y1=max(trait))
	
	} else {
		
		x		<- sim.traits[[1]]    # specify one of the simulations number
		tree	<- sim.obj[[1]]
		
		q		<- seq(0.00, xuu, step.size)
		tgram	<- matrix(nrow=length(tree$edge.length), ncol=length(q)+1)
		
		for (edge in 1:nrow(node.height)){
		
			trait_edge	<- x[tree$edge[edge, ]]
			trange		<- node.height[edge, 1:2]
			reg			<- as.numeric(names(tree$maps[[edge]]))	#annotated regime on edges of simmap
			re			<- sample(reg, 1)+1  # transform character states into index and randomly draw a regime in case of multiple per branch       
			
			tgram[edge,length(q) + 1]	<- re
			
			for(s in 1:length(q)){   # # first branching time is root, starts at 2-second btime
				xpoint <- q[s]
				if (xpoint > trange[1] & xpoint <trange[2]){	
					v0 <- as.numeric(trait_edge[1])   #starting trait
					v2 <- as.numeric(trait_edge[2])
					t0 <- trange[1] 
					
					if ((names(trait_edge[2]))%in% tree$tip.label) {
						t2 <- 1	
					} else { 
						t2 <- trange[2]
					}
					t1 <- xpoint
					v1 <- theta[re] + (v0 - theta[re]) * sinh(alpha[re] * (t2-t1)) / sinh(alpha[re] * (t2-t0)) + (v2-theta[re]) * sinh(alpha[re] * (t1-t0)) / sinh(alpha[re] * (t2-t0))
					tgram[edge, s] <- v1
				 } else { 
					if (xpoint==round(trange[1], 2)) { 
						tgram[edge, s] <- as.numeric(trait_edge[1])  # if xpoint (slice q) is equal to starting 
					} else { 
						if (xpoint==round(trange[2], 2)){
						tgram[edge, s] <- as.numeric(trait_edge[2])
						} else { 
							tgram[edge, s] <- tgram[edge, s]
						}
					}
				 }
			}# close loop over time
		}# close loop over edges
		class(tgram) <- "numeric"
		
		plot(0, 0.95, type="n", ylim=c(yll, yuu),xlim=c(xll, xuu),  ylab="Estimated Trait range", xlab="Time")
		polygon(c(slices, rev(slices)), c(trait_q95, rev(trait_q5)), col = "#C0C0C0AA", border = NA)

		for (r in 1:dim(tgram)[1]){
		
			if (length(which( !is.na(tgram[r,]), arr.ind=TRUE))==1){
				points(y=tgram[r,], x=q, pch=20, col="black")
			} else {
				col=ifelse(tgram[r,length(q)+1]==1,"#d8b365","#5ab4ac")
				y=tgram[r,1:length(q)]
				lines(y=y,x=q, col=col,lwd=2)
			}
		}

		segments(x0=xuu, x1=xuu, y0=min(trait), y1=max(trait), lwd=2)
	
	}
	
  
}