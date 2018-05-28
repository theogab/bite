# input: [make.simmap output, single run]
# does: transform simmap$mapped.edges from lengths into relative probabilities (e.g. sums to 1 for each row)
relSim <- function(simmap){

	foo<-function(x){
		x/sum(x)
	}
	simmap$mapped.edge <- t(apply(simmap$mapped.edge, 1, FUN=foo))
	return(simmap)

}


# input: alpha, T - total tree length, s - branch length
# does: transform branch length according to OU-model; see Butler and King 2004, appendix eq. A5
vcvMat <- function(alpha, T, s){

	v.mat <- exp(-2 * alpha * (T-s)) * (1-exp(-2 * alpha * s))
	return(v.mat)
}


# input: i - lineage number (e.g. starting node), j - regime, e - output from exp.mat, n.sp - number of species
# does: traverse lineage and sums already exponentiated branch length; does it repeatedly for each regime column
travLineage <- function(i, j, e, n.sp, t.sum=0){
	
	add = e[which(e[, 2] == i), j] 
	if ( i == n.sp + 1){ 
		return(t.sum)
	} else {
		t.sum <- t.sum + add 
		pa <- e[e[, 2] == i][1]
		t.sum <- travLineage(pa, j, e, n.sp, t.sum)
	} 
}

# input: tree - phylo.object , map - output from relSim function [matrix], alpha, T.len - total tree length
# does: take nodes, assign age to each node, exponentiate node age and calculate the difference  between each from-to pair; see Butler and King 2004, appendix eq. A4 and A7 
expMat <- function(tree, map, alpha, T.len){
	
	
	d <- cbind(tree$edge,tree$edge.length)
	br <- branching.times(tree)
	#names(br)[1]<-as.character(as.numeric(names(br)[2])-1) #<- THIS MUST BE CODED IN A BETTER WAY TO ENSURE NAMES OF THE BRANCHING TIMES
	d <- cbind(d, T.len - br[match(d[, 1], names(br))], T.len - br[match(d[, 2], names(br))])
	d[which(is.na(d[, 5])), 5] <- T.len
	e <- exp(-alpha * T.len) * (exp(alpha * d[, 5]) - exp(alpha * d[, 4]))
	e <- cbind(tree$edge, exp(-alpha * T.len), e * map)
	return(e)
	
}


# input: n.reg - number of regimes (without theta0), e - output from expMat, n.sp - number of species
# does: construct a weight matrix;  see Butler and King 2004, appendix eq. A7 
weigthMat <- function(n.reg, e, n.sp ){
	
	w.mat <- matrix(nrow=n.sp, ncol=(n.reg + 1), 0)
	w.mat[, 1]<- e[1, 3] #get the zero regime
	
	for (j in 1:n.reg){
		for (i in 1:n.sp){
			w.mat[i, j+1]=travLineage(i, j+3, e, n.sp) 
		}
	}
	
	return(w.mat)
}



# input: pars - c(alpha,sig.sq,theta0,theta1...thetaN), sigma.val, tree and map (output from relSim)
# does: calculate log-likelihood; see Butler and King 2004, appendix eq. A8 and A9 
likOU <- function(pars, x, tree, map){ #x - vertical, theta - horizontal, pars - alpha, sig.sq, thetas
	
	 # extract variables
	 alpha  <- pars[1]
	 sig.sq <- pars[2]
	 n.reg  <- length(pars)-3
	 theta  <- pars[3:length(pars)]
	 t.vcv  <- vcv(tree)
	 T.len  <- t.vcv[1, 1]
	 n      <- dim(t.vcv)[1]
	 
	 # calculate matricies
	 t.vcv <- sig.sq/(2 * alpha) * (vcvMat(alpha, T.len, t.vcv))
	 e     <- expMat(tree, map, alpha, T.len)
	 w     <- weigthMat(n.reg, e, n)
	 DET   <- determinant(t.vcv, logarithm=T)

	 log.lik.OU <- try((-n * log(2 * pi)/2 - (as.numeric(DET$modulus))/2 - (t(x - w%*%theta)%*%ginv(t.vcv)%*%(x - w%*%theta))/2),silent=T)

	 #print(log.lik)
	 if (is.na(log.lik.OU) | (class( log.lik.OU) == "try-error" )) {
		return(-Inf)
	 } else {
		return(log.lik.OU)
	 }

}



