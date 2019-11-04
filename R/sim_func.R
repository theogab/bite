sim_pet <- function(phy, map, model, pars, ntips, bounds){
  
  # initialization
  norm_pars <- norm_func(model)
  x.val <- numeric(phy$Nnode + ntips)
  x.val[ntips+1] <- pars[["the0"]]
  names(x.val) <- c(phy$tip.label, as.character(ntips + 1:phy$Nnode))
  
  for(i in order(phy$edge[,1])){
    
    # target <- phy$edge[i,2]  
    # anc <- phy$edge[i,1]
    
    if("WN" %in% model & phy$edge[i,1] > ntips + 1){
      map[[i]] <- c(map[[which(phy$edge[,2] == phy$edge[i,1])]], map[[i]])
    }
    
    x.val[phy$edge[i,2]] <- x.val[phy$edge[i,1]]
    for(j in 1:length(map[[i]])){
      x <- x.val[phy$edge[i,2]]
      t <- map[[i]][j]
      reg <- as.numeric(names(map[[i]])[j])
      n.pars <- norm_pars(x, pars, t)
      x.val[phy$edge[i,2]] <- reflect(rnorm(1, n.pars[reg,1], sqrt(n.pars[reg,2])), bounds) 
    }
  }
  
  return(x.val)
}

norm_func <- function(model){
  if(any(c("BM", "WN") %in% model)){
    out <- function(x, pars, t){
      sig2 <- pars[["sig"]]
      cbind(x, sig2*t)
    }
  }
  
  if("OU" %in% model){
    out <- function(x, pars, t){
      the <- pars[["the"]]
      sig2 <- pars[["sig"]]
      alp <- pars[["alp"]]
      cbind(the + (the - x) * exp(-alp * t), (sig2/(2*alp)) * (1 - exp(-2 * alp * t)))
    }
  }
  return(out)
}

## function to respect bounds
reflect <- function(yy, bounds) {
  while (yy < bounds[1] || yy > bounds[2]) {
    if (yy < bounds[1]) 
      yy <- 2 * bounds[1] - yy
    if (yy > bounds[2]) 
      yy <- 2 * bounds[2] - yy
  }
  return(yy)
}