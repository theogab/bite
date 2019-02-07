sim_pet <- function(phy, map, model, pars, ntips){
  
  # initialization
  norm_pars <- norm_func(model)
  x.val <- numeric(phy$Nnode + ntips)
  x.val[ntips+1] <- pars[[1]][1]
  names(x.val) <- c(phy$tip.label, as.character(ntips + 1:phy$Nnode))
  
  for(i in order(phy$edge[,1])){
    
    # target <- phy$edge[i,2]  
    # anc <- phy$edge[i,1]
    
    if(model %in% c("WN", "WNM") & phy$edge[i,1] > ntips + 1){
      map[[i]] <- c(map[[which(phy$edge[,2] == phy$edge[i,1])]], map[[i]])
    }
    
    x.val[phy$edge[i,2]] <- x.val[phy$edge[i,1]]
    for(j in 1:length(map[[i]])){
      x <- x.val[phy$edge[i,2]]
      t <- map[[i]][j]
      reg <- as.numeric(names(map[[i]])[j])
      n.pars <- norm_pars(x, pars, t)
      x.val[phy$edge[i,2]] <- rnorm(1, n.pars[reg,1], sqrt(n.pars[reg,2]))
    }
  }
  
  return(x.val)
}

norm_func <- function(model){
  if(model %in% c("BM", "BMM", "WN", "WNM")){
    out <- function(x, pars, t){
      sig2 <- pars[[2]]
      cbind(x, sig2*t)
    }
  }
  
  if(model %in% c("OU", "OUM")){
    out <- function(x, pars, t){
      the <- pars[[1]][-1]
      sig2 <- pars[[2]]
      alp <- pars[[3]]
      cbind(the + (the - x) * exp(-alp * t), (sig2/(2*alp)) * (1 - exp(-2 * alp * t)))
    }
  }
  return(out)
}