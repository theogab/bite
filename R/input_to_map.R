#' @import stats
## Convert simmap or nodelabels in map object

input_to_map <- function(phy, simmap = NULL, ndlabels = NULL, map = NULL, nreg = NULL){
  
  n <- length(phy$tip.label)
  ne <- length(phy$edge.length)
  
  if(is.null(phy$Nnode)){
    S <- cbind(rep(0,ne), rep(0,ne))
    beta <- matrix(0,ne)
  } else {
    nodetime <- c(rep(max(branching.times(phy)), n), max(branching.times(phy)) - branching.times(phy))
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    pp <- prop.part(phy)
    ppt <- lapply(1:(length(phy$tip.label)+phy$Nnode), function(t) if(t <= n) t else pp[[t-n]])
    
    ## Map not provided
    if(is.null(simmap) & is.null(map) & is.null(ndlabels)){
      S <- matrix(0, ne, 2)
      beta <- matrix(1, ne)
      colnames(beta) <- ""
      gamma <- matrix(0, nrow = ne, ncol = n)
      colnames(gamma) <- phy$tip.label
      for(i in 1:ne){
        S[i,] <- c(nodetime[e1[i]], nodetime[e2[i]])
        gamma[i,ppt[[e2[i]]]] <- 1
      }
    }
    
    ## Stochastic mapping from phytools provided
    if(!is.null(simmap)){
      reg <- unique(do.call(c, lapply(simmap, names)))
      nepochs <- length(do.call(c, simmap))
      S <- matrix(0, nepochs, 2)
      beta <- matrix(0, nepochs,length(reg))
      colnames(beta) <- reg
      gamma <- matrix(0, nrow = nepochs, ncol = n)
      colnames(gamma) <- phy$tip.label
      j <- 1
      for(i in 1:ne){
        x <- simmap[[i]]
        S[j,] <- c(nodetime[e1[i]], nodetime[e1[i]] + x[1])
        beta[j,which(reg %in% names(x)[1])] <- 1 
        gamma[j,ppt[[e2[i]]]] <- 1
        j <- j + 1
        if(length(x) > 1){
          for(k in 2:length(x)){
            S[j,] <- c(S[j-1,2], S[j-1,2] + x[k])
            beta[j,which(reg %in% names(x)[k])] <- 1 
            gamma[j,ppt[[e2[i]]]] <- 1
            j <- j + 1
          }
        }
      }
    }
    
    ## Mapping provided with nodelabels
    if(!is.null(ndlabels)){
      S <- matrix(0, ne, 2)
      reg <- unique(ndlabels)
      beta <- matrix(0, ne,length(reg))
      colnames(beta) <- reg
      gamma <- matrix(0, nrow = ne, ncol = n)
      colnames(gamma) <- phy$tip.label
      for(i in 1:ne){
        S[i,] <- c(nodetime[e1[i]], nodetime[e1[i]] + phy$edge.length[i])
        beta[i,which(reg %in% ndlabels[e1[i]-n])] <- 1
        gamma[i,ppt[[e2[i]]]] <- 1
      }
    }
    
    ## Mapping provided with matrix
    if(!is.null(map)){
      nepochs <- sum(map>0)
      S <- matrix(0, nepochs, 2)
      reg <- colnames(map)
      beta <- matrix(0, nepochs,length(reg))
      colnames(beta) <- reg
      gamma <- matrix(0, nrow = nepochs, ncol = n)
      colnames(gamma) <- phy$tip.label
      j <- 1
      for(i in 1:ne){
        x <- map[i,]
        S[j,] <- c(nodetime[e1[i]], nodetime[e1[i]] + x[x>0][1])
        beta[j,which(reg %in% names(x)[x>0][1])] <- 1
        gamma[j,ppt[[e2[i]]]] <- 1
        x <- x[x>0]
        j <- j + 1
        if(length(x) > 1){
          for(k in 2:length(x)){
            if(x[k]>0){
              S[j,] <- c(max(S[j-1,]), max(S[j-1,]) + x[k])
              beta[j,reg %in% which(names(x)[k])]
              gamma[j,ppt[[e2[i]]]] <- 1
              j <- j + 1
            }
          }
        }
      }
    }
  }
  
  return(list(S = S, beta = beta, gamma=gamma))
  
}

map_to_simmap<- function(phy, map){
  
  S <- map$S
  beta <- map$beta
  n <- length(phy$tip.label)
  st <- colnames(map$beta)
  nodetime <- c(rep(max(branching.times(phy)),n), max(branching.times(phy)) - branching.times(phy))
  e1 <- phy$edge[, 1]
  e2 <- phy$edge[, 2]
  
  ## maps
  phy$maps <- list()
  j <- 1
  for(i in 1:length(e1)){
    phy$maps[[i]] <- numeric(0)
    while(round(S[j,2],5) != round(nodetime[e2[i]],5)){
      phy$maps[[i]] <- c(phy$maps[[i]], S[j,2]-S[j,1])
      names(phy$maps[[i]])[length(phy$maps[[i]])] <- colnames(beta)[beta[j,]==1]
      j <- j + 1
    }
    phy$maps[[i]] <- c(phy$maps[[i]], S[j,2]-S[j,1])
    names(phy$maps[[i]])[length(phy$maps[[i]])] <- colnames(beta)[beta[j,]==1]
    j <- j + 1
  }
  
  ## mapped.edge
  phy$mapped.edge <- t(sapply(phy$maps, function(x){
    sapply(st, function(a){
      sum(x[a], na.rm = T)
    })
  }))
  
  if(length(st) == 1){
    phy$mapped.edge <- t(phy$mapped.edge)
  }
  
  rownames(phy$mapped.edge) <- sprintf("%s,%s", phy$edge[,1], phy$edge[,2])
  class(phy) <- c("phylo", "simmap")
  return(phy)
}