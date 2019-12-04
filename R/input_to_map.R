#' @import stats
## Convert simmap or nodelabels in map object

input_to_map <- function(phy, simmap = NULL, ndlabels = NULL, map = NULL, nreg = NULL){
  
  n <- length(phy$tip.label)
  
  if(is.null(phy$Nnode)){
    newmap <- lapply(phy$tip.label, function(sp) rbind(1:nreg, rep(0,nreg), rep(0,nreg)))
  } else {
    nodetime <- max(branching.times(phy)) - branching.times(phy)
    newmap <- list()
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    
    ## Map not provided
    if(is.null(simmap) & is.null(map) & is.null(ndlabels)){
      for(i in 1:nrow(phy$edge)){
        newmap[[e2[i]]] <- matrix(c(1,nodetime[e1[i]-n], ifelse(e2[i] > n, nodetime[e2[i]-n], max(branching.times(phy)))), ncol = 1)
      }
    }
    
    ## Stochastic mapping from phytools provided
    if(!is.null(simmap)){
      reg <- unique(do.call(c, lapply(simmap, names)))
      for(i in 1:nrow(phy$edge)){
        x <- simmap[[i]]
        newmap[[e2[i]]] <- matrix(c(which(names(x)[1] == reg), nodetime[e1[i]-n], nodetime[e1[i]-n] + x[1]), ncol = 1)
        if(length(x) > 1){
          for(j in 2:length(x)){
            newmap[[e2[i]]] <- cbind(newmap[[e2[i]]], c(which(names(x)[j] == reg), newmap[[e2[i]]][3,j-1], newmap[[e2[i]]][3,j-1] + x[j]))
          }
        }
      }
    }
    
    ## Mapping provided with nodelabels
    if(!is.null(ndlabels)){
      reg <- unique(ndlabels)
      for(i in 1:nrow(phy$edge)){
        newmap[[e2[i]]] <- matrix(c(which(ndlabels[phy$edge[i,1]-n] == reg), nodetime[phy$edge[i,1]-n], nodetime[phy$edge[i,1]-n] + phy$edge.length[i]), ncol = 1)
      }
    }
    
    ## Mapping provided with matrix
    if(!is.null(map)){
      reg <- colnames(map)
      for(i in 1:nrow(phy$edge)){
        x <- map[i,]
        newmap[[e2[i]]] <- matrix(c(which(names(x)[x>0][1] == reg), nodetime[phy$edge[i,1]-n], nodetime[phy$edge[i,1]-n] + x[x>0][1]), ncol = 1)
        x <- x[x>0]
        if(length(x) > 1){
          for(j in 2:length(x)){
            newmap[[e2[i]]] <- cbind(newmap[[e2[i]]], c(which(names(x)[j] == reg), newmap[[e2[i]]][3,j-1], newmap[[e2[i]]][3,j-1] + x[j]))
          }
        }
      }
    }
  }
  
  
  return(newmap)
}

map_to_simmap<- function(phy, map){
  n <- length(phy$tip.label)
  st <- as.character(unique(do.call(cbind,map)[1,]))
  
  ## maps
  phy$maps <- lapply(map, function(x){
    out <- x[3,] - x[2,]
    names(out) <- x[1,]
    return(out)
  })[phy$edge[,2]]
  
  ## mapped.edge
  phy$mapped.edge <- t(sapply(phy$maps, function(x){
    sapply(st, function(a){
      sum(x[a], na.rm = T)
    })
  }))
  rownames(phy$mapped.edge) <- sprintf("%s,%s", phy$edge[,1], phy$edge[,2])
  class(phy) <- c("phylo", "simmap")
  return(phy)
}

