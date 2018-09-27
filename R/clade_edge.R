## Select all branches that share a comon ancestor

clade.edge <- function(phy,node){
  x <- phy$edge[phy$edge[, 1] == node, 2]
  y <- which(phy$edge[,1] == node)
  repeat {
    xx <- x
    x <- sort(unique(c(x, phy$edge[, 2][phy$edge[, 1] %in% x])))
    y <- c(y,which(is.element(phy$edge[,1],x)))
    if (identical(x, xx)) break
  }
  y <- as.numeric(levels(as.factor(y)))
  return(y)
}