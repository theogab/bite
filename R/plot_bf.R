#' @title Plots summary of Bayes Factors calculations 
#' @description Lolipop plot representing values of BF (Bayes Factor) scores for different models 
#' @param m.liks matrix of marginal likelihoods (see details)
#' @param thr value of BF threshold for model selection (default is 2), Several thresholds can be given in the form of a vector
#' @param dir string giving the direction of the plot ("vertical", "horizontal")
#' @param col color of the lines and dots. Could be of size one or more. The first element applies to the best model, the second to the ones below thr[1] and so on...
#' @param col.thr color for the threshold area (must be at the sime length as thr)
#' @param ax.lab label of the axis.
#' @param main an overall title for the plot
#' @param rank logical: should models be ranked by BF scores? (default is TRUE)
#' @param dec logical: should models be displayed in decreasing order? (default is TRUE)
#' @param group.pattern Regular expression given a pattern to be matched in names(m.liks). The values matching that pattern will be grouped in the plot
#' @param cex size of the dots and width of the lines. Could be of size one or more. The first element applies to the best model, the second to the ones below thr[1] and so on...
#' @param space only evaluated if group.pattern != NULL, numeric(2) indicating the space between groups and the space between members of the same groups
#' @param mod.lab name of the models. If mod.lab = NULL the names of m.liks are used
#' @param srt.lab,adj.lab rotation and justification of model names, see srt and adj in \code{\link[graphics]{par}}
#' 
#' @author Theo Gaboriau
#' @export
#' @examples
#' ## Fake marginal likelihood data
#' m.liks <- c(20 ,33, 56, 51, 55, 12)
#' names(m.liks) <- c("MBM-VBM", "MBM-VWN", "MBM-VOU", "MOU-VBM", "MOU-VWN", "MOU-VOU")
#' 
#' #Does not go well with default margin sizes
#' par(mar = c(5,1,4,6))
#' plot_bf(m.liks)
#' plot_bf(m.liks, thr = c(2,6), col.thr = c("#a2c5ac", "#ade1e5"),
#'  col = c("#d32f23", "#468189","#2e86ab","#000000"), cex = c(1.2,1,0.8,0.8))
#' plot_bf(m.liks, group.pattern = "MBM", rank = FALSE)
#' plot_bf(m.liks, group.pattern = "MOU", rank = TRUE)
#' plot_bf(m.liks, group.pattern = c("VWN", "VOU", "VBM"), rank = TRUE)
#' par(mar = c(6,5,1,2))
#' plot_bf(m.liks, dir = "horizontal", srt.lab = -60, adj.lab = c(0,0.8))
#' 
#' @encoding UTF-8

plot_bf <- function(m.liks, thr = 2, dir = c("vertical", "horizontal"), col = c("#d32f23","#2e86ab","#000000"),
                    col.thr = c("#a6e1fa"), ax.lab = "log(BF)", main = "", rank = T, dec = T,
                    group.pattern = NULL, cex = c(1.2,1,1), mod.lab = NULL, srt.lab = 0, adj.lab = 0, space = c(1.2,0.8)){
  
  BF <- 2*(max(m.liks, na.rm = T) - m.liks)
  
  if(dir[1] == "vertical"){
    plot(0, ylim = c(0,length(BF)), xlim = c(0,max(BF)*(11/10)), xaxt = "n", yaxt = "n", bty = "n", xlab = ax.lab, ylab = "", col = "#000000", main = main, type = "n")
    for(i in 1:length(thr)){
      if(i == 1) low = 0 else low = thr[i-1]
      polygon(c(low,low,thr[i],thr[i]), c(-length(BF)/10, length(BF)+1, length(BF)+1, -length(BF)/10), border = NA, col = col.thr[i])
    }
    
    a <- axis(1, labels = F, tick = F)
    for(k in a){
      lines(c(k,k),c(-length(BF)/10,length(BF)+1), col = ifelse(k == 0, "#000000", "#cecece"), lty = ifelse(k == 0, 2,1), lwd = par()$lwd*ifelse(k == 0, 2,1))
    }
    axis(1, at = a, labels = a)
    
    if(!is.null(group.pattern)){
      find.groups <- lapply(group.pattern, function(pattern) grepl(pattern, names(BF)))
      if(sum(colSums(do.call(rbind, find.groups)) == 0) > 0){
        find.groups[[length(find.groups)+1]] <- colSums(do.call(rbind, find.groups)) == 0
      }
      groups <- lapply(find.groups, which)
      y <- cumsum(unlist(lapply(groups, function(n) c(space[1], rep(space[2],length(n)-1)))))
      y <- (y/max(y+space[1])) * length(BF)
    } else {
      y <- 1:length(BF)
      groups <- list(1:length(BF))
    }
    
    if(rank){
      it <- unlist(lapply(groups, function(x) x[order(BF[x], decreasing = T)]))
      if(!dec){
        y <- rev(y)
      }
    } else {
      it <- unlist(groups)
    }
    
    for(i in 1:length(it)){
      
      ## Define plotting parameters for every lolipop
      k <- length(cex) - sum(BF[it[i]] <= c(0,thr))
      cex.k <- cex[ifelse(k <= 0, 1, k)]
      k <- length(col) - sum(BF[it[i]] <= c(0,thr))
      col.k <- col[ifelse(k <= 0, 1, k)]
     
      ## Plot lolipop
      points(BF[it[i]], y[i], pch = 16, cex = cex.k, col = col.k)
      lines(c(BF[it[i]], par()$usr[2] + max(BF)*1/100), c(y[i],y[i]), lwd = par()$lwd*cex.k, col = col.k)
      
      ## Plot legend
      if(is.null(mod.lab)) lab.k <- names(BF)[it[i]]
      else lab.k <- mod.lab[it[i]]
      text(x = par()$usr[2] + max(BF)*1/100, y = y[i], labels = lab.k, las = 1, cex = cex.k, col = col.k, xpd = T, pos = 4, srt = srt.lab)
    }
  } else if(dir[1] == "horizontal") {
    
    plot(0, xlim = c(.5,length(BF)), ylim = c(-max(BF)*(11/10),0), xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = ax.lab, col = "#FFFFFF", main = main, type = "n")
    for(i in 1:length(thr)){
      if(i == 1) low = 0 else low = thr[i-1]
      polygon(c(-length(BF)/10, length(BF)+1, length(BF)+1, -length(BF)/10), c(low,low,-thr[i],-thr[i]), border = NA, col = col.thr[i])
    }

    a <- axis(2, labels = F, tick = F)
    for(k in a){
      lines(c(-length(BF)/10,length(BF)+1),c(k,k), col = ifelse(k == 0, "#000000", "#cecece"), lty = ifelse(k == 0, 2,1), lwd = par()$lwd*ifelse(k == 0, 2,1))
    }
    axis(2, at = a, labels = -a)
    
    if(!is.null(group.pattern)){
      find.groups <- lapply(group.pattern, function(pattern) grepl(pattern, names(BF)))
      if(sum(colSums(do.call(rbind, find.groups)) == 0) > 0){
        find.groups[[length(find.groups)+1]] <- colSums(do.call(rbind, find.groups)) == 0
      }
      groups <- lapply(find.groups, which)
      x <- cumsum(unlist(lapply(groups, function(n) c(space[1], rep(space[2],length(n)-1)))))
      x <- (x/max(x+space[1])) * length(BF)
    } else {
      x <- 1:length(BF)
      groups <- list(1:length(BF))
    }
    
    if(rank){
      it <- unlist(lapply(groups, function(x) x[order(BF[x], decreasing = T)]))
      if(dec){
        x <- rev(x)
      }
    } else {
      it <- unlist(groups)
    }
    
    for(i in 1:length(it)){
      
      ## Define plotting parameters for every lolipop
      k <- length(cex) - sum(BF[it[i]] <= c(0,thr))
      cex.k <- cex[ifelse(k <= 0, 1, k)]
      k <- length(col) - sum(BF[it[i]] <= c(0,thr))
      col.k <- col[ifelse(k <= 0, 1, k)]
      
      ## Plot lolipop
      points(x[i], -BF[it[i]], pch = 16, cex = cex.k, col = col.k)
      lines(c(x[i],x[i]), c(-BF[it[i]], par()$usr[3]), lwd = par()$lwd*cex.k, col = col.k)
      
      ## Plot legend
      if(is.null(mod.lab)) lab.k <- names(BF)[it[i]]
      else lab.k <- mod.lab[it[i]]
      text(x = x[i], y = par()$usr[3] - max(BF)/100, labels = lab.k, cex = cex.k, col = col.k, xpd = T, srt = srt.lab, adj = adj.lab)
    }
  } else {stop(sprintf("The direction '%s' is not supported", dir[1]))}
  
    
}

