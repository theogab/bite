#' @title Plots summary of Bayes Factors calculations 
#' @description Lolipop plot representing values of BF (Bayes Factor) scores for different models 
#' @param m.liks matrix of marginal likelihoods (see details)
#' @param thr value of BF threshold for model selction (default is 2)
#' @param col.thr color for the threshold area
#' @param ylab name of the models. If ylab = "default" the dimnames of m.liks are used
#' @param rank should models be ranked by BF scores? (default is TRUE)
#' @param group only evaluated if rank == F, character indicating whether the models should be grouped regarding the evolutionnary model on "mean" or "var" 
#' @param space only evaluated if rank == F, numeric(2) indicating the space between groups and the space between members of the same groups
#' @param cex size of the dots. Could be of size one two or three. The first element applies to models that are above thr, the last element applies to the best model, the middle element applies to models that are below thr
#' @param lwd.l thickness of the lines. Could be of size one two or three. The first element applies to models that are above thr, the last element applies to the best model, the middle element applies to models that are below thr
#' @param col color of the lines and dots. Could be of size one two or three. The first element applies to models that are above thr, the last element applies to the best model, the middle element applies to models that are below thr
#' @param col.axis color of axis annotation. Could be of size one two or three. The first element applies to models that are above thr, the last element applies to the best model, the middle element applies to models that are below thr
#' @author Theo Gaboriau
#' @export
#' @examples
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' 
#' ## Run a MCMC chain with thermodynamic Integration
#' mean_models <- c("BM", "OU")
#' var_models <- c("BM", "OU", "BMM", "OUM")
#' mliks <- matrix(0, nrow = 2, ncol = 4, dimnames = list(mean_models, var_models))
#' 
#' for(mn in mean_models){
#'  for(vr in var_models){
#'   my.jive <- make_jive(Anolis_tree, Anolis_traits,  model.var=vr, model.mean=mn)
#'   logfile <- sprintf("my.jive_MCMC_M-%s_V-%s.log", mn, vr)
#'   mcmc_jive(my.jive, log.file=logfile, ncat=10, sampling.freq=10, print.freq=100, ngen=5000, burnin=500)
#'   res <- read.table(logfile, header = T, sep = "\t")
#'   mliks[mn,vr] <- marginal.lik(res)
#'  }
#' }
#' 
#' plot_bf(mliks)
#' plot_bf(mliks, thr = 4)
#' plot_bf(mliks, thr = 4, rank = F, group = "mean")
#' plot_bf(mliks, dir = "horizontal", srt.lab = -60, adj.lab = c(.3,-1))

plot_bf <- function(m.liks, thr = 2, dir = c("vertical", "horizontal"), col = c("#000000","#2e86ab","#d32f23"),
                    col.thr = c("#a6e1fa"), ax.lab = "Bayes Factor", main = "", mod.lab = "default", rank = T,
                    group = c("mean", "var"), cex = c(1,1,2), lwd.l = c(1,1,2), srt.lab = 0, adj.lab = 0,
                    col.axis =  c("#000000","#2e86ab","#d32f23"), cex.axis = c(.8, .8, 1), space = c(1.2,0.8)){
  
  BF <- max(m.liks, na.rm = T) - m.liks
  
  if(dir[1] == "vertical"){
    plot(0, ylim = c(0,length(BF)), xlim = c(0,max(BF)*(11/10)), xaxt = "n", yaxt = "n", bty = "n", xlab = ax.lab, ylab = "", col = "#000000", main = main)
    polygon(c(0,0,thr,thr), c(-length(BF)/10, length(BF)+1, length(BF)+1, -length(BF)/10), border = NA, col = col.thr)
    
    a <- axis(1, labels = F, tick = F)
    for(k in a){
      lines(c(k,k),c(-length(BF)/10,length(BF)+1), col = ifelse(k == 0, "#000000", "#cecece"), lty = ifelse(k == 0, 2,1), lwd = par()$lwd*ifelse(k == 0, 2,1))
    }
    axis(1, at = a, labels = a, cex.axis = cex.axis[1], col = col.axis[1])
    
    if(rank){
      it <- sapply(BF[order(BF, decreasing = T)], function(x) which(BF == x, arr.ind = T))
      y <- 1:ncol(it)
    } else if(group == "mean"){
      it <- rbind(rep(1:nrow(BF), each = ncol(BF)), rep(1:ncol(BF),nrow(BF)))
      y <- cumsum(rep(c(space[1], rep(space[2],ncol(BF)-1)),nrow(BF)))
    } else if(group == "var"){
      it <- rbind(rep(1:nrow(BF), ncol(BF)), rep(1:ncol(BF),each = nrow(BF)))
      y <- cumsum(rep(c(space[1], rep(space[2], nrow(BF)-1)), ncol(BF)))
    } else {
      stop("sort must be 'mean' or 'var'")
    }
    
    for(k in 1:ncol(it)){
      
      i <- it[1,k]
      j <- it[2,k]
      
      if(length(cex) == 1) 
        cex.k <- cex
      if(length(cex) == 2) 
        cex.k <- ifelse(BF[i,j] == 0, cex[2], cex[1])
      if(length(cex) == 3)
        cex.k <- ifelse(BF[i,j] == 0, cex[3], ifelse(BF[i,j] <= thr, cex[2], cex[1]))
      if(length(lwd.l) == 1)
        lwd.k <- lwd.l
      if(length(lwd.l) == 2)
        lwd.k <- ifelse(BF[i,j] == 0, lwd.l[2], lwd.l[1])
      if(length(lwd.l) == 3)
        lwd.k <- ifelse(BF[i,j] == 0, lwd.l[3], ifelse(BF[i,j] <= thr, lwd.l[2], lwd.l[1]))
      if(length(col) == 1)
        col.k <- col
      if(length(col) == 2)
        col.k <- ifelse(BF[i,j] == 0, col[2], col[1])
      if(length(col) == 3)
        col.k <- ifelse(BF[i,j] == 0, col[3], ifelse(BF[i,j] <= thr, col[2], col[1]))
      if(length(col.axis) == 1)
        col.axis.k <- col.axis
      if(length(col.axis) == 2)
        col.axis.k <- ifelse(BF[i,j] == 0, col.axis[2], col.axis[1])
      if(length(col.axis) == 3)
        col.axis.k <- ifelse(BF[i,j] == 0, col.axis[3], ifelse(BF[i,j] <= thr, col.axis[2], col.axis[1]))
      if(length(cex.axis) == 1)
        cex.axis.k <- cex.axis
      if(length(cex.axis) == 2)
        cex.axis.k <- ifelse(BF[i,j] == 0, cex.axis[2], cex.axis[1])
      if(length(cex.axis) == 3)
        cex.axis.k <- ifelse(BF[i,j] == 0, cex.axis[3], ifelse(BF[i,j] <= thr, cex.axis[2], cex.axis[1]))
      
      points(BF[i,j], y[k], pch = 16, cex = cex.k, col = col.k)
      lines(c(BF[i,j], max(BF)*11/10), c(y[k],y[k]), lwd = lwd.k, col = col.k)
      
      if(mod.lab == "default")
        lab.k <- paste0(c("M", "V"), c(rownames(BF)[i], colnames(BF)[j]), collapse = " - ")
      else 
        lab.k <- ylab[k]
      text(x = max(BF)*11/10, y = y[k], labels = lab.k, las = 1, cex = cex.axis.k, col = col.axis.k, xpd = T, pos = 4, srt = srt.lab)
    }
  } else if(dir[1] == "horizontal") {
    
    plot(0, xlim = c(.5,length(BF)), ylim = c(-max(BF)*(11/10),0), xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = ax.lab, col = "#FFFFFF", main = main)
    polygon(c(-length(BF)/10, length(BF)+1, length(BF)+1, -length(BF)/10), c(0,0,-thr,-thr), border = NA, col = col.thr)
    
    a <- axis(2, labels = F, tick = F)
    for(k in a){
      lines(c(-length(BF)/10,length(BF)+1),c(k,k), col = ifelse(k == 0, "#000000", "#cecece"), lty = ifelse(k == 0, 2,1), lwd = par()$lwd*ifelse(k == 0, 2,1))
    }
    axis(2, at = a, labels = -a, cex.axis = cex.axis[1], col = col.axis[1])
    
    if(rank){
      it <- sapply(BF[order(BF, decreasing = T)], function(x) which(BF == x, arr.ind = T))
      x <- 1:ncol(it)
    } else if(group == "mean"){
      it <- rbind(rep(1:nrow(BF), each = ncol(BF)), rep(1:ncol(BF),nrow(BF)))
      x <- cumsum(rep(c(space[1], rep(space[2],ncol(BF)-1)),nrow(BF)))
    } else if(group == "var"){
      it <- rbind(rep(1:nrow(BF), ncol(BF)), rep(1:ncol(BF),each = nrow(BF)))
      x <- cumsum(rep(c(space[1], rep(space[2], nrow(BF)-1)), ncol(BF)))
    } else {
      stop("sort must be 'mean' or 'var'")
    }
    
    for(k in 1:ncol(it)){
      
      i <- it[1,k]
      j <- it[2,k]
      
      if(length(cex) == 1) 
        cex.k <- cex
      if(length(cex) == 2) 
        cex.k <- ifelse(BF[i,j] == 0, cex[2], cex[1])
      if(length(cex) == 3)
        cex.k <- ifelse(BF[i,j] == 0, cex[3], ifelse(BF[i,j] <= thr, cex[2], cex[1]))
      if(length(lwd.l) == 1)
        lwd.k <- lwd.l
      if(length(lwd.l) == 2)
        lwd.k <- ifelse(BF[i,j] == 0, lwd.l[2], lwd.l[1])
      if(length(lwd.l) == 3)
        lwd.k <- ifelse(BF[i,j] == 0, lwd.l[3], ifelse(BF[i,j] <= thr, lwd.l[2], lwd.l[1]))
      if(length(col) == 1)
        col.k <- col
      if(length(col) == 2)
        col.k <- ifelse(BF[i,j] == 0, col[2], col[1])
      if(length(col) == 3)
        col.k <- ifelse(BF[i,j] == 0, col[3], ifelse(BF[i,j] <= thr, col[2], col[1]))
      if(length(col.axis) == 1)
        col.axis.k <- col.axis
      if(length(col.axis) == 2)
        col.axis.k <- ifelse(BF[i,j] == 0, col.axis[2], col.axis[1])
      if(length(col.axis) == 3)
        col.axis.k <- ifelse(BF[i,j] == 0, col.axis[3], ifelse(BF[i,j] <= thr, col.axis[2], col.axis[1]))
      if(length(cex.axis) == 1)
        cex.axis.k <- cex.axis
      if(length(cex.axis) == 2)
        cex.axis.k <- ifelse(BF[i,j] == 0, cex.axis[2], cex.axis[1])
      if(length(cex.axis) == 3)
        cex.axis.k <- ifelse(BF[i,j] == 0, cex.axis[3], ifelse(BF[i,j] <= thr, cex.axis[2], cex.axis[1]))
      
      points(x[k], -BF[i,j], pch = 16, cex = cex.k, col = col.k)
      lines(c(x[k],x[k]), c(-BF[i,j],par()$usr[3]), lwd = lwd.k, col = col.k)
      
      if(mod.lab == "default")
        lab.k <- paste0(c("M", "V"), c(rownames(BF)[i], colnames(BF)[j]), collapse = " - ")
      else 
        lab.k <- ylab[k]
      text(y = par()$usr[3] - max(BF)/100, x = x[k], labels = lab.k, cex = cex.axis.k,
           col = col.axis.k, xpd = T, srt = srt.lab, adj = adj.lab)
    }
  } else {stop(sprintf("The direction '%s' is not supported", dir[1]))}
  
    
}

