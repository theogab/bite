#' @title Plot posterior probabilities of paramenter value under multiple regimes estimated with the beast implementation of JIVE, OU, BM and WN.
#' @description This function plots the phylogenetic tree along with mean posterior probabilities of the chosen parameter. 
#' Options are included to plot node support and age bars
#'
#' @param mcc A character containing the path to the Maximum credibility tree from the Beast analysis and extracted from treeAnnotator
#' @param post A logical specifying whether the mean posterior value of a parameter should be plot along the tree
#' @param post.var A character containing the name of the parameter (only evaluated if post == TRUE)
#' @param post.cex A numeric giving the size of the dots (only evaluated if post == TRUE)
#' @param post.col Color of the dots (only evaluated if post == TRUE)
#' @param post.alpha A numeric giving the transparency of the dots (only evaluated if post == TRUE)
#' @param post.border Color of the dots' border (only evaluated if post == TRUE)
#' @param post.lwd Line width of the dots' border (only evaluated if post == TRUE)
#' @param leg A logical specifying whether a legend should be plotted (only evaluated if post == TRUE)
#' @param leg.frac A vector giving the proportion of the plot window taken by the legend (only evaluated if post == TRUE)
#' @param leg.lab Label of the legend (only evaluated if post == TRUE)
#' @param sup A logical specifying whether node support information should be displayed
#' @param sup.var A character containing the name of the support variable (only evaluated if sup == TRUE)
#' @param sup.cex A numeric giving the size of node support information (only evaluated if sup == TRUE)
#' @param bars A logical specifying whether node age bars should be plotted
#' @param bar.var A character containing the name of the node age variable. MIN and MAX should be available for this variable (only evaluated if bars == TRUE)
#' @param bar.col Color of the bars (only evaluated if bars == TRUE)
#' @param bar.alpha A numeric giving the transparency of the bars (only evaluated if bars == TRUE)
#' @param bar.thc A numeric giving the thickness of the bars (only evaluated if bars == TRUE)
#' @param bar.border Color of the bars' border (only evaluated if bars == TRUE)
#' @param bar.lwd Line width of the bars' border (only evaluated if bars == TRUE)
#' @param ... additional parameters that can be parsed to \code{\link[ape]{plot.phylo}}
#' @export
#' @import ape
#' @author Theo Gaboriau
#' @return plot 


plot_post_beast <- function(mcc, post = T, post.var = "rate", post.cex = par("cex"), post.col = "#ee964b",
                          post.alpha = 1, post.border = "#000000", post.lwd = par("lwd"), leg = T,
                          leg.frac = c(.2,.5), leg.lab = "rate", sup = T, sup.var = "posterior", sup.cex = par("cex")/2,
                          bars = F, bar.var = "height_95%_HPD", bar.col = "#2e86ab", bar.alpha = 0.8,
                          bar.thc = par("cex")/10, bar.border = "#000000", bar.lwd = par("lwd"),  ...){
  
  tree <- read.nexus(mcc)
  n <- length(tree$tip.label)
  data <- beast_data(mcc, tree)
  
  plot(tree, type = "phylogram", plot = FALSE, ...)
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  mrg <- par("mar")
  if(pp$direction == "leftwards") par(mar = mrg[c(1,4,3,2)]) 
  if(pp$direction == "upwards") par(mar = mrg[c(2,3,4,1)]) 
  if(pp$direction == "downwards") par(mar = mrg[c(4,1,2,3)]) 
  
  if(post) label.offset = pp$label.offset + post.cex*2
  else label.offset = pp$label.offset
  
  if(post & leg){
    if(pp$direction == "rightwards") par(new = TRUE, fig = c(leg.frac[1],1,0,1), bty = "n", xpd = NA)
    if(pp$direction == "leftwards") par(new = TRUE, fig = c(0,1-leg.frac[1],0,1), bty = "n", xpd = NA)
    if(pp$direction == "upwards") par(new = TRUE, fig = c(0,1,leg.frac[1],1), bty = "n", xpd = NA)
    if(pp$direction == "downwards") par(new = TRUE, fig = c(0,1,0,1-leg.frac[1]), bty = "n", xpd = NA)
  }
  
  plot(tree, type = "phylogram", plot = TRUE, label.offset = label.offset, ...)
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  if(bars){
    bar.col <- adjustcolor(bar.col, alpha.f = bar.alpha)
    if(!is.na(bar.border)) bar.border <- adjustcolor(bar.border, alpha.f = bar.alpha)
    cc <- cbind(data[,sprintf("%s_MIN", bar.var)], data[,sprintf("%s_MAX", bar.var)])
    for(i in 1:tree$Nnode + n){
      if(pp$direction == "rightwards") rect(max(pp$xx)-cc[data$nodes == i,1], pp$yy[i]-bar.thc/2,
                                            max(pp$xx)-cc[data$nodes == i,2], pp$yy[i]+bar.thc/2,
                                            col = bar.col, lwd = bar.lwd, border = bar.border)
      if(pp$direction == "leftwards") rect(min(pp$xx)+cc[data$nodes == i,2], pp$yy[i]-bar.thc/2,
                                            min(pp$xx)+cc[data$nodes == i,1], pp$yy[i]+bar.thc/2,
                                            col = bar.col, lwd = bar.lwd, border = bar.border)
      if(pp$direction == "upwards") rect(pp$xx[i]-bar.thc/2, max(pp$yy)-cc[data$nodes == i,1],
                                          pp$xx[i]+bar.thc/2, max(pp$yy)-cc[data$nodes == i,2],
                                          col = bar.col, lwd = bar.lwd, border = bar.border)
      if(pp$direction == "downwards") rect(pp$xx[i]-bar.thc/2, min(pp$yy)+cc[data$nodes == i,2],
                                           pp$xx[i]+bar.thc/2, min(pp$yy)+cc[data$nodes == i,1],
                                           col = bar.col, lwd = bar.lwd, border = bar.border)
    }
  }
  
  if(post){
    if(length(post.col) == 1) grad <- adjustcolor(colorRampPalette(c("#FFFFFF", post.col))(101), alpha.f = post.alpha)
    else grad <- adjsutcolor(colorRampPalette(post.col)(101), alpha.f = post.alpha)
    val <- data[,post.var]
    sc.val <- round((val-min(val, na.rm = T))/(max(val, na.rm = T) - min(val, na.rm = T)) * 100, 0)
    symbols(pp$xx, pp$yy, circles = rep(post.cex/diff(par("usr")[3:4])*10, length(pp$xx)), bg = grad[sc.val + 1], fg = post.border, inches = F, add = T)
  }
  
  if(sup){
    v <- round(data[,sup.var],2)
    for(i in 1:tree$Nnode + n){
      if(pp$direction == "rightwards") text(pp$xx[i], pp$yy[i], labels = v[data$nodes == i], cex = sup.cex, adj = c(-1,-1))
      if(pp$direction == "leftwards") text(pp$xx[i], pp$yy[i], labels = v[data$nodes == i], cex = sup.cex, adj = c(2,2))
      if(pp$direction == "upwards") text(pp$xx[i], pp$yy[i], labels = v[data$nodes == i], cex = sup.cex, adj = c(2,-1))
      if(pp$direction == "downwards") text(pp$xx[i], pp$yy[i], labels = v[data$nodes == i], cex = sup.cex, adj = c(-1,2))
    }          
  }
  
  if(post & leg){
    if(pp$direction == "rightwards") par(new = TRUE, fig = c(0,leg.frac[1],leg.frac[2],1))
    if(pp$direction == "leftwards") par(new = TRUE, fig = c(1-leg.frac[1],1,leg.frac[2],1))
    if(pp$direction == "upwards") par(new = TRUE, fig = c(leg.frac[2],1,0,leg.frac[1]))
    if(pp$direction == "downwards") par(new = TRUE, fig = c(leg.frac[2],1,1-leg.frac[1],1))
    if(any(pp$direction %in% c("rightwards", "leftwards"))){
      plot(rep(0,101),seq(min(val, na.rm = T), max(val, na.rm = T), length = 101), pch = 15,
           col = grad, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
      co <- c(0,min(val, na.rm = T),0,max(val, na.rm = T)) + c(-par("cxy"),par("cxy"))*.5/2
      rect(co[1],co[2],co[3],co[4])
      axis(4, pos = co[3])
      mtext(leg.lab, side = 3, at = 0, adj = 0)
    }
    if(any(pp$direction %in% c("upwards", "downwards"))){
      plot(seq(min(val, na.rm = T), max(val, na.rm = T), length = 101), rep(0,101), pch = 15, col = grad, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
      co <- c(min(val, na.rm = T),0,max(val, na.rm = T),0) + c(-par("cxy"),par("cxy"))*.5/2
      rect(co[1],co[2],co[3],co[4])
      axis(1, pos = co[2])
      mtext(leg.lab, side = 1, at = co[3], adj = 0)
    }
  }
  
  par(fig = c(0,1,0,1), mar = mrg)
}
