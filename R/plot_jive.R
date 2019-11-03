#' @title plot input data from a jive object
#' @description This function plots the phylogenetic tree, the trait data and the map used as an input for a jive analysis
#' @param jive a jive object built with the function \code{\link{make_jive}} 
#' @param col.map a vector of mode character indicating colors of the edges for the map ploting. It should be of same size than the number of regimes
#' @param col a character indicating the color of the vioplots
#' @param show.tip.label a logical indicating whether to show the tip labels on the phylogeny (defaults to TRUE, i.e. the labels are shown).
#' @param direction a character string specifying the direction of the tree. Two values are possible: "rightwards" (the default) and "upwards".
#' @param trait.lab a charachter specifying the axis label for the traits
#' @param trait.lim a vector of mode numeric indicating the limits for trait ploting
#' @param srt.label an integer indicating the string rotation in degrees for the tip labels
#' @param tip.color the colours used for the tip labels, eventually recycled.
#' @param ... additional parameters that can be passed to \code{\link{vioplot}}
#' @export
#' @import phytools ape vioplot sm
#' @author Theo Gaboriau
#' @return plot 
#' @examples
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' 
#' my.jive <- make_jive(Anolis_tree, Anolis_traits, 
#'  model.var="OU", model.mean="BM")
#' par(cex.lab = .8, cex.axis = .8, las = 1, mgp = c(2,0.5,0))
#' plot_jive(jive = my.jive, show.tip.label = TRUE, 
#' trait.lab = "Snout to vent length (cm)", srt.label = 0)
#' my.jive <- make_jive(Anolis_tree, Anolis_traits, Anolis_map,
#'  model.var=c("OU", "theta"), model.mean="BM")
#' par(cex.lab = .8, cex.axis = .8, las = 1, mgp = c(2,0.5,0))
#' plot_jive(jive = my.jive, show.tip.label = TRUE,
#'  trait.lab = "Snout to vent length (cm)", srt.label = 70, direction = "upwards")
#'  



plot_jive <- function(jive, col.map = NULL, col = "lightgrey", show.tip.label = T, show.models = T, direction = "rightwards",
                      trait.lab = "", trait.lim = NULL, srt.label = 0, tip.color, ...){
  
  tree <- jive$data$tree
  map <- jive$data$map
  traits <- jive$data$traits
  
  if(is.null(col.map)){
    st <- max(do.call(cbind,map)[1,])
    col.map <- palette()[1:st]
    if("simmap" %in% class(tree)) names(col.map) <- colnames(tree$mapped.edge)
    if (length(st) > 1) {
      cat("no colors provided. using the following legend:\n")
      print(col.map)
    }
  }

  if(direction == "upwards"){
    par(mar = c(1,4,1,2))
    root.len <- max(branching.times(tree))
    if(show.tip.label){
      ylim = c(0, 3*root.len)
    } else {
      ylim = c(0, 2*root.len)
    }
    
    if("simmap" %in% class(tree)){
      plotSimmap(tree, direction = "upwards", ftype = "off", ylim = ylim, mar = par()$mar, colors = col.map)
    } else {
      plot(tree, direction = "upwards", show.tip.label = F, y.lim = ylim, edge.col = col.map[unlist(sapply(tree$edge[,2], function(i){
        x <- map[[i]]
        x[1,which.max(x[3,]-x[1,])]
      }))])
    }
    init.usr <- par()$usr
    init.mar <- par()$mar
    par(new = TRUE, fig = c(0,1,ifelse(rep(show.tip.label,2),c(0.35,0.7),c(0.5,1))), bty = "n")
    plot(c(1,length(tree$tip.label)), range(unlist(traits), na.rm = T), type = "n", xaxt = "n", ylab = trait.lab,
         xlab = "", ylim = ifelse(rep(!is.null(trait.lim),2), trait.lim, range(unlist(traits),na.rm = T)))
    pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    
    for(i in 1:length(tree$tip.label)){
      sp <- tree$tip.label[i]
      vioplot(traits[[sp]][!is.na(traits[[sp]])], add = T, at = pp$xx[i], col = col, ...)
    } 
    
    par(fig=c(0,1,0,1), usr = init.usr, mar = init.mar)
    if(show.tip.label){
      text(x = pp$xx[1:length(tree$tip.label)], y = 0.7*par()$usr[4], labels = gsub("_", " ", tree$tip.label), srt = srt.label, adj = 0, xpd = NA)
    }
    
  } else if (direction == "rightwards"){
    if(show.models){
      par(mar = c(4,1,4,1))
    } else {
      par(mar = c(4,1,2,1))
    }
    
    root.len <- max(branching.times(tree))
    if(show.tip.label){
      xlim = c(0, 3*root.len)
    } else {
      xlim = c(0, 2*root.len)
    }
    
    if("simmap" %in% class(tree)){
      plotSimmap(tree, direction = "rightwards", ftype = "off", xlim = xlim, mar = par()$mar, colors = col.map)
    } else {
      plot(tree, direction = "rightwards", show.tip.label = F, x.lim = xlim, edge.col = col.map[unlist(sapply(tree$edge[,2], function(i){
        x <- map[[i]]
        x[1,which.max(x[3,]-x[1,])]
      }))])
    }
    
    if(show.models){
      mm <- strsplit(jive$prior.mean$name, " ")[[1]][1]
      mv <- strsplit(jive$prior.var$name, " ")[[1]][1]
      pars.m <- sapply(strsplit(names(jive$prior.mean$hprior), "[.]"), function(pr){
        if(length(pr)==2){
          if(pr[2] == "sig") out <- "sigma^2"
          if(pr[2] == "the") out <- "theta[0]"
        }
        if(length(pr) == 3){
          if(pr[2] == "sig") out <- "sigma^2"
          if(pr[2] == "the") out <- "theta[0]"
        }
        return(out)
      })
      mtext(eval(parse(text = paste("substitute(a~~list(", paste(pars.m, collapse = ","),"), list(a=sprintf('Mean prior model: %s' , jive$prior.mean$name)))"))), line = 2, at = 0, adj = 0, cex = 0.6)
      mtext(eval(parse(text = paste("substitute(a~~list(", paste(pars.v, collapse = ","),"), list(a=sprintf('LogVariance prior model: %s' , jive$prior.var$name)))"))), line = 2, at = 0, adj = 0, cex = 0.6)
    }
    
    init.usr <- par()$usr
    init.mar <- par()$mar
    par(new = TRUE, fig = c(ifelse(rep(show.tip.label,2),c(0.35,0.7),c(0.5,1)),0,1), bty = "n")
    plot(range(unlist(traits), na.rm = T), c(1,length(tree$tip.label)), type = "n", yaxt = "n", xlab = trait.lab,
         ylab = "", xlim = ifelse(rep(!is.null(trait.lim),2), trait.lim, range(unlist(traits),na.rm = T)))
    
    pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    
    for(i in 1:length(tree$tip.label)){
      sp <- tree$tip.label[i]
      vioplot(traits[[sp]][!is.na(traits[[sp]])], add = T, at = pp$yy[i], horizontal = T, col = col, ...)
    } 
    
    par(new = T, fig=c(0,1,0,1), usr = init.usr, mar = init.mar)
    if(show.tip.label){
      text(y = pp$yy[1:length(tree$tip.label)], x = 0.7*par()$usr[2], labels = gsub("_", " ", tree$tip.label), srt = srt.label, adj = 0, xpd = NA)
    }
    
} else stop(sprintf("%s: unknown direction", direction))
  
}

