#' @title plot input data from a jive object
#' @description This function plots the phylogenetic tree, the trait data and the map used as an input for a jive analysis
#' @param jive a jive object built with the function \code{\link{make_jive}} 
#' @param col.map a vector of mode character indicating colors of the edges for the map ploting. It should be of same size than the number of regimes
#' @param col a character indicating the color of the vioplots
#' @param show.tip.label a logical indicating whether to show the tip labels on the phylogeny (defaults to TRUE, i.e. the labels are shown).
#' @param show.models a logical indicating whether to show details about model specification in the jive object.
#' @param direction a character string specifying the direction of the tree. Two values are possible: "rightwards" (the default) and "upwards".
#' @param trait.lab a charachter specifying the axis label for the traits
#' @param trait.lim a vector of mode numeric indicating the limits for trait ploting
#' @param srt.label an integer indicating the string rotation in degrees for the tip labels
#' @param c.reg a real number indicating where to plot the names of the regimes. The names are not plotted if c.reg == NULL.
#' @param tip.color the colours used for the tip labels, eventually recycled.
#' @param ... additional parameters that can be passed to \code{\link{vioplot}}
#' @export
#' @import phytools ape vioplot sm
#' @author Theo Gaboriau
#' @return plot 
#' @encoding UTF-8
#' @examples
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' 
#' colnames(Anolis_map) <- c("Hispaniola", "Cuba")
#' my.jive <- make_jive(Anolis_tree, Anolis_traits, 
#'  model.var="OU", model.mean="BM")
#' par(cex.lab = .8, cex.axis = .8, las = 1, mgp = c(2,0.5,0))
#' plot_jive(jive = my.jive, show.tip.label = TRUE, 
#' trait.lab = "Snout to vent length (cm)", srt.label = 0, c.reg = 2)
#' my.jive <- make_jive(Anolis_tree, Anolis_traits, Anolis_map,
#'  model.var=c("OU", "theta"), model.mean="BM")
#' par(cex.lab = .8, cex.axis = .8, las = 1, mgp = c(2,0.5,0))
#' plot_jive(jive = my.jive, show.tip.label = TRUE, c.reg = 2,
#'  trait.lab = "Snout to vent length (cm)", srt.label = 70, direction = "upwards")
#'  



plot_jive <- function(jive, col.map = NULL, col = "lightgrey", show.tip.label = T, show.models = T, direction = "rightwards",
                      trait.lab = "x", trait.lim = NULL, srt.label = 0, c.reg = NULL, tip.color = "#000000", ...){
  
  tree <- jive$data$tree
  map <- jive$data$map
  traits <- jive$data$traits
  
  st <- ncol(map$beta)
  tree <- map_to_simmap(tree, map)
  if(is.null(col.map)){
    col.map <- palette()[1:st]
    if (length(st) > 1) {
      cat("no colors provided. using the following legend:\n")
      print(col.map)
    }
  }
  if(st > 1) names(col.map) <- colnames(tree$mapped.edge)
  
  if(show.models){
    mm <- strsplit(jive$prior.mean$name, " ")[[1]]
    mv <- strsplit(jive$prior.var$name, " ")[[1]]
    ## Convert pars name into expressions to plot mathematical notations
    pars.m <- sapply(strsplit(names(jive$prior.mean$hprior), "[_]"), function(pr){
      if(pr[1] == "sigma^2") out <- ifelse("sigma" %in% mm, sprintf("sigma[m%s]^2", pr[2]), "sigma[m]^2")
      if(pr[1] == "alpha") out <- ifelse("alpha" %in% mm, sprintf("alpha[m%1$s]", pr[2]), "alpha[m]")
      if(pr[1] == "theta") out <- ifelse("theta" %in% mm, sprintf("theta[m%s]", pr[2]), "theta[m]")
      if(pr[1] == "root") out <- "theta[m0]"
      return(out)
    })
    pars.v <- sapply(strsplit(names(jive$prior.var$hprior), "[_]"), function(pr){
      if(pr[1] == "sigma^2") out <- ifelse("sigma" %in% mv, sprintf("sigma[v%s]^2", pr[2]), "sigma[v]^2")
      if(pr[1] == "alpha") out <- ifelse("alpha" %in% mv, sprintf("alpha[v%1$s]", pr[2]), "alpha[v]")
      if(pr[1] == "theta") out <- ifelse("theta" %in% mv, sprintf("theta[v%s]", pr[2]), "theta[v]")
      if(pr[1] == "root") out <- "theta[v0]"
      return(out)
    })
  }

  if(direction == "upwards"){
    if(show.models){
      par(mar = c(3,4,1,2))
    } else {
      par(mar = c(1,4,1,2))
    }
    root.len <- max(branching.times(tree))
    if(show.tip.label){
      ylim = c(0, 3*root.len)
    } else {
      ylim = c(0, 2*root.len)
    }
  
    if(srt.label == 0) srt.label = 90
    
    if(st > 1){
      plotSimmap(tree, direction = "upwards", ftype = "off", ylim = ylim, mar = par()$mar, colors = col.map, fsize = par("cex"))
      if(!is.null(c.reg)){
        text(x = rep(0, st), y = c.reg - seq(0, max(branching.times(tree))/10,length.out = st), labels = paste(1:st,jive$data$reg, sep = ":"),
             col = col.map, xpd = NA, cex = 1, adj = 0)
      }
    } else {
      plot(tree, direction = "upwards", y.lim = ylim, mar = par()$mar, edge.color = col.map, show.tip.label = F)
    }
    
    
    if(show.models){
      mtext(eval(parse(text = paste("substitute(a~~group('{', list(", paste(pars.m, collapse = ","),"), '}'), list(a=sprintf('Mean prior model: %s' , mm[1])))"))), side = 1, line = 0, at = 0, adj = 0, cex = 0.6)
      mtext(eval(parse(text = paste("substitute(a~~group('{', list(", paste(pars.v, collapse = ","),"), '}'), list(a=sprintf('LogVariance prior model: %s' , mv[1])))"))), side = 1, line = 1, at = 0, adj = 0, cex = 0.6)
    }
    
    init.usr <- par()$usr
    init.mar <- par()$mar
    par(new = TRUE, fig = c(0,1,ifelse(rep(show.tip.label,2),c(0.35,0.7),c(0.5,1))), bty = "n", xpd = NA)
    plot(c(1,length(tree$tip.label)), range(unlist(traits), na.rm = T), type = "n", xaxt = "n", ylab = trait.lab,
         xlab = "", ylim = ifelse(rep(!is.null(trait.lim),2), trait.lim, range(unlist(traits),na.rm = T)))
    pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    
    for(i in 1:length(tree$tip.label)){
      sp <- tree$tip.label[i]
      vioplot(traits[[sp]][!is.na(traits[[sp]])], add = T, at = pp$xx[i], col = col, ...)
    } 
    
    par(fig=c(0,1,0,1), usr = init.usr, mar = init.mar)
    if(show.tip.label){
      text(x = pp$xx[1:length(tree$tip.label)], y = 0.7*par()$usr[4],
           labels = gsub("_", " ", tree$tip.label), srt = srt.label,
           adj = 0, xpd = NA, col = tip.color)
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
    
    
    
    if(st > 1){
      plotSimmap(tree, direction = "rightwards", ftype = "off", xlim = xlim, mar = par()$mar, colors = col.map, fsize = par("cex"))
      if(!is.null(c.reg)){
        text(x = rep(0, st), y = c.reg - 0:(st-1), labels = paste(1:st,jive$data$reg, sep = ":"),
             col = col.map, xpd = NA, cex = 1, adj = 0)
      }
    } else {
      plot(tree, direction = "rightwards", x.lim = xlim, mar = par()$mar, edge.color = col.map, show.tip.label = F)
    }
    
    if(show.models){
      mtext(eval(parse(text = paste("substitute(a~~group('{', list(", paste(pars.m, collapse = ","),"), '}'), list(a=sprintf('Mean prior model: %s' , mm[1])))"))), line = 2, at = 0, adj = 0, cex = 0.6)
      mtext(eval(parse(text = paste("substitute(a~~group('{', list(", paste(pars.v, collapse = ","),"), '}'), list(a=sprintf('LogVariance prior model: %s' , mv[1])))"))), line = 1, at = 0, adj = 0, cex = 0.6)
    }
    
    init.usr <- par()$usr
    init.mar <- par()$mar
    par(new = TRUE, fig = c(ifelse(rep(show.tip.label,2),c(0.35,0.7),c(0.5,1)),0,1), bty = "n", xpd = NA)
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

