#' @title Plots estimates of species traits distribution
#' @description Density plot representing estimated species trait distributions under a jive model.
#' This function plots the mean or median density distribution and the HPD distributions assuming that the trait is normally distributed
#' @param phy phylogenetic tree provided as either a simmap or a phylo object
#' @param traits trait data used to perform the jive analysis. This has to be of the same form as the one used in \code{\link{make_jive}}
#' @param map map used to perform the jive analysis. This has to be of the same form as the one used in \code{\link{make_jive}}
#' @param mcmc.log the output file of a \code{\link{mcmc_bite}} run
#' @param tip An integer giving the row corresponding to the species to be plotted. If tip == NA, the posterior distribution of every tip is plotted along with the phylogenetic tree
#' @param burnin The size of the burnin in number of iterations or the proportion of iteration you want to remove
#' @param conf A number of [0,1] giving the confidence level desired.
#' @param stat A character giving the function to be used to estimate species mean and variance from the posterior distributions. Must be one of be "mean" and "median"
#' @param trait.lab a charachter specifying the axis label for the traits
#' @param col color of the density filling. Must be of size two for estimates and HPD. If col and border are NULL, two random colors are assigned
#' @param border color of denstiy contour. Could be of size one or two. If border = NULL, no borders are plotted
#' @param lab logical indicating whether to show species name in the plot. Only evaluated if tip =! NA
#' @param lolipop a logical specifying wether the sample positions should be presented as lolipops
#' @param cex.tip size of the tips
#' 
#' @author Theo Gaboriau
#' @export
#' @examples
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#' 
#' # Run a simple MCMC chain
#' my.jive <- make_jive(Anolis_tree, Anolis_traits,  model.var="OU", model.mean="BM")
#' mcmc_bite(my.jive, log.file="my.jive_mcmc.log", sampling.freq=10, print.freq=100, ngen=50000) 
#'
#' ## import the results in R
#' logfile <- "my.jive_mcmc.log"
#' res <- read.csv(logfile, header = TRUE, sep = "\t")
#' 
#' plot_pvo(Anolis_tree, Anolis_traits, mcmc.log = res)
#' 
#' @encoding UTF-8


plot_pvo <- function(phy, traits, map = NULL, mcmc.log, tip = NA, burnin = 0.1, conf = 0.95, stat = "median", trait.lab = "x",
                     col = NULL, lab = T, lolipop = c(0.4, 0.4), cex.tip = par("cex"), ...){
  
  ### Main function
  plot_func <- function(traits, rmid, rhpd, label, trait.lab, col, lolipop, lab, plot.tree, minmax = NA, cy = 0, nreg){
    if(!plot.tree){
      plot(c(rhpd[,"x"], rmid[,"x"]),c(rhpd[,"y"], rmid[,"y"]), type = "n", main = ifelse(lab, gsub("_", " ", label), ""), 
           ylab = "density", xlab = trait.lab)      
    } 
    
    col <- c(col, adjustcolor(col, .5))
    polygon(rhpd[,"x"], rhpd[,"y"]+cy, col = col[1], border = NA)
    polygon(rmid[,"x"], rmid[,"y"]+cy, col = col[2], border = NA)
    
    points(traits[[label]]+cy,rep(max(rmid[,"y"])/10+cy,length(traits[[label]])), pch = 16, cex = lolipop[1])
    for(i in 1:length(traits[[label]])){
      lines(rep(traits[[label]][i],2)+cy, c(0, max(rmid[,"y"])/10)+cy, lwd = lolipop[2])
    }
  }
  
  ### dealing with traits
  traits <- sapply(phy$tip.label, function(sp){
    if(any(traits[,1]  == sp)){
      traits[traits[,1]  == sp,2]
    } else {
      NA  
    }
  }, USE.NAMES = T, simplify = F)
  
  if(burnin < 1) burnin <- burnin * nrow(mcmc.log)
  
  if(is.na(tip)){
    
    ## map info
    if(!is.null(map)){
      reg.names <- colnames(map)
      map <- input_to_map(phy, map = map)
      phy <- map_to_simmap(phy, map)
    }
    if(!is.null(phy$nodelabels)){
      map <- input_to_map(phy, ndlabels = phy$nodelabels)
      reg.names <- unique(phy$nodelabels)
      phy <- map_to_simmap(phy, map)
    }
    
    
    ## get densities
    dens <- lapply(phy$tip.label, function(label){
      
      chain <- as.mcmc(mcmc.log[(burnin+1):nrow(mcmc.log),sprintf(c("%s_m", "%s_v"), label)])
      
      sam <- sample(1:nrow(chain), 1e6, replace = T)
      rhpd <- rnorm(1e6, chain[sam,1], sqrt(chain[sam,2]))
      hpd <- HPDinterval(as.mcmc(rhpd), prob = conf)
      if(stat == "median") mid <- apply(chain,2,median)
      else if(stat == "mean") mid <- apply(chain,2,mean)
      else stop(sprintf("%s: unknown stat"))
      
      rmid <- density(rnorm(1e6, mid[1], sqrt(mid[2])))
      rhpd <- density(rhpd[rhpd >= hpd[1] & rhpd <= hpd[2]])
      return(list(rmid = cbind(x = rmid$x, y = rmid$y), rhpd = cbind(x = rhpd$x, y = rhpd$y)))
    })
    
    names(dens) <- phy$tip.label
    minmax <- apply(do.call(rbind, do.call(rbind, dens)), 2, range)
    n <- length(phy$tip.label)
    
    mrg <- par("mar")
    par(fig = c(0, 0.4, 0, 1), mar = c(mrg[1:3],0))
    if("simmap" %in% class(phy)){
      nreg <- ncol(phy$mapped.edge)
      if(is.null(col)){
        col <- setNames(palette()[1:nreg+1], colnames(phy$mapped.edge))
      } else {
        names(col) <- colnames(phy$mapped.edge)
      }
      plotSimmap(phy, ftype = "off", colors = col, mar = par("mar"))
    } else {
      plot(phy, show.tip.label = F, ...)
      nreg <- 1
    }
    pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    par(fig = c(0.4,0.7,0,1), mar = c(mrg[1],0,mrg[3],0), new = TRUE, xpd = NA)
    plot(0,xlim = minmax[,1], ylim = c(minmax[1,2], minmax[2,2]*12/10*(n-1)),
         yaxt = "n", xlab = trait.lab, ylab = "", bty = "n", type = "n", ...)
    i <- 1
    for(label in phy$tip.label){
      cy <- (pp$yy[i]-1)*minmax[2,2]*12/10
      lines(c(minmax[1,1]-minmax[1,1]*1/4, minmax[2,1]), rep(cy,2), lty = 2, col = "#7a6563", xpd = NA)
      if("simmap" %in% class(phy)){
        reg <- colnames(phy$mapped.edge)
        br <- phy$maps[[which(phy$edge[,2] == which(phy$tip.label == label))]]
        col.reg <- col[reg %in% names(br)[length(br)]]
      } else {
        col.reg <- col
      }
      plot_func(traits, dens[[label]]$rmid, dens[[label]]$rhpd, label, trait.lab, col.reg, lolipop, lab = T, plot.tree = T, minmax,cy, nreg)
      i <- i + 1
    }
    par(fig = c(0.7,1,0,1), mar = c(mrg[1],0,mrg[3:4]), new = TRUE)
    plot(0,xlim = c(0,1), ylim = c(1,n), yaxt = "n", ylab = "", bty = "n", xaxt = "n", xlab = "", type = "n")
    text(rep(0,n), pp$yy[1:16], gsub("_", " ", phy$tip.label), adj = 0, cex = cex.tip)
    par(mar = mrg, fig = c(0,1,0,1))
  } else {
    
    if(is.numeric(tip)) label <- names(traits)[tip]
    else label <- tip
    
    chain <- as.mcmc(mcmc.log[(burnin+1):nrow(mcmc.log),sprintf(c("%s_m", "%s_v"), label)])
    
    sam <- sample(1:nrow(chain), 5e6, replace = T)
    rhpd <- rnorm(1e6, chain[sam,1], sqrt(chain[sam,2]))
    hpd <- HPDinterval(as.mcmc(rhpd), prob = conf)
    if(stat == "median") mid <- apply(chain,2,median)
    else if(stat == "mean") mid <- apply(chain,2,mean)
    else stop(sprintf("%s: unknown stat"))
    
    rmid <- density(rnorm(1e6, mid[1], sqrt(mid[2])))
    rhpd <- density(rhpd[rhpd >= hpd[1] & rhpd <= hpd[2]])
    plot_func(traits, rmid, rhpd, tip, burnin, conf, stat, trait.lab, col, border, legend, lolipop, show.lab = T,1)
  }
  
}
