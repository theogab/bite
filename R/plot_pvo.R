#' @title Plots estimates of species traits distribution
#' @description Density plot representing estimated species trait distributions under a jive model.
#' This function plots the mean or median density distribution and the HPD distributions assuming that the trait is normally distributed
#' @param phy phylogenetic tree provided as either a simmap or a phylo object
#' @param traits trait data used to perform the jive analysis. This has to be of the same form as the one used in \code{\link{make_jive}}
#' @param map map used to perform the jive analysis. This has to be of the same form as the one used in \code{\link{make_jive}}
#' @param mcmc.log the output file of a \code{\link{mcmc_bite}} run
#' @param tip A string giving the species to be plotted. If tip == NA, the posterior distribution of every tip is plotted along with the phylogenetic tree
#' @param burnin The size of the burnin in number of iterations or the proportion of iteration you want to remove
#' @param conf A number of [0,1] giving the confidence level desired.
#' @param stat A character giving the function to be used to estimate species mean and variance from the posterior distributions. Must be one of be "mean" and "median"
#' @param trait.lab a charachter specifying the axis label for the traits
#' @param col color of the density filling. Must be of size two for estimates and HPD. If col and border are NULL, two random colors are assigned
#' @param lab logical indicating whether to show species name in the plot. Only evaluated if tip =! NA
#' @param lolipop size and width of the lolipops representing samples
#' @param cex.tip size of the tips
#' @param var.f alternative distribution used to model intraspecific variation of the form function(n, pars). The function must return n samples from the given distribution.
#' @param ... Additional parameters that can be parsed to plot
#' @author Theo Gaboriau
#' @export
#' @examples
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' # Run a simple MCMC chain
#' my.jive <- make_jive(Anolis_tree, Anolis_traits[,-3],  model.priors=list(mean="BM", logvar = "OU"))
#' bite_ex <- tempdir()
#' logfile <- sprintf("%s/my.jive_mcmc.log", bite_ex)
#' mcmc_bite(my.jive, log.file=logfile, sampling.freq=1, print.freq=1, ngen=500) 
#' # import the results in R
#' res <- read.csv(logfile, header = TRUE, sep = "\t")
#'  plot_pvo(phy = Anolis_tree, traits = Anolis_traits, tip = NA, mcmc.log = res)
#'
#' @encoding UTF-8


plot_pvo <- function(phy, traits, map = NULL, mcmc.log, tip = NA, burnin = 0.1, conf = 0.95, stat = "median", trait.lab = "x",
                     col = NULL, lab = TRUE, lolipop = c(0.4, 0.4), cex.tip = par("cex"), var.f = NULL, ...){
  
  ### Main function
  plot_func <- function(traits, rmid, rhpd, label, trait.lab, col, lolipop, lab, plot.tree, minmax = NA, cy = 0, nreg =1){
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
  }, USE.NAMES = TRUE, simplify = FALSE)
  
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
    
    if(is.null(var.f)){
      var.f <- function(n, pars){
        rnorm(n, pars[,1], sqrt(exp(pars[,2])))
      }
    }
    
    ## get densities
    dens <- lapply(phy$tip.label, function(label){
      
      chain <- as.mcmc(mcmc.log[(burnin+1):nrow(mcmc.log),grepl(label, colnames(mcmc.log))])
      
      sam <- sample(1:nrow(chain), 1e4, replace = TRUE)
      rhpd <- var.f(1e4, chain[sam,])
      hpd <- HPDinterval(as.mcmc(rhpd), prob = conf)
      if(stat == "median") mid <- matrix(apply(chain,2,median), ncol = 2)
      else if(stat == "mean") mid <- matrix(apply(chain,2,mean), ncol = 2)
      else stop(sprintf("%s: unknown stat"))
      
      rmid <- density(var.f(1e4, mid))
      rhpd <- density(rhpd[rhpd >= hpd[1] & rhpd <= hpd[2]])
      return(list(rmid = cbind(x = rmid$x, y = rmid$y), rhpd = cbind(x = rhpd$x, y = rhpd$y)))
    })
    
    names(dens) <- phy$tip.label
    minmax <- apply(do.call(rbind, do.call(rbind, dens)), 2, range)
    n <- length(phy$tip.label)
    
    oldpar <- par(no.readonly = T)
    on.exit(par(oldpar))
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
      if(is.null(col)){
        col <- "lightgrey"
      }
      plot(phy, show.tip.label = FALSE, ...)
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
      plot_func(traits = traits, rmid = dens[[label]]$rmid, rhpd = dens[[label]]$rhpd, label = label,
                trait.lab = trait.lab, col = col.reg,lolipop =  lolipop, lab = TRUE, plot.tree = TRUE, 
                minmax = minmax,cy = cy, nreg = nreg)
      i <- i + 1
    }
    par(fig = c(0.7,1,0,1), mar = c(mrg[1],0,mrg[3:4]), new = TRUE)
    plot(0,xlim = c(0,1), ylim = c(1,n), yaxt = "n", ylab = "", bty = "n", xaxt = "n", xlab = "", type = "n")
    text(rep(0,n), pp$yy[1:16], gsub("_", " ", phy$tip.label), adj = 0, cex = cex.tip)
    par(mar = mrg, fig = c(0,1,0,1))
  } else {
    
    label <- tip
    
    chain <- as.mcmc(mcmc.log[(burnin+1):nrow(mcmc.log),grepl(label, colnames(mcmc.log))])
    
    sam <- sample(1:nrow(chain), 1e4, replace = TRUE)
    rhpd <- rnorm(1e4, chain[sam,1], sqrt(chain[sam,2]))
    hpd <- HPDinterval(as.mcmc(rhpd), prob = conf)
    if(stat == "median") mid <- apply(chain,2,median)
    else if(stat == "mean") mid <- apply(chain,2,mean)
    else stop(sprintf("%s: unknown stat"))
    
    rmid <- density(rnorm(1e4, mid[1], sqrt(mid[2])))
    rmid <- cbind(x = rmid$x, y = rmid$y)
    rhpd <- density(rhpd[rhpd >= hpd[1] & rhpd <= hpd[2]])
    rhpd <- cbind(x = rhpd$x, y = rhpd$y)
    if(is.null(col)){
      col <- "lightgrey"
    }
    plot_func(traits = traits, rmid = rmid, rhpd = rhpd, label = label, trait.lab = trait.lab,
              col = col, lolipop =lolipop, lab = TRUE, plot.tree = FALSE)
  }
  
}
