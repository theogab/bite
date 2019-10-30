#' @title Plots estimates of species traits distribution
#' @description Density plot representing estimated species trait distributions under a jive model.
#' This function plots the mean or median density distribution and the HPD distributions assuming that the trait is normally distributed
#' @param phy phylogenetic tree provided as either a simmap or a phylo object
#' @param traits trait data used to perform the jive analysis. This has to be of the same form as the one used in \code{\link{make_jive}}
#' @param mcmc.log the output file of a \code{\link{mcmc_bite}} run
#' @param tip An integer giving the row corresponding to the species to be plotted
#' @param burnin The size of the burnin in number of iterations or the proportion of iteration you want to remove
#' @param legend A character giving the species name
#' @param conf A number of [0,1] giving the confidence level desired.
#' @param stat A character giving the function to be used to estimate species mean and variance from the posterior distributions. Must be one of be "mean" and "median"
#' @param trait.lab a charachter specifying the axis label for the traits
#' @param col color of the density filling. Must be of size two for estimates and HPD. If col and border are NULL, two random colors are assigned
#' @param border color of denstiy contour. Could be of size one or two. If border = NULL, no borders are plotted
#' @param lolipop a logical specifying wether the sample positions should be presented as lolipops
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
#' plot_pvo(Anolis_tree, Anolis_traits, res, 5)


plot_pvo <- function(phy, traits, mcmc.log, tip = 1, burnin = 0, conf = 0.95, stat = "median", trait.lab = "trait",
                          col = NULL, border = NULL, legend = F, lolipop = T,...){
  
  
  ### dealing with traits
  traits <- sapply(phy$tip.label, function(sp){
    if(any(traits[,1]  == sp)){
      traits[traits[,1]  == sp,2]
    } else {
      NA  
    }
  }, USE.NAMES = T, simplify = F) 
  
  if(is.numeric(tip)) label <- names(traits)[tip]
  else label <- tip
  
  if(burnin < 1) burnin <- burnin * nrow(mcmc.log)
  chain <- as.mcmc(mcmc.log[(burnin+1):nrow(mcmc.log),sprintf(c("%s_m", "%s_v"), label)])
  
  sam <- sample(1:nrow(chain), 5e6, replace = T)
  rhpd <- rnorm(1e6, chain[sam,1], sqrt(chain[sam,2]))
  hpd <- HPDinterval(as.mcmc(rhpd), prob = conf)
  if(stat == "median") mid <- apply(chain,2,median)
  else if(stat == "mean") mid <- apply(chain,2,mean)
  else stop(sprintf("%s: unknown stat"))
  
  rmid <- density(rnorm(1e6, mid[1], sqrt(mid[2])))
  rhpd <- density(rhpd[rhpd >= hpd[1] & rhpd <= hpd[2] ])
  
  label <- gsub("_", " ", label)
  plot(c(rhpd$x, rmid$x),c(rhpd$y, rmid$y), type = "n", main = label, ylab = "density", xlab = trait.lab)             
  
  if(is.null(border) & is.null(col)){
    col <- c("#2e86ab", "#001c55")
  } 
  if(is.null(border)){
    polygon(rhpd$x, rhpd$y, col = col[1], border = NA)
    polygon(rmid$x, rmid$y, col = col[2], border = NA)
  } else if(is.null(col)){
    lines(rhpd$x, rhpd$y, col = border[1])
    lines(rmid$x,rmid$y, col = border[2]) 
  } else {
    polygon(rhpd$x, rhpd$y, col = col[1], border = NA)
    lines(rhpd$x, rhpd$y, col = border[1])
    polygon(rmid$x, rmid$y, col = col[2], border = NA)
    lines(rmid$x,rmid$y, col = border[2]) 
  }
  
  if(lolipop){
    points(traits[[tip]],rep(max(rmid$y)/10,length(traits[[tip]])), pch = 16)
    for(i in 1:length(traits[[tip]])){
      lines(rep(traits[[tip]][i],2), c(0, max(rmid$y)/10))
    }
  } else {
    points(traits[[tip]],rep(par()$usr[3]+2/100*par()$usr[4],length(traits)), pch = 4)
  }
  
  
}
