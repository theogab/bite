#' @title Plots estimates of species traits distribution
#' @description Density plot representing estimated species trait distributions under a jive model.
#' This function plots the mean or median density distribution and the HPD distributions assuming that the trait is normally distributed
#' @param traits trait data used to perform the jive analysis. This has to be of the same form as the one used in \code{\link{make_jive}}
#' @param mcmc.log the output file of a \code{\link{jive_mcmc}} run
#' @param tip An integer giving the row corresponding to the species to be plotted
#' @param label A character giving the species name
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
#' data(Anolis_map)
#' 


plot_pvo <- function(traits, mcmc.log, tip = 1, label = "tip 1", conf = 0.95, stat = "median", trait.lab = "trait",
                          col = NULL, border = NULL, legend = F, lolipop = T,...){
  
  chain <- as.mcmc(mcmc.log[,sprintf(c("%s_m", "%s_v"), rownames(traits)[tip])])
  
  hpd <- HPDinterval(chain)
  if(stat == "median") mid <- apply(chain,2,median)
  else if(stat == "mean") mid <- apply(chain,2,mean)
  else stop(sprintf("%s: unknown stat"))
  
  rmid <- density(rnorm(1e6, mid[1], sqrt(mid[2])))
  rhpd <- density(rnorm(1e6, sample(hpd[1,], 1e6, T), sqrt(sample(hpd[2,], 1e6, T))))
  
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
    points(traits[tip,],rep(max(rmid$y)/10,ncol(traits)), pch = 16)
    for(i in 1:ncol(traits)){
      lines(rep(traits[tip,i],2), c(0, max(rmid$y)/10))
    }
  } else {
    points(traits[tip,],rep(par()$usr[3]+2/100*par()$usr[4],ncol(traits)), pch = 4)
  }
  
  
}