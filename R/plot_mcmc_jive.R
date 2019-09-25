#' @title Plot \code{\link{jive_mcmc}} results
#' @description This function plots the trace and/or density of each mcmc sample
#' @details  
#' @param mcmc.log Any mcmc sample with the saved iterations in rows and the variables in columns
#' @param type character taken in c("trace", "density"). If both are specified, they are plotted side by side in the same graphical device
#' @param burnin The size of the burnin in number of iterations or the proportion of iteration you want to remove
#' @param variable The name or number of the variable to plot
#' @param label Full variable name to be plotted
#' @param ... Other graphical parameters to parse to \code{\link{par()}}
#' @export
#' @import coda
#' @author Theo Gaboriau
#' 
#' @examples
#' 
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' 
#' ## Run a simple MCMC chain
#' my.jive <- make_jive(Anolis_tree, Anolis_traits,  model.var="OU", model.mean="BM", root.station = T)
#' mcmc_jive(my.jive, log.file="my.jive_mcmc.log", sampling.freq=10, print.freq=100, ngen=50000) 
#'
#' ## import the results in R
#' logfile <- "my.jive_mcmc.log"
#' res <- read.csv(logfile, header = T, sep = "\t")
#' 
#' ## plot the results
#' plot_mcmc_jive(res, burnin = 0.2, variable = "prior.mean", cex.est = .7)
#' 

plot_mcmc_jive <- function(mcmc.log, type = c("trace", "density"), burnin = 0, variable = "log.lik",
                           label = ifelse(is.numeric(variable),colnames(mcmc.log)[variable], variable), col = "#000000", cex.est = 1, bty = "n", ...){
  
  par(mar =c(5,4,2,1))
  if(length(type) == 2){
    par(mfrow = c(1,2))
  }
  
  if(is.numeric(variable)) variable <- colnames(mcmc.log)[variable]
  
  if(burnin>0){
    temp <- unique(mcmc.log$temperature)
    if(burnin < 1){
      burnin <- sapply(temp, function(t) quantile(mcmc.log$iter[mcmc.log$temperature == t], burnin))
    } else {
      burnin <- sapply(temp, function(t) min(mcmc.log$iter[mcmc.log$temperature == t])) + burnin
    }
    burn <- unlist(lapply(1:length(temp), function(t) mcmc.log[mcmc.log$temperature == temp[t],"iter"] <= burnin[t]))
  } else {
    burn <- rep(F, nrow(mcmc.log))
  }
  
  x <- as.mcmc(mcmc.log[!burn,])
  ess <- round(effectiveSize(x[,variable]),2)
  hpd <- HPDinterval(x[,variable])
  
  if("trace" %in% type){
    plot(mcmc.log[,"iter"], mcmc.log[,variable], type = "n", ylab = label, xlab = "Iterations", las = 1, ylim = range(mcmc.log[!burn,variable]), bty = bty, ...)
    for(t in unique(mcmc.log$temperature)){
      lines(mcmc.log[burn & mcmc.log$temperature == t,"iter"], mcmc.log[burn & mcmc.log$temperature == t,variable], col = adjustcolor(col, .5))
      lines(mcmc.log[!burn & mcmc.log$temperature == t,"iter"], mcmc.log[!burn & mcmc.log$temperature == t,variable], col = col)
    }
    mtext(sprintf("Estimated sample size: ESS = %s", ess), 3, at = 0, adj = 0, cex = cex.est)
  }

  if("density" %in% type){
    dens <- density(mcmc.log[!burn,variable])
    whpd <- dens$x >= hpd[,1] & dens$x <= hpd[,2]
    plot(dens$x, dens$y, type = "l", las = 1, xlab = label, ylab = "Density", col = col, bty = bty, ...)
    polygon(c(hpd[,1], dens$x[whpd], hpd[,2]), c(0, dens$y[whpd], 0), col = adjustcolor(col, .5), border = NA)
    mtext(sprintf("HPD [%s,%s] ; mean = %s", round(hpd[1], 2), round(hpd[2], 2), round(mean(x[,variable]), 2)), 3, adj = 0, cex = cex.est)
  }

}

