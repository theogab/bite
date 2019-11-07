#' @title Plot trace and density from a log file
#' @description This function plots the trace and/or density of each mcmc sample.
#' @param mcmc.log Any mcmc sample with the saved iterations in rows and the variables in columns
#' @param type Character taken in c("trace", "density"). If both are specified, they are plotted side by side in the same graphical device
#' @param burnin The size of the burnin in number of iterations or the proportion of iteration you want to remove
#' @param variable The name or number of the variable to plot. If is.na(variable), all columns of mcmc.log will be plotted except "iter" and "temperature"
#' @param label Full variable name to be plotted
#' @param cex.est The magnification to be used for estimates display
#' @param kp.burn Logical specifying whether the plot window should adjust to the pre-burnin values (only evaluated if "trace" in type)
#' @param col,bty,... Other graphical parameters to parse to \code{\link[graphics]{par}}
#' @export
#' @import coda
#' @author Theo Gaboriau
#' @encoding UTF-8
#' @examples
#' 
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' 
#' ## Run a simple MCMC chain
#' my.jive <- make_jive(Anolis_tree, Anolis_traits,  model.var="OU", model.mean="BM")
#' mcmc_bite(my.jive, log.file="my.jive_mcmc.log", sampling.freq=10, print.freq=100, ngen=50000) 
#'
#' ## import the results in R
#' logfile <- "my.jive_mcmc.log"
#' res <- read.csv(logfile, header = TRUE, sep = "\t")
#' 
#' ## plot the results
#' plot_mcmc_bite(res, burnin = 0.2, variable = NA, cex.est = .7)
#' plot_mcmc_bite(res, burnin = 0.2, variable = "prior.mean", cex.est = .7)
#' 

plot_mcmc_bite <- function(mcmc.log, type = c("trace", "density"), burnin = 0, variable = NA,
                           label = NA, col = "#000000", cex.est = 1, bty = "n", kp.burn = FALSE,  ...){
  
  plot_func <- function(mcmc.log, type, burnin, variable, label, col, cex.est, bty, kp.burn,  ...){
    if(length(type) == 2){
      mrg <- par("mar")
      plot.new()
      par(fig=c(0,0.7,0,1), mar = c(mrg[1:3], 0.5))
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
      plot(mcmc.log[,"iter"], mcmc.log[,variable], type = "n", ylab = label, xlab = "Iterations",
           ylim = range(mcmc.log[if(kp.burn) burn|!burn else !burn,variable]), bty = bty, ...)
      for(t in unique(mcmc.log$temperature)){
        lines(mcmc.log[burn & mcmc.log$temperature == t,"iter"], mcmc.log[burn & mcmc.log$temperature == t,variable], col = adjustcolor(col, .5))
        lines(mcmc.log[!burn & mcmc.log$temperature == t,"iter"], mcmc.log[!burn & mcmc.log$temperature == t,variable], col = col)
      }
      mtext(sprintf("ESS = %s ; HPD [%s,%s] ; mean = %s", ess, round(hpd[1], 2), round(hpd[2], 2), round(mean(x[,variable]), 2)), 3, at = 0, adj = 0, cex = cex.est)
    }
    
    if("density" %in% type){
      if("trace" %in% type){
        par(fig=c(0.7,1,0,1), mar = c(mrg[1], 0.5, mrg[3:4]), new = TRUE)
        dens <- density(mcmc.log[!burn,variable])
        whpd <- dens$x >= hpd[,1] & dens$x <= hpd[,2]
        plot(dens$y, dens$x, type = "l", xlab = "Density", col = col, bty = bty, yaxt = "n", ylab = "",
             ylim = range(mcmc.log[if(kp.burn) burn|!burn else !burn,variable]), ...)
        polygon(c(0, dens$y[whpd], 0), c(hpd[,1], dens$x[whpd], hpd[,2]), col = adjustcolor(col, .5), border = NA)
        par(fig = c(0,1,0,1), mar = mrg)
      } else {
        dens <- density(mcmc.log[!burn,variable])
        whpd <- dens$x >= hpd[,1] & dens$x <= hpd[,2]
        plot(dens$x, dens$y, type = "l", las = 1, xlab = label, ylab = "Density", col = col, bty = bty, ...)
        polygon(c(hpd[,1], dens$x[whpd], hpd[,2]), c(0, dens$y[whpd], 0), col = adjustcolor(col, .5), border = NA)
        mtext(sprintf("HPD [%s,%s] ; mean = %s", round(hpd[1], 2), round(hpd[2], 2), round(mean(x[,variable]), 2)), 3, adj = 0, cex = cex.est)
      }
    }
  }
  
  if(!"temperature" %in% colnames(mcmc.log)) mcmc.log <- cbind(mcmc.log, temperature = rep(1, nrow(mcmc.log)))
  
  if(is.na(variable)){
    for(i in 1:ncol(mcmc.log)){
      var <- colnames(mcmc.log)[i]
      if(!var %in% c("iter", "temperature")){
        if(is.na(label)) lab <- var
        plot_func(mcmc.log, type = type, burnin = burnin, variable = var,
                  label = lab, col = col, cex.est = cex.est, bty = bty, kp.burn = kp.burn,  ...)
      }
    }
  } else {
    if(is.na(label)) label <- ifelse(is.numeric(variable), colnames(mcmc.log)[variable], variable)
    plot_func(mcmc.log, type = type, burnin = burnin, variable = variable,
              label = label, col = col, cex.est = cex.est, bty = bty, kp.burn = kp.burn,  ...)
  }
}

