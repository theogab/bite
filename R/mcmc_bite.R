#' @title MCMC algorithm
#' @description Implements Markov chain Monte Carlo sampling for trait evolution models
#' 
#' @details This function runs MCMC sampling on jive object \code{\link{make_jive}} or objects describing other models of the bite package.
#' The jive object contains both the dataset and set of model to be used in MCMC. This function implements both a conventional MCMC
#' and an MCMC with thermodynamic integration. The latter option is turned off by default and can be changed by
#' setting ncat to values > 1. The recommended ncat for TI is 10. When setting ncat > 1, make sure to specify burning.
#' As a rule of thumb set burning to 1/10 fraction of ngen. 
#' 
#' @param model an object of class "jive" or other objects from the bite package (see details)
#' @param log.file name of the output file that will store the log of MCMC chain
#' @param sampling.freq sampling frequency of the MCMC chain (how often chain will be saved into output file
#' @param print.freq printing frequency of the MCMC chain (how often chain will be printed in the R console). Setting it to 0 will suppress every printed message					
#' @param ncat number of classes for thermodynamic integration (see details)
#' @param beta.param beta value to define classes for thermodynamic integration (see details)
#' @param ngen number of generation in MCMC chain
#' @param burnin a burning phase of MCMC chain (has to be specified for thermodynamic integration)
#' @export
#' @author Theo Gaboriau, Anna Kostikova, Daniele Silvestro, and Simon Joly
#' @return Generates a log file in the users filespace at the path defined by log.file
#' @examples
#'  ## Load test data
#'  data(Anolis_traits)
#'  data(Anolis_tree)
#'  data(Anolis_map)
#'  
#'  ## Run a simple MCMC chain
#'  set.seed(300)
#'  my.jive <- make_jive(phy = Anolis_tree, traits = Anolis_traits[,-3],
#'   model.priors = list(mean = "BM", logvar= c("OU", "root")))
#'  bite_ex <- tempdir()
#'  logfile <- sprintf("%s/my.jive_mcmc.log", bite_ex)
#'  mcmc_bite(model = my.jive, log.file=logfile,
#'  sampling.freq=10, print.freq=0, ngen=1000) 
#'   
#'  ## Run an MCMC chain with thermodynamic integration
#'  logfile <- sprintf("%s/my.jive_mcmc_TI.log", bite_ex)
#'  mcmc_bite(my.jive, log.file=logfile, ncat=10, 
#'   sampling.freq=10, print.freq=100, ngen=1000, burnin=10) 
#' @encoding UTF-8


mcmc_bite <- function(model, log.file = "bite_mcmc.log", sampling.freq = 1000, print.freq = 1000, 
                      ncat = 1, beta.param = 0.3, ngen = 5000000, burnin = 0, sample.pars = T,
                      continue = F, ls.file = "ls_bite_mcmc.log")
{
  
  # General syntax
  # 0 : used for claculations
  # 1 : not used for calculations
  
  ## define the chain length for each category
  it <- ngen/ncat
  
  ## burnin
  if(burnin < 1) burnin <- burnin*it
  
  ## get the heating parameter for a chain - scaling classes
  if (ncat > 1) {
    beta.class <- heat_par(ncat, beta.param) 
  } else {
    beta.class <- 1
  }
  
  
  ## initial conditions
  if(print.freq > 0){
    cat("setting initial conditions\n")
  }
  
  # likelihood level
  if(continue){ ## get the previous log file to start from the last sampled point
    ## retrieve data from previous samples
    res <- read.table(log.file, header = T)
    if(sample.pars){
      last.line <- last.samp <- res[nrow(res),]
    } else {
      last.samp <- read.table(ls.file)
      names(last.samp) <- model$header
      last.line <- res[nrow(res),]
    }
    rm(res)
    
    ## Set conditions at the last sample
    pars.lik0 <- list()
    for(p in 1:length(model$priors)){
      pars.lik0[[p]] <- unlist(last.samp[grep(paste0(names(model$priors)[p], "_"), names(last.samp))])
      names(pars.lik0[[p]]) <-  gsub(paste0(names(model$priors)[p], "_"), "", names(pars.lik0[[p]]))
    }
    names(pars.lik0) <- names(model$lik$init)
  } else {
    ## Get defined initial conditions
    pars.lik0 <- model$lik$init
  }
  
  lik0 <- model$lik$model(pars.lik0, model$data$traits, model$data$counts)
  
  # prior level
  pars.priors0 <- list()
  priors0 <- c()
  hpriors0 <- list()
  
  
  if(continue){
    iter <-last.line$iter 
    for(p in 1:length(model$priors)){
      priors0[p] <- unlist(last.line[paste("prior",names(model$priors)[p], sep = ".")])
      pars.priors0[[p]] <- unlist(last.line[grep(paste0(names(model$priors)[p], "."), names(last.line))])
      names(pars.priors0[[p]]) <- names(model$priors[[p]]$init)
      hpriors0[[p]] <- unlist(mapply(do.call, model$priors[[p]]$hprior, lapply(pars.priors0[[p]], list))[1,])
    }
  } else {
    for(p in 1:length(model$priors)){
      pars.priors0[[p]] <- model$priors[[p]]$init
      priors0[p] <- model$priors[[p]]$value
      hpriors0[[p]] <- unlist(mapply(do.call, model$priors[[p]]$hprior, lapply(pars.priors0[[p]], list))[1,])
    }
    iter <- 1
  }
  
  
  # mcmc parameters
  if(print.freq > 0){
    cat("generation\tposterior\n")
    if(!continue){
      if(sample.pars){
        cat(paste(model$header, collapse = "\t"), "\n", append = FALSE, file = log.file)
      } else {
        cat(paste(model$header[-c(grep("mean_", model$header), grep("logvar_", model$header))], collapse = "\t"), "\n", append = FALSE, file = log.file)
      }
    }
  }
  
  if(continue){
    bet <- last.line$temperature
    it.beta <- which(abs(beta.class - bet) < 1e-6)
  } else {
    it.beta <- 1
    bet <- beta.class[it.beta]
  }
  
  if(ncat > 1) cat("beta = ", bet, "\n")
  # making sure the update frequencies sum to 1
  update.freq <- c(model$lik$update.freq, sapply(model$priors, function(x) x$update.freq))
  update.freq <- cumsum(update.freq/sum(update.freq))
  proposals <- c(0,0,0) # 1st number: update means and variance; 2nd: update mean priors, 3rd: update variance prior
  proposals.accepted <- c(0,0,0)
  
  # posterior level
  post0 <- (sum(lik0) + sum(priors0 * bet) + sum(unlist(hpriors0)))
  

  ## Start iterations
  for (i in iter:(it*ncat)) {
    
    # test whether we update parameters, or hyper parameters
    r <- min(which(runif(1) <= update.freq))
    proposals[r] <- proposals[r] + 1
    
    if (r == 1) # update lik.pars, calculate new lik
    { 
      ind <- sample(1:model$data$n, model$lik$n.u, replace = FALSE)
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      lik1 <- lik0
      pars.lik1 <- pars.lik0
      priors1 <- priors0
      hasting.ratio <- 0
      
      for(p in 1:length(model$priors)){
        tmp <- model$lik$prop[[p]](i = pars.lik0[[p]][ind], d = model$lik$ws[[p]][ind], u)
        pars.lik0[[p]][ind] <- tmp$v
        hasting.ratio <- hasting.ratio + tmp$lnHastingsRatio
        priors0[[p]] <- model$priors[[p]]$model(x = pars.lik0[[p]], n = model$data$n, pars = pars.priors0[[p]],
                                                Pi = model$priors[[p]]$Pi, par.n = 0, 
                                                data = model$priors[[p]]$data, map = model$priors[[p]]$map)$loglik
      }
      
      lik0 <- model$lik$model(pars.lik0, model$data$traits, model$data$counts)
      
    } else { # update priors.pars calculate new priors and new hpriors
      p <- r - 1
      par.n <- sample(1:length(model$priors[[p]]$prop), 1) # one parameter at a time for now (update simultaneously sigmas/thetas for multiple regimes?)
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      pars.priors1 <- pars.priors0 # ancient parameter values are kept
      priors1 <- priors0 # ancient prior value is kept
      hpriors1 <- hpriors0 # ancient hprior value is kept
      tmp <- model$priors[[p]]$prop[[par.n]](i = pars.priors0[[p]][par.n], d = model$priors[[p]]$ws[par.n], u) #update with proposal function
      pars.priors0[[p]][par.n] <- tmp$v
      
      # calculate new var/covar and expectation
      mat1 <- model$priors[[p]]$data
      mat0 <- try(model$priors[[p]]$model(x = pars.lik0[[p]], n = model$data$n, pars = pars.priors0[[p]],
                                              Pi = model$priors[[p]]$Pi, par.n = par.n, 
                                              data = model$priors[[p]]$data,
                                              map = model$prior[[p]]$map), silent = TRUE)
      if(any(grepl("Error", mat0))){
        priors0[p] <- -Inf
      } else {
        model$priors[[p]]$data <- mat0$data
        # calculate prior and hprior
        priors0[p] <- mat0$loglik
        hpriors0[[p]] <- unlist(mapply(do.call, model$priors[[p]]$hprior, lapply(pars.priors0[[p]], list))[1,])
      }
      
      hasting.ratio <- tmp$lnHastingsRatio
      
    }
    
    # Posterior calculation
    post1 <- post0
    post0 <- (sum(lik0) + sum(priors0 * bet) + sum(unlist(hpriors0)))
    
    # acceptance probability (log scale)
    if(any(is.infinite(c(lik0, priors0, unlist(hpriors0))))){
      pr <- -Inf
    } else {
      pr <- post0 - post1 + hasting.ratio
    }
    
    if (pr >= log(runif(1))){ # count acceptance
      proposals.accepted[r] <- proposals.accepted[r] + 1
    } else {# cancel changes
      post0 <- post1
      if (r == 1){
        pars.lik0 <- pars.lik1
        lik0 <- lik1
        priors0 <- priors1
      } else {
        pars.priors0 <- pars.priors1
        priors0 <- priors1
        hpriors0 <- hpriors1
        model$priors[[p]]$data <- mat1
      }
    }
    
    # log to file with frequency sampling.freq
    if (i %% sampling.freq == 0 & i >= burnin) {
      if(sample.pars){
        cat(paste(c(i, post0, sum(lik0), priors0, unlist(sapply(1:length(model$priors), function(p) c(pars.priors0[[p]], pars.lik0[[p]]))), sum(proposals.accepted)/i, bet), collapse = "\t"), "\n",
            append=TRUE, file=log.file) 
      } else {
        cat(paste(c(i, post0, sum(lik0), priors0, unlist(sapply(1:length(model$priors), function(p) pars.priors0[[p]])), sum(proposals.accepted)/i, bet), collapse = "\t"), "\n",
            append=TRUE, file=log.file)
        cat(paste(c(i, post0, sum(lik0), priors0, unlist(sapply(1:length(model$priors), function(p) c(pars.priors0[[p]], pars.lik0[[p]]))), sum(proposals.accepted)/i, bet), collapse = "\t"), "\n",
            append=FALSE, file=ls.file) 
      }
    }
    
    # Print to screen
    if(print.freq > 0){
      if (i %% print.freq == 0) {
        cat(i,'\t',post0,'\n') 
      }
    }
    
    # change beta value if the length of the category is reach 
    if(i%%it == 0 & i < ngen){
      it.beta = it.beta+1
      bet <- beta.class[it.beta]
      cat("beta = ", bet, "\n")
    }
    
  } # end of for 
  
  if(print.freq > 0){
   # Calculate acceptance rate
    acceptance.results <- proposals.accepted / proposals
    names(acceptance.results) <- names(proposals) <- c("Likelihood parameters",sprintf("prior.%s",names(model$priors)))
    cat("\nEffective proposal frequency\n")
    print(proposals/ngen)
    cat("\nAcceptance ratios\n")
    print(acceptance.results) 
  }
  
  
}



