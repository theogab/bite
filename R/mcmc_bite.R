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
#' @param print.freq printing frequency of the MCMC chain (how often chain will be printed in the R console)					
#' @param ncat number of classes for thermodynamic integration (see details)
#' @param beta.param beta value to define classes for thermodynamic integration (see details)
#' @param ngen number of generation in MCMC chain
#' @param burnin a burning phase of MCMC chain (has to be specified for thermodynamic integration)
#' @export
#' @author Theo Gaboriau, Anna Kostikova, Daniele Silvestro, and Simon Joly
#' @return none
#' @examples
#'
#' ## Load test data
#' data(Anolis_traits)
#' data(Anolis_tree)
#' data(Anolis_map)
#' 
#' ## Run a simple MCMC chain
#' set.seed(300)
#' my.jive <- make_jive(Anolis_tree, Anolis_traits, map = Anolis_map,
#'  model.var=c("OU", "root", "theta", "alpha"), model.mean="BM")
#' mcmc_bite(my.jive, log.file="my.jive_MCMC.log",
#'  sampling.freq=10, print.freq=100, ngen=5000) 
#'
#'  my.jive <- make_jive(Anolis_tree, Anolis_traits, map = Anolis_map,
#'   model.var=c("OU", "root", "theta"), model.mean="BM")
#' mcmc_bite(my.jive, log.file="my.jive_MCMC.log",
#'  sampling.freq=10, print.freq=100, ngen=5000) 
#'
#' ## Run an MCMC chain with thermodynamic integration
#' mcmc_bite(my.jive, log.file="my.jive_MCMC_TI.log", ncat=10, 
#' sampling.freq=10, print.freq=100, ngen=5000, burnin=500) 
#' 
#' @encoding UTF-8


mcmc_bite <- function(model, log.file = "bite_mcmc.log", sampling.freq = 1000, print.freq = 1000, 
                      ncat = 1, beta.param = 0.3, ngen = 5000000, burnin = 0)
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
  cat("setting initial conditions\n")
  # likelihood level
  m.sp0 <- model$lik$init$m.sp
  v.sp0 <- model$lik$init$v.sp
  lik0 <- model$lik$model(m.sp0, v.sp0, model$data$traits, model$data$counts)
  
  # prior mean level
  pars.m0 <- do.call(c, model$prior.mean$init)
  prior.mean0 <- calc_prior(n = model$data$n, mat = model$prior.mean$data, x = m.sp0)
  hprior.mean0 <- unlist(mapply(do.call, model$prior.mean$hprior, lapply(pars.m0, list))[1,])
  
  # prior var level
  pars.v0 <- do.call(c, model$prior.var$init)
  prior.var0 <- calc_prior(n = model$data$n, mat = model$prior.var$data, x = log(v.sp0))
  hprior.var0 <- unlist(mapply(do.call, model$prior.var$hprior, lapply(pars.v0, list))[1,])
  
  # mcmc parameters
  cat("generation\tposterior\n")
  cat(paste(model$header, collapse = "\t"), "\n", append = FALSE, file = log.file)
  it.beta <- 1
  bet <- beta.class[it.beta]
  if(ncat > 1) cat("beta = ", bet, "\n")
  # making sure the update frequencies sum to 1
  update.freq <- c(model$lik$update.freq, model$prior.mean$update.freq, model$prior.var$update.freq)
  update.freq <- cumsum(update.freq/sum(update.freq))
  proposals <- c(0,0,0) # 1st number: update means and variance; 2nd: update mean priors, 3rd: update variance prior
  proposals.accepted <- c(0,0,0)
  
  # posterior level
  post0 <- (sum(lik0) + (prior.mean0 * bet) + (prior.var0 * bet) + sum(c(hprior.mean0, hprior.var0)))
  
  
  ## Start iterations
  for (i in 1:(it*ncat)) {
    
    # test whether we update lik (r[1]), prior_mean (r[2]) or prior_var (r[3])
    r <- runif(1) <= update.freq
    r[2] <- r[2] & !r[1]
    r[3] <- r[3] & !r[2] & !r[1]
    proposals[r] <- proposals[r] + 1
    
    if (r[1]) # update m.sp, v.sp, calculate new lik
    { 
      ind <- sample(1:length(m.sp0), 1, replace = F) # 5 is just for now number of species updated
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      m.sp1 <- m.sp0
      v.sp1 <- v.sp0
      lik1 <- lik0
      prior.mean1 <- prior.mean0
      prior.var1 <- prior.var0
      
      tmp <- model$lik$prop$m.sp(i = m.sp0[ind], d = model$lik$ws$m.sp[ind], u)
      m.sp0[ind] <- tmp$v
      prior.mean0 <- calc_prior(n = model$data$n, mat = model$prior.mean$data, x = m.sp0)
        
      tmp <- model$lik$prop$v.sp(i = v.sp0[ind], d = model$lik$ws$v.sp[ind], u)
      v.sp0[ind] <- tmp$v
      prior.var0 <- calc_prior(n = model$data$n, mat = model$prior.var$data, x = log(v.sp0))
      
      lik0 <- model$lik$model(m.sp0, v.sp0, model$data$traits, model$data$counts)
      hasting.ratio <- sum(tmp$lnHastingsRatio)
      
    }
    
    if (r[2]) # update parameters of the mean prior, calculate new prior_mean and new hprior_mean
    {
      par.n <- sample(1:length(model$prior.mean$prop), 1) # one parameter at a time for now (update simultaneously sigmas/thetas for multiple regimes?)
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      pars.m1 <- pars.m0 # ancient parameter values are kept
      prior.mean1 <- prior.mean0 # ancient prior_mean value is kept
      hprior.mean1 <- hprior.mean0 # ancient hprior_mean value is kept
      tmp <- model$prior.mean$prop[[par.n]](i = pars.m0[par.n], d = do.call(c,model$prior.mean$ws)[par.n], u) #update with proposal function
      pars.m0[par.n] <- tmp$v
      # calculate new var/covar and expectation
      mat.mean1 <- model$prior.mean$data
      mat.mean0 <- try(model$prior.mean$model(n = model$data$n, n.p = par.n,
                                           pars = pars.v0, tree = model$data$tree,
                                           map = model$prior.mean$map, t.vcv = model$data$vcv, nr = model$prior.mean$nr), silent = T)
      if(any(grepl("Error", mat.mean0))){
        prior.mean0 <- Inf
      } else {
        model$prior.mean$data <- lapply(1:3, function(k) if(mat.mean0[[k]][[1]]) mat.mean0[[k]][[2]] else mat.mean1[[k]]) # keep only updated ones
        # calculate prior and hprior
        prior.mean0 <- calc_prior(n = model$data$n, mat = model$prior.mean$data, x = m.sp0)
        hprior.mean0 <- unlist(mapply(do.call, model$prior.mean$hprior, lapply(pars.m0, list))[1,])
      }
      hasting.ratio <- tmp$lnHastingsRatio
    } 
    
    if (r[3]) # update parameters of the var prior, calculate new prior_var and new hprior_var
    {
      par.n <- sample(1:length(model$prior.var$prop), 1) # one parameter at a time for now (update simultaneously sigmas/thetas for multiple regimes?)
      u = runif(1) # parameter of the multiplier proposal (ignored if proposal == "SlidingWin")
      pars.v1 <- pars.v0 # ancient parameter values are kept
      prior.var1 <- prior.var0 # ancient prior_var value is kept
      hprior.var1 <- hprior.var0 # ancient hprior_var value is kept
      tmp <- model$prior.var$prop[[par.n]](i = pars.v0[par.n], d = do.call(c, model$prior.var$ws)[par.n], u) #update with proposal function
      pars.v0[par.n] <- tmp$v 
      # calculate new var/covar and expectation
      mat.var1 <- model$prior.var$data
      mat.var0 <- try(model$prior.var$model(n = model$data$n, n.p = par.n,
                                       pars = pars.v0, tree = model$data$tree,
                                       map = model$prior.var$map, t.vcv = model$data$vcv, nr = model$prior.var$nr), silent = T)
      if(any(grepl("Error", mat.var0))){
        prior.var0 <- Inf
      } else {
        model$prior.var$data <- lapply(1:3, function(k) if(mat.var0[[k]][[1]]) mat.var0[[k]][[2]] else mat.var1[[k]]) # keep only updated ones
        # calculate prior and hprior
        prior.var0 <- calc_prior(n = model$data$n, mat = model$prior.var$data, x = log(v.sp0))
        hprior.var0 <- unlist(mapply(do.call, model$prior.var$hprior, lapply(pars.v0, list))[1,])
      }
      hasting.ratio <- tmp$lnHastingsRatio
    }
    
    # Posterior calculation
    post1 <- post0
    post0 <- (sum(lik0) + (prior.mean0 * bet) + (prior.var0 * bet) + sum(c(hprior.mean0, hprior.var0)))
    
    # acceptance probability (log scale)
    if(any(is.infinite(c(lik0, prior.mean0, prior.var0, hprior.mean0, hprior.var0)))){
      pr <- -Inf
    } else {
      pr <- post0 - post1 + hasting.ratio
    }
    
    if (pr >= log(runif(1))) # count acceptance
    {
      proposals.accepted[r] <- proposals.accepted[r] + 1
    } else # cancel changes
    {
      post0 <- post1
      
      if (r[1]){
        m.sp0 <- m.sp1
        v.sp0 <- v.sp1
        lik0 <- lik1
        prior.mean0 <- prior.mean1
        prior.var0 <- prior.var1
      }
      
      if (r[2]){
        pars.m0 <- pars.m1
        prior.mean0 <- prior.mean1
        hprior.mean0 <- hprior.mean1
        model$prior.mean$data <- mat.mean1
      }
      
      if (r[3]){
        pars.v0 <- pars.v1
        prior.var0 <- prior.var1
        hprior.var0 <- hprior.var1
        model$prior.var$data <- mat.var1
      }
    }
    
    # log to file with frequency sampling.freq
    if (i %% sampling.freq == 0 & i >= burnin) {
      cat(paste(c(i, post0, sum(lik0), prior.mean0, prior.var0, pars.m0, pars.v0, m.sp0, v.sp0, sum(proposals.accepted)/i, bet), collapse = "\t"), "\n",
       append=TRUE, file=log.file) 
    }
    
    # Print to screen
    if (i %% print.freq == 0) {
      cat(i,'\t',post0,'\n') 
    }
    
    # change beta value if the length of the category is reach 
    if(i%%it == 0 & i < ngen){
      it.beta = it.beta+1
      bet <- beta.class[it.beta]
      cat("beta = ", bet, "\n")
    }
    
  } # end of for 
  
  # Calculate acceptance rate
  names(proposals) <- c("means & variances","prior.means","prior.variances")
  cat("\nEffective proposal frequency\n")
  print(proposals/ngen)
  acceptance.results <- proposals.accepted / proposals
  names(acceptance.results) <- c("means & variances","prior.means","prior.variances")
  cat("\nAcceptance ratios\n")
  print(acceptance.results)
  
}



