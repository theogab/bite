% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control_jive.R
\encoding{UTF-8}
\name{control_jive}
\alias{control_jive}
\title{Control tuning parameters of the jive algorithm}
\usage{
control_jive(jive, level = c("lik", "prior"), intvar = NULL,
  pars = NULL, window.size = NULL, initial.values = NULL,
  proposals = NULL, hyperprior = NULL, update.freq = NULL)
}
\arguments{
\item{jive}{a jive object obtained from \code{\link{make_jive}}}

\item{level}{character taken in c("lik", "prior") to specify on which level of the jive model, the control will operate (see details)}

\item{intvar}{character taken in names(jive$priors) giving the variable to be edited (see details)}

\item{pars}{vector of character taken in names(jive$priors[[intvar]]$init) giving the names of the hyper parameter to be edited}

\item{window.size}{initial window size for proposals during the mcmc algorithm. matrix or vector depending on the value of level and nreg (see details)}

\item{initial.values}{starting parameter values of the mcmc algorithm. matrix or vector depending on the value of level and nreg (see details)}

\item{proposals}{vector of characters taken in c("slidingWin", "slidingWinAbs", "logSlidingWinAbs","multiplierProposal", "multiplierProposalLakner","logNormal", "absNormal") to control proposal methods during mcmc algorithm (see details)}

\item{hyperprior}{list of hyperprior functions that can be generated with \code{\link{hpfun}}function. Ignored if level == "lik" (see details)}

\item{update.freq}{numeric giving the frequency at which parameters should be updated.}
}
\value{
A JIVE (of class "JIVE" and "list") object to parse into mcmc_bite function (see \code{\link{make_jive}})
}
\description{
This function modifies a jive object to tune the jive mcmc algorithm. The output will be different regarding which level of the jive model the user wants to tune ($lik, $priors). This function allows tuning of : initial window size for proposals, starting parameter value, proposal methods, Hyperpriors and update frequencies
}
\details{
If level == "lik" changes will be applied to the likelihood level of the algorithm. intvar is giving the variable on whic the changes will be operated
window.size and initial.values must be entered as a vector of length equal to the number of species. proposal must be a character

If level == "prior" changes will be applied to the prior level of the algorithm. intvar is giving the variable on which te change will be operated.
window.size and initial.values must be entered as a vector of size equal to the number of parameters or equal to the length of pars. 

Note that if you want to change the tuning at the three levels of the algorithm, you will have to use the control_jive function three times

proposals
Has to be one the following : "slidingWin" for Sliding window proposal unconstrained at maximum, "multiplierProposal", for multiplier proposal

Hyperprior
list of hyperpriror functions (see \code{\link{hpfun}}). User must provide a list of size equal to the number of parameters or equal to the length of pars
}
\examples{

data(Anolis_traits)
data(Anolis_tree)
 
## Create a jive object
my.jive <- make_jive(Anolis_tree, Anolis_traits[,-3],
model.priors = list(mean = "BM", logvar= c("OU", "root")))

## change starting values for the species means
my.jive$lik$init #default values
new.init <- rep(40,16)
my.jive <- control_jive(my.jive, level = "lik", intvar = "mean", initial.values = new.init)
my.jive$lik$init #mean initial values changed

 ## change hyperpriors for prior.mean
 plot_hp(my.jive) #default values
 new.hprior <- list(hpfun("Gamma", hp.pars = c(2,6)), hpfun("Uniform", c(20,80)))
 my.jive <- control_jive(my.jive, level = "prior", intvar = "mean", hyperprior = new.hprior)
 plot_hp(my.jive) #mean initial values changed
}
\author{
Theo Gaboriau
}
