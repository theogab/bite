\name{marginal.lik}
\alias{marginal.lik}
\title{Marginal likelihood estimation}
\usage{
	marginal.lik(file = "jive_mcmc_OU1.log")
}
\arguments{
	\item{file}{log file recorded by \code{\link{jiveMCMC}} function}
}
\value{
	The Marginal likelihood of the model as estimated by thermodynamic integration .
}
\description{
	Calculate the marginal likelihood my thermodynamic integration.
}
\details{
This function estimates the marginal likelihood my thermodynamic integration (Lartillot and Philippe, 2006)
from the output of jiveMCMC.
}
\examples{
## Load test data
data(traitsOU1)
data(treeOU1)
my.jive <- jiveMake(treeOU1, traitsOU1,  model.var="OU1", model.mean="BM", model.lik="Multinorm")
jiveMCMC(my.jive, log.file="OU1_log.log", sampling.freq=10, print.freq=100, ngen=5000)
marginal.lik(log.file="OU1_log.log")
}
\author{
Simon Joly
}

