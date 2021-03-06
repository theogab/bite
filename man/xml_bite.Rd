% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xml_bite.R
\encoding{UTF-8}
\name{xml_bite}
\alias{xml_bite}
\title{Write xml file with model}
\usage{
xml_bite(model, xml, out = sprintf("\%s_edited.xml", gsub(".xml", "",
  xml)))
}
\arguments{
\item{model}{an object of class "model" (see details)}

\item{xml}{name of the output file that will store the log of MCMC chain}

\item{out}{where to write the edited xml}
}
\value{
no return value, called for side effects

No return value: Modifies the .xml file given in xml in the user's filespace.
}
\description{
Modifies a .xml file from beauti to include a model model in the Beast 2 analysis
}
\details{
This function takes a .xml file generated with Beauti and a model object generated with \code{\link{make_jive}}
Only model objects that use models supported by the Beast implementation of model ("BM", c("BM", "sigma"), "WN", "OU", c("OU", "theta"), c("OU", "root"), c("OU", "root", "theta"))
}
\author{
Theo Gaboriau
}
