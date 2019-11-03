#' @title Write xml file with jive
#' @description Modifies a xml file from beauti to include a jive model in the Beast 2 analysis
#' 
#' @details This function takes a xml file generated with Beauti and a jive object generated with \code{\link{make_jive}}
#' Only jive objects that use models supported by the Beast implementation of jive ("BM", c("BM", "sigma"), "WN", "OU", c("OU", "theta"), c("OU", "root"), c("OU", "root", "theta"))
#' 
#' @param jive an object of class "jive" (see details)
#' @param xml name of the output file that will store the log of MCMC chain
#' @param out where to write the edited xml
#' @import xml2
#' @export
#' @author Theo Gaboriau
#' @return none


xml_jive <- function(jive, xml, out = sprintf("%s_edited.xml", gsub(".xml", "", xml))){
  
  x <- read_xml(xml)
  
  treeid <- xml_attr(xml_find_first(x, "//tree"), "id")
  
  spnames <- unique(sapply(xml_find_all(x, "//taxon[@spec='Taxon']"), xml_attr, "id"))
  
  ## Likelihood
  likelihood_xml(x, jive$lik, jive$data$traits, spnames)
  
  ## Check
  forbiden.models <- list(c("WN", "sigma"), c("OU", "alpha"), c("OU", "sigma"))
  if(any(sapply(forbiden.models, function(m){
    all(sapply(m, grepl, jive$prior.mean$name))
  }))) stop(sprintf("%s model is not supported by the Beast implementation", jive$prior.mean$name))
  if(any(sapply(forbiden.models, function(m){
    all(sapply(m, grepl, jive$prior.var$name))
  }))) stop(sprintf("%s model is not supported by the Beast implementation", jive$prior.var$name))
  
  ## Prior
  prior_xml(x, jive$prior.var, vari = "LogVar", treeid, spnames)
  prior_xml(x, jive$prior.mean, vari = "Mean", treeid, spnames)
  
  ## Hyperprior
  hyperprior_xml(x, jive$prior.var, vari = "LogVar")
  hyperprior_xml(x, jive$prior.mean, vari = "Mean")
  
  ## operators
  operator_xml(x, jive$lik, vari = "lik")
  operator_xml(x, jive$prior.var, vari = "LogVar", treeid, spnames)
  operator_xml(x, jive$prior.mean, vari = "Mean", treeid, spnames)
  
  pars <- xml_attr(xml_find_all(x, "//*[@spec[starts-with(., 'parameter')] and @id[starts-with(., 'Jive')]]"), "id")
  
  ## States
  state <- xml_find_first(x, "//state[@id = 'state']")
  for(p in pars){
    rm <- (!any(sapply(c("theta", "sigma", "alpha"), grepl, jive$prior.mean$name)) & grepl("MeanAssignments", p)) |
      (!any(sapply(c("theta", "sigma", "alpha"), grepl, jive$prior.var$name)) & grepl("LogVarAssignments", p) | 
         (grepl("OU", jive$prior.mean$name) & !grepl("root", jive$prior.mean$name) & grepl("JiveMeanRootValue", p)) |
         (grepl("OU", jive$prior.var$name) & !grepl("root", jive$prior.var$name) & grepl("JiveLogVarRootValue", p)))
    if(!rm){
      xml_add_child(state,"stateNode", idref = p, .where = 0)
    }
  } 
  
  ## Logger
  log <- xml_find_first(x, "//logger[@id = 'tracelog']")
  for(p in pars){
    rm <- (!any(sapply(c("theta", "sigma", "alpha"), grepl, jive$prior.mean$name)) & grepl("MeanAssignments", p)) |
      (!any(sapply(c("theta", "sigma", "alpha"), grepl, jive$prior.var$name)) & grepl("LogVarAssignments", p) | 
         (grepl("OU", jive$prior.mean$name) & !grepl("root", jive$prior.mean$name) & grepl("JiveMeanRootValue", p)) |
         (grepl("OU", jive$prior.var$name) & !grepl("root", jive$prior.var$name) & grepl("JiveLogVarRootValue", p)))
    if(!rm){
      xml_add_child(log,"log", idref = p, .where = 3)
    }
  } 
  
  write_xml(x, out)
  
}

