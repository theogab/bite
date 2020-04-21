#' @title Write xml file with model
#' @description Modifies a .xml file from beauti to include a model model in the Beast 2 analysis
#' 
#' @details This function takes a .xml file generated with Beauti and a model object generated with \code{\link{make_jive}}
#' Only model objects that use models supported by the Beast implementation of model ("BM", c("BM", "sigma"), "WN", "OU", c("OU", "theta"), c("OU", "root"), c("OU", "root", "theta"))
#' 
#' @param model an object of class "model" (see details)
#' @param xml name of the output file that will store the log of MCMC chain
#' @param out where to write the edited xml
#' @import xml2
#' @export
#' @return no return value, called for side effects
#' @author Theo Gaboriau
#' @return No return value: Modifies the .xml file given in xml in the user's filespace. 
#' @encoding UTF-8


xml_bite <- function(model, xml, out = sprintf("%s_edited.xml", gsub(".xml", "", xml))){
  
  x <- read_xml(xml)
  
  treeid <- xml_attr(xml_find_first(x, "//tree"), "id")
  
  spnames <- unique(sapply(xml_find_all(x, "//taxon[@spec='Taxon']"), xml_attr, "id"))
  
  ## Likelihood
  likelihood_xml(x, model$lik, model$data$traits, spnames)
  
  ## Check
  forbiden.models <- list(c("WN", "sigma"), c("OU", "alpha"), c("OU", "sigma"))
  if(any(sapply(forbiden.models, function(m){
    all(sapply(m, grepl, model$prior.mean$name))
  }))) stop(sprintf("%s model is not supported by the Beast implementation", model$prior.mean$name))
  if(any(sapply(forbiden.models, function(m){
    all(sapply(m, grepl, model$prior.var$name))
  }))) stop(sprintf("%s model is not supported by the Beast implementation", model$prior.var$name))
  
  ## Prior
  prior_xml(x, model$prior.var, vari = "LogVar", treeid, spnames)
  prior_xml(x, model$prior.mean, vari = "Mean", treeid, spnames)
  
  ## Hyperprior
  hyperprior_xml(x, model$prior.var, vari = "LogVar")
  hyperprior_xml(x, model$prior.mean, vari = "Mean")
  
  ## operators
  operator_xml(x, model$lik, vari = "lik")
  operator_xml(x, model$prior.var, vari = "LogVar", treeid, spnames)
  operator_xml(x, model$prior.mean, vari = "Mean", treeid, spnames)
  
  pars <- xml_attr(xml_find_all(x, "//*[@spec[starts-with(., 'parameter')] and @id[starts-with(., 'model')]]"), "id")
  
  ## States
  state <- xml_find_first(x, "//state[@id = 'state']")
  for(p in pars){
    rm <- (!any(sapply(c("theta", "sigma", "alpha"), grepl, model$prior.mean$name)) & grepl("MeanAssignments", p)) |
      (!any(sapply(c("theta", "sigma", "alpha"), grepl, model$prior.var$name)) & grepl("LogVarAssignments", p) | 
         (grepl("OU", model$prior.mean$name) & !grepl("root", model$prior.mean$name) & grepl("modelMeanRootValue", p)) |
         (grepl("OU", model$prior.var$name) & !grepl("root", model$prior.var$name) & grepl("modelLogVarRootValue", p)))
    if(!rm){
      xml_add_child(state,"stateNode", idref = p, .where = 0)
    }
  } 
  
  ## Logger
  log <- xml_find_first(x, "//logger[@id = 'tracelog']")
  for(p in pars){
    rm <- (!any(sapply(c("theta", "sigma", "alpha"), grepl, model$prior.mean$name)) & grepl("MeanAssignments", p)) |
      (!any(sapply(c("theta", "sigma", "alpha"), grepl, model$prior.var$name)) & grepl("LogVarAssignments", p) | 
         (grepl("OU", model$prior.mean$name) & !grepl("root", model$prior.mean$name) & grepl("modelMeanRootValue", p)) |
         (grepl("OU", model$prior.var$name) & !grepl("root", model$prior.var$name) & grepl("modelLogVarRootValue", p)))
    if(!rm){
      xml_add_child(log,"log", idref = p, .where = 3)
    }
  } 
  
  cat(gsub("&amp;", "&", paste(xml_root(x))), append = FALSE, file = out)
  
}

