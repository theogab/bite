## xml utilities
likelihood_xml <- function(x, jive, traits, spnames){
  
  ## Find where to put Jive Likelihood
  dist <- xml_find_first(x, "//distribution[@id='likelihood' and @spec='util.CompoundDistribution']")
  spi <- sapply(spnames, grep, names(traits))
  if(any(is.na(spi))) stop("Species name do not match between Jive object and xml object")
  
  ## Likelihood node
  anc <- xml_add_child(dist, "distribution", id="JiveLikelihood", spec="contraband.JiveLikelihood", .where = 0)
  
  ## Data
  xml_add_child(anc, "sampleData", id="JiveSampleData", spec="contraband.ManyValuesOneContTrait",
                traitValues=paste(sapply(spnames, function(sp){
                  sprintf("%s=%s", sp, paste(traits[[spi[sp]]], collapse = ","))    
                }), collapse = "&#124;")
  )
  
  ## log variances
  xml_add_child(anc, "logSigmaSqs", id="JiveLogVars", spec="parameter.RealParameter", paste(log(jive$init$v.sp[spi]), collapse = " "))
  
  ## means
  xml_add_child(anc, "mus", id="JiveMeans", spec="parameter.RealParameter", paste(jive$init$m.sp[spi], collapse = " "))
  
}


prior_xml <- function(x, jive ,vari = c("Mean", "LogVar"), treeid, spnames){
  
  ## Find where to put Jive priors
  dist <- xml_find_first(x, "//distribution[@id='prior' and @spec='util.CompoundDistribution']")
  model <- jive$name
  nreg <- as.numeric(gsub(".*?\\[(.*?)\\].*", "\\1", model, perl = T))

  ## prior node
  if(grepl("WN", model)){
    anc <- xml_add_child(dist, "distribution", id = sprintf("Jive%sPrior", vari), spec = "contraband.WNLikelihood", .where = 0)
  }
  if(grepl("BM", model)){
    anc <- xml_add_child(dist, "distribution", id = sprintf("Jive%sPrior", vari), spec = "contraband.BMMVNShiftLikelihoodOneTrait", .where = 0)
  }
  if(grepl("OU", model)){
    anc <- xml_add_child(dist, "distribution", id = sprintf("Jive%sPrior", vari), spec = "contraband.OUMVNLikelihoodOneTrait", eqDist=ifelse(grepl("root", model), "false", "true"), useRootMetaData=ifelse(grepl("root", model), "false", "true"), .where = 0)
  }
  
  # tree
  if(grepl("OU", model) | grepl("BM", model)){
    xml_add_child(anc, "tree", idref = treeid)
  }
  
  # traits
  traits <- xml_add_child(anc, "oneTraitData",id=sprintf("Jive%sStash", vari), spec="contraband.OneValueContTraits", spNames=paste(spnames, collapse = ","), nTraits="1")
  xml_add_child(traits, "traitValues",idref = sprintf("Jive%ss", vari))
  
  # mapping
  if(grepl("BM", model)){
    rate <- xml_add_child(anc, ifelse(any(model %in% "sigma"), sprintf("Jive%srateManager", vari) ,"rateManager"), id = sprintf("Jive%srateManager", vari), spec="contraband.TreeToVCVMat", coalCorrection="false")
    if(grepl("sigma", model)){
      mod <- xml_add_child(rate, "branchRateModel", id=sprintf("Jive%sCatClock", vari), spec="contraband.RandomLocalColorModel", scaling="false", includeRoot="true")
      xml_add_child(mod, "indicators", id=sprintf("Jive%sShiftIndicators", vari), spec="parameter.BooleanParameter", value= paste(c(rep("true", nreg), rep("false", 2*length(spnames)-2-nreg)), collapse = " "))
      if(jive$hprior$bm.sig.1(NA)[[2]][[1]] == "Uniform"){ # Uniform hyperprior must be specified here
        xml_add_child(mod, "colors", id = sprintf("Jive%sEvolRate", vari), spec="parameter.RealParameter",
                      lower = jive$hprior$bm.sig.1(NA)[[2]][[2]][1],
                      upper = jive$hprior$bm.sig.1(NA)[[2]][[2]][2], 
                      value = jive$init$bm.sig[1])
      } else {
        xml_add_child(mod, "colors", id = sprintf("Jive%sEvolRate", vari), spec="parameter.RealParameter", lower = "0.0", value=paste(jive$init$bm.sig ,collapse = " "))
      }
      xml_add_child(mod, "tree", idref = treeid)
    } else {
      mod <- xml_add_child(rate, "branchRateModel", id=sprintf("Jive%sRateCatClock",vari), spec="contraband.RateCategoryClockModel", nCat="1")
      if(jive$hprior$bm.sig(NA)[[2]][[1]] == "Uniform"){ # Uniform hyperprior must be specified here
        xml_add_child(mod, "rates", id=sprintf("Jive%sEvolRate", vari), spec="parameter.RealParameter",
                      lower = jive$hprior$bm.sig(NA)[[2]][[2]][1],
                      upper = jive$hprior$bm.sig(NA)[[2]][[2]][2], 
                      value=paste(jive$init$bm.sig ,collapse = " "))
      } else {
        xml_add_child(mod, "rates", id=sprintf("Jive%sEvolRate", vari), spec="parameter.RealParameter", lower="0.0", jive$init$bm.sig)
      }
      xml_add_child(mod, "rateCatAssign", id=sprintf("Jive%sAssignments", vari), spec="parameter.IntegerParameter", lower="0", upper="0", paste(rep(0, 2*length(spnames)-1), collapse = " "))
      xml_add_child(mod, "tree", idref = treeid)
    }
    xml_add_child(rate, "tree", idref = treeid)
  }
  
  if(grepl("OU", model)){
    optimum <- xml_add_child(anc, "optimumManager", id="OptimumManager", spec="contraband.TreeToVCVMat", coalCorrection="false")
    if(grepl("theta", model)){
      mod <- xml_add_child(optimum, "branchRateModel", id=sprintf("Jive%sCatClock", vari), spec="contraband.RandomLocalColorModel", scaling="false", includeRoot="true")
      xml_add_child(mod, "indicators", id=sprintf("Jive%sShiftIndicators", vari), spec="parameter.BooleanParameter", value=paste(rep("true", nreg), rep("false", 2*length(spnames)-2-nreg), collapse = " "))
      if(jive$hprior$ou.the.1(NA)[[2]][[1]] == "Uniform"){ # Uniform hyperprior must be specified here
        xml_add_child(mod, "colors", id=sprintf("Jive%sTheta", vari), spec="parameter.RealParameter",
                      lower = as.character(jive$hprior$ou.the.1(NA)[[2]][[2]][1]),
                      upper = as.character(jive$hprior$ou.the.1(NA)[[2]][[2]][2]),
                      value=paste(ifelse(grepl("root", model), jive$init$ou.the[2], jive$init$ou.the[1]) ,collapse = " "))
      } else {
        xml_add_child(mod, "colors", id=sprintf("Jive%sTheta", vari), spec="parameter.RealParameter", value=paste(ifelse(grepl("root", model), jive$init$ou.the[1], jive$init$ou.the[2]) ,collapse = " "))
      }
      xml_add_child(mod, "tree", idref = treeid)
    } else {
      mod <- xml_add_child(optimum, "branchRateModel", id=sprintf("Jive%sRateCatClock", vari), spec="contraband.RateCategoryClockModel", nCat="1")
      if(jive$hprior$ou.the.1(NA)[[2]][[1]] == "Uniform"){ # Uniform hyperprior must be specified here
        xml_add_child(mod, "rates", id=sprintf("Jive%sTheta", vari), spec="parameter.RealParameter",
                      lower = as.character(jive$hprior$ou.the.1(NA)[[2]][[2]][1]),
                      upper = as.character(jive$hprior$ou.the.1(NA)[[2]][[2]][2]),
                      value=paste(ifelse(grepl("root", model), jive$init$ou.the[2], jive$init$ou.the[1]) ,collapse = " "))
      } else {
        xml_add_child(mod, "rates", id=sprintf("Jive%sTheta", vari), spec="parameter.RealParameter", ifelse(grepl("root", model), jive$init$ou.the[2], jive$init$ou.the[1]))
      }
      xml_add_child(mod, "rateCatAssign", id=sprintf("Jive%sAssignments", vari), spec="parameter.IntegerParameter", lower="0", upper="0", paste(rep(0, 2*length(spnames)-1), collapse = " "))
      xml_add_child(mod, "tree", idref=treeid)
    }
    xml_add_child(optimum, "tree", idref = treeid)
  }
  
  
  # parameters
  if(grepl("WN", model)){
    if(jive$hprior$wn.sig(NA)[[2]][[1]] == "Uniform"){
      xml_add_child(anc, "sigmaSqs", id = sprintf("Jive%sEvolRate", vari), spec="parameter.RealParameter",
                    lower = jive$hprior$wn.sig(NA)[[2]][[1]][1],
                    upper = jive$hprior$wn.sig(NA)[[2]][[1]][2],
                    jive$init$wn.sig)
    } else {
      xml_add_child(anc, "sigmaSqs", id = sprintf("Jive%sEvolRate", vari), spec="parameter.RealParameter", jive$init$wn.sig)
    }
    xml_add_child(anc, "mus", id=sprintf("Jive%sRootValue", vari), spec="parameter.RealParameter", jive$init$wn.the)
  }
  
  if(grepl("BM", model)){
    xml_add_child(anc, "rootValue", id=sprintf("Jive%sRootValue", vari), dimension="1", spec="parameter.RealParameter", jive$init$bm.the)
  }
  
  if(grepl("OU", model)){
    xml_add_child(anc, "rootValue", id=sprintf("Jive%sRootValue", vari), dimension="1", spec="parameter.RealParameter", jive$init$ou.the[1])
    xml_add_child(anc, "alpha", id=sprintf("Jive%sAlpha", vari), dimension="1", spec="parameter.RealParameter", jive$init$ou.sv)
    xml_add_child(anc, "sigmasq", id=sprintf("Jive%sEvolRate", vari), dimension="1", spec="parameter.RealParameter", jive$init$ou.sig)
  }
  
  # mapping
  if(grepl("WN", model)){
    xml_add_child(anc, "normalAssignments", id=sprintf("Jive%sAssignments", vari), spec="parameter.IntegerParameter", lower="0", upper="0", paste(rep(0, length(spnames)), collapse = " "))
  }
  
}

hyperprior_xml <- function(x, jive, vari = c("Mean", "LogVar")){
  
  ## Find where to put Jive hyperpriors
  dist <- xml_find_first(x, "//distribution[@id='prior' and @spec='util.CompoundDistribution']")
  model <- jive$name
  nreg <- grep("\\[(.*?)\\]", model)
  
  if(grepl("WN", model)){
    pars <- cbind(r = c("wn.sig", "wn.the"), beast = sprintf(c("Jive%sEvolRate", "Jive%sRootValue"), vari))
  }
  
  if(grepl("BM", model)){
    if(grepl("sigma", model)){
      pars <- cbind(r = c("bm.sig.1", "bm.the"), beast = sprintf(c("Jive%sEvolRate", "Jive%sRootValue"), vari))
    } else {
      pars <- cbind(r = c("bm.sig", "bm.the"), beast = sprintf(c("Jive%sEvolRate", "Jive%sRootValue"), vari))
    }
  }
  
  if(grepl("OU", model)){
    pars <- cbind(r = c("ou.sv.1", "ou.sig.1", "ou.the.1"), beast = sprintf(c("Jive%sAlpha", "Jive%sEvolRate", "Jive%sTheta"), vari))
    if(grepl("root", model)){
      pars <- rbind(pars, c("ou.the.0", sprintf("Jive%sRootValue", vari)))
    }
  }
  
  for(i in 1:nrow(pars)){
    hp <- jive$hprior[[pars[i,1]]](NA)[[2]]
    if(hp[[1]] == "Normal"){
      normal_hp(dist, pars[i,2], hp[[2]][1], hp[[2]][2])
    }
    if(hp[[1]] == "Lognormal"){
      lognormal_hp(dist, pars[i,2], hp[[2]][1], hp[[2]][2])
    }
    if(hp[[1]] == "Gamma"){
      gamma_hp(dist, pars[i,2], hp[[2]][1], hp[[2]][2])
    }
    if(hp[[1]] == "Exponential"){
      exp_hp(dist, pars[i,2], hp[[2]])
    }
  }
  
}

normal_hp <- function(dist, par.id, mn, sig){
  anc <- xml_add_child(dist, "distribution", id = sprintf("%sHyperPrior", par.id), spec="beast.math.distributions.Prior", x=sprintf("@%s", par.id))
  xml_add_child(anc, "distr", id= sprintf("%sNormal", par.id), spec="beast.math.distributions.Normal", offset="0.0", mean=as.character(mn), sigma=as.character(sig))
}

exp_hp <- function(dist, par.id, mn){
  anc <- xml_add_child(dist, "distribution", id = sprintf("%sHyperPrior", par.id), spec="beast.math.distributions.Prior", x=sprintf("@%s", par.id))
  xml_add_child(anc, "distr", id= sprintf("%sExponential", par.id), spec="beast.math.distributions.Exponential", offset="0.0", mean=as.character(mn))
}

lognormal_hp <- function(dist, par.id, mn, sig){
  anc <- xml_add_child(dist, "distribution", id = sprintf("%sHyperPrior", par.id), spec="beast.math.distributions.Prior", x=sprintf("@%s", par.id))
  xml_add_child(anc, "distr", id= sprintf("%sLogNormal", par.id), spec="beast.math.distributions.LogNormalDistributionNormal", offset="0.0", S=as.character(sig), M=as.character(mn))
}

gamma_hp <- function(dist, par.id, alpha, beta){
  anc <- xml_add_child(dist, "distribution", id = sprintf("%sHyperPrior", par.id), spec="beast.math.distributions.Prior", x=sprintf("@%s", par.id))
  xml_add_child(anc, "distr", id= sprintf("%sGamma", par.id), spec="beast.math.distributions.Gamma", alpha = as.character(alpha), beta = as.character(beta), mode = "ShapeScale")
}


operator_xml <- function(x, jive ,vari = c("Mean", "LogVar"), treeid, spnames){
  ## Find where to put Jive priors
  op <- xml_find_first(x, "//operator")
  
  if(vari == "lik"){
    xml_add_sibling(op, "operator", id = "JiveLogVarsRandomWalk", spec = "RealRandomWalkOperator", parameter = "@JiveLogVars", windowSize = as.character(mean(jive$ws$v.sp)), weight = "20.0", .where = "before")
    xml_add_sibling(op, "operator", id = "JiveMeansRandomWalk", spec = "RealRandomWalkOperator", parameter = "@JiveMeans", windowSize = as.character(mean(jive$ws$m.sp)), weight = "20.0", .where = "before")
  } else {
    model <- jive$name
    nreg <- grep("\\[(.*?)\\]", model)
    
    if(grepl("WN", model)){
      xml_add_sibling(op, "operator", id = sprintf("Jive%sEvolRateScaler", vari), spec = "ScaleOperator", parameter = sprintf("@Jive%sEvolRate", vari), scaleFactor = as.character(jive$ws$wn.sig), weight = "3.0", .where = "before")
      xml_add_sibling(op, "operator", id = sprintf("Jive%sRootValueRandomWalk", vari), spec = "RealRandomWalkOperator", parameter = sprintf("@Jive%sRootValue", vari), windowSize = as.character(jive$ws$wn.the), weight = "20.0", .where = "before")
    }
    
    if(grepl("BM", model)){
      if(grepl("sigma", model)){
        xml_add_sibling(op, "operator", id = sprintf("Jive%sEvolRateScaler", vari), spec = "ScaleOperator", parameter = sprintf("@Jive%sEvolRate", vari), scaleFactor = as.character(jive$ws$bm.sig[1]), weight = "3.0", .where = "before")
        shift <- xml_add_sibling(op, "operator", id = sprintf("Jive%sShiftIndicatorMove", vari), spec = "contraband.operators.BitMoveOperator",  weight = "40.0", parameter = sprintf("@Jive%sShiftIndicators", vari), k = as.character(nreg), .where = "before")
        xml_add_child(shift, "tree", idref = treeid)
      } else {
        xml_add_sibling(op, "operator", id = sprintf("Jive%sEvolRateScaler", vari), spec = "ScaleOperator", parameter = sprintf("@Jive%sEvolRate", vari), scaleFactor = as.character(jive$ws$bm.sig), weight = "3.0", .where = "before")
      }
      xml_add_sibling(op, "operator", id = sprintf("Jive%sRootValueRandomWalk", vari), spec = "RealRandomWalkOperator", parameter = sprintf("@Jive%sRootValue", vari), windowSize = as.character(jive$ws$bm.the), weight = "20.0", .where = "before")
    }
    
    if(grepl("OU", model)){
      xml_add_sibling(op, "operator", id = sprintf("Jive%sThetaRandomWalk", vari), spec = "RealRandomWalkOperator", parameter = sprintf("@Jive%sTheta", vari), windowSize = as.character(jive$ws$ou.the[ifelse(grepl("root", model), 2, 1)]), weight = "20.0", .where = "before")
      if(grepl("theta", model)){
        shift <- xml_add_sibling(op, "operator", id = sprintf("Jive%sShiftIndicatorMove", vari), spec = "contraband.operators.BitMoveOperator",  weight = "40.0", parameter = sprintf("@Jive%sShiftIndicators", vari), k = as.character(nreg), .where = "before")
        xml_add_child(shift, "tree", idref = treeid)
      }
      xml_add_sibling(op, "operator", id = sprintf("Jive%sAlphaScaler", vari), spec = "ScaleOperator", parameter = sprintf("@Jive%sAlpha", vari), scaleFactor = as.character(jive$ws$ou.sv), weight = "3.0", .where = "before")
      xml_add_sibling(op, "operator", id = sprintf("Jive%sEvolRateScaler", vari), spec = "ScaleOperator", parameter = sprintf("@Jive%sEvolRate", vari), scaleFactor = as.character(jive$ws$ou.sig), weight = "3.0", .where = "before")
      if(grepl("root", model)){
        xml_add_sibling(op, "operator", id = sprintf("Jive%sRootValueRandomWalk", vari), spec = "RealRandomWalkOperator", parameter = sprintf("@Jive%sRootValue", vari), windowSize = as.character(jive$ws$ou.the[1]), weight = "20.0", .where = "before")
      }
    }
  }
}

