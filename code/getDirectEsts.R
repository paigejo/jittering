# library(haven)
library(survey)

# Gets direct estimates for all areas in the input cluster level dataset.
# clustDat: a cluster level data.frame with elements (in any order):
#    clustID
#    y
#    n
#    urban
#    samplingWeight
#    area
#    other optional elements that shouldn't be divided by n at individual level
# divideWeight: if TRUE, divides clustDat$samplingWeight by n for each individual in the cluster, 
#               if FALSE, gives individuals in the cluster the same sampling weight as the cluster
# signifs: significance levels of the CIs to generate
# stratByUrbanRural: if customStratVar isn't input, whether to stratify by 
#                    urban/rural x area or just area
# customStratVarName: name of a variable in the cluster level dataset to use for 
#                     stratification if stratification by area x urban/rural 
#                     isn't desired
# References:
# https://www.jstor.org/stable/pdf/26408229.pdf (method for svyglm)
# https://www.stat.umn.edu/geyer/5601/notes/sand.pdf (notes on sandwich estimation)
getDirectEsts = function(clustDat, divideWeight=TRUE, signifs=c(.5, .8, .9, .95), 
                         stratByUrbanRural=TRUE, customStratVarName=NULL) {
  
  clustIDs = clustDat$clustID
  indivDat = extendDataset(clustDat, clustIDs=clustIDs, divideWeight=divideWeight)
  
  # pick stratification variable. Either custom, area, or area x urban/rural
  if(!is.null(customStratVarName)) {
    # use a custom stratVar if user requests
    stratVar = indivDat[[customStratVarName]]
  } else {
    # add urban/rural stratification if necessary
    if(stratByUrbanRural) {
      indivDat$regionRural <- with(indivDat, interaction(area, urbanRural), drop=TRUE)
      stratVar = indivDat$regionRural
    } else {
      stratVar = indivDat$area
    }
  }
  
  res = defineSurvey(indivDat, 
                     stratVar=stratVar,
                     useSamplingWeights = TRUE, 
                     signifs=signifs)
  
  res
}

# helper functions -----
# function to extend full DHS dataset to binary form (i.e. from cluster level to individual level)
# clustDat: a dataset of cluster level data.frame with elements:
#    y
#    n
#    urban
#    samplingWeight
#    other optional elements that shouldn't be divided by n at individual level
# clustIDs: cluster IDs. v001 in DHS data
# divideWeight: if TRUE, divides clustDatRow$samplingWeight by n for each individual in the cluster, 
#               if FALSE, gives individuals in the cluster the same sampling weight as the cluster
extendDataset <- function(clustDat, clustIDs, divideWeight=TRUE){
  extendRow = function(r) {
    extendDataRow(data.frame(clustDat[r,]), clustID=clustIDs[r], divideWeight=divideWeight)
  }
  
  do.call("rbind", lapply(1:nrow(clustDat), extendRow))
}

# function to help extend DHS dataset to binary form (i.e. from cluster level to individual level)
# clustDatRow: a row of cluster level data.frame with elements:
#    y
#    n
#    urban
#    samplingWeight
#    other optional elements that shouldn't be divided by n at individual level
# clustID: cluster ID. v001 in DHS data
# divideWeight: if TRUE, divides clustDatRow$samplingWeight by n for each individual in the cluster, 
#               if FALSE, gives individuals in the cluster the same sampling weight as the cluster
extendDataRow <- function(clustDatRow, clustID, divideWeight=TRUE){
  
  # add extra columns for ageMonth, ageGrpD, v001, v002
  n = clustDatRow$n
  # tmp = data.frame(clustDatRow[c(1, 6:16)])
  tmp = clustDatRow
  tmp$clustID = clustID
  
  clustID = rep(clustID, n)
  
  # All 25 households are sampled
  hhID = 1:n
  
  y = c(rep(0,n-clustDatRow$y), rep(1, clustDatRow$y))
  if(clustDatRow["urban"][1,1]){
    urbanRural = rep("urban", n)
  } else {
    urbanRural = rep("rural", n)
  }
  # area = rep(clustDatRow$area, n)
  tmp$y = NULL
  tmp$n = NULL
  res = merge(data.frame(y, clustID, hhID, urbanRural), tmp, by="clustID")
  
  if(divideWeight)
    res$samplingWeight = res$samplingWeight / n
  
  return(res)
}

# - a function that reads in a glm or svyglm - #
# - object and returns the estimate and SE - #
# - specifics in the supplementary materials - #
## This function returns summary statistics about the estimate
get.est<-function(glm.ob, signifs=c(.95)) {
  
  beta<-summary(glm.ob)$coef[,1]
  
  est <-expit(beta)
  logit.var <- vcov(glm.ob)[1,1]
  
  lowerSignifs = (1-signifs)/2
  upperSignifs = 1-(1-signifs)/2
  
  # compute CI intervals
  lower <- logit(est)+qnorm(c(lowerSignifs))*sqrt(logit.var)
  upper <- logit(est)+qnorm(c(upperSignifs))*sqrt(logit.var)
  
  # calculate normal approximate on probability scale
  probEst = logitNormMeanSimple(cbind(logit(est), sqrt(logit.var)))
  # probVar = logitNormSqMeanSimple(cbind(logit(est), sqrt(logit.var))) - probEst^2
  probVar = logitNormVarSimple(cbind(logit(est), sqrt(logit.var)))
  
  out = c(probEst,probVar, logit(est),logit.var,lower, upper)
  names(out) = c("est", "var", "logit.est", "logit.var", 
                 paste("logit.lower", 100*signifs, sep=""), 
                 paste("logit.upper", 100*signifs, sep=""))
  return(out)
}

# -- a function to subset the design based on a region and time period -- #
# -- and then run the svyglm function and return the get.est() results -- #

## First line in function allows you to subset your data and ALSO the specified
## svydesign object into area (usually v024 variable in DHS) 
## and time (per5 is a variable we construct for the 5-year periods in the Stata step)
## Second line fits the survey-weighted glm

region.time.HTDat<-function(dataobj, svydesign, area, nationalEstimate, signifs=.95) {
  
  if(!nationalEstimate) {
    thisArea=area
    tmp<-subset(svydesign, (area==thisArea))
    
    tt2 <- tryCatch(glmob<-svyglm(y~1,
                                  design=tmp,family=quasibinomial, maxit=50), 
                    error=function(e) e, warning=function(w) w)
  } else {
    thisUrban = area == 1
    tmp<-subset(svydesign, (urban==thisUrban))
    tt2 <- tryCatch(glmob<-svyglm(y~1,
                                  design=tmp,family=quasibinomial, maxit=50), 
                    error=function(e) e, warning=function(w) w)
  }
  
  if(is(tt2, "warning")){
    if(grepl("agegroups", tt2)){
      res <- get.est(glmob, signifs=signifs)
      res = c(res, 2)
    } else {
      res = c(rep(NA, 5), 3)
    }
    return(res)
  }
  if(is(tt2,"error")){
    res = c(rep(NA, 5), 1)
    return(res)
  } else {
    res <- get.est(glmob, signifs=signifs)
    res = c(res, 0)
    return(res)
  }
}

defineSurvey <- function(dat_obj, stratVar, useSamplingWeights=TRUE, nationalEstimate=FALSE, 
                            getContrast=nationalEstimate, signifs=.95){
  
  options(survey.lonely.psu="adjust")
  
  # --- setting up a place to store results --- #
  regions <- sort(unique(dat_obj$area))
  regions_num  <- 1:length(regions)
  
  if(!nationalEstimate) {
    results = matrix(nrow=length(regions), ncol=6 + length(signifs)*2)
    colnames(results) = c("area", "est", "var", "logit.est", "logit.var", 
                          paste("logit.lower", 100*signifs, sep=""), 
                          paste("logit.upper", 100*signifs, sep=""), 
                          "converge")
    results = data.frame(results)
    results$area=regions
  }
  else {
    results = matrix(nrow=2, ncol=6 + length(signifs)*2)
    colnames(results) = c("urban", "est", "var", "logit.est", "logit.var", 
                          paste("lower", 100*signifs, sep=""), 
                          paste("upper", 100*signifs, sep=""), 
                          "converge")
    results = data.frame(results)
    results$urban=c(TRUE, FALSE)
  }
  
  if(useSamplingWeights){
    dat_obj$wt <- dat_obj$samplingWeight
  } else {
    dat_obj$wt <- NULL
  }
  
  if(is.null(stratVar)){
    # --- setting up the design object --- #
    ## NOTE: -the clustID denote
    ##        one stage cluster design (clustID is cluster)
    ##       -This call below specifies our survey design
    # TODO: check if weights should be NULL here
    my.svydesign <- svydesign(id= ~clustID,
                              strata =NULL,
                              weights=NULL, data=dat_obj)
  } else {
    ## not in all surveys does v022 contain the correct sampling strata
    ## Thus, the correct vector has to be provided externally
    dat_obj$strat <- stratVar
    
    # --- setting up the design object --- #
    ## NOTE: -the clustID denote
    ##        one stage cluster design (clustID is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~clustID,
                              strata=~strat, nest=T, 
                              weights=~wt, data=dat_obj)
  }
  
  for(i in 1:nrow(results)){
    if(!nationalEstimate) {
      results[i, -1] <- region.time.HTDat(dataobj=dat_obj, svydesign=my.svydesign, 
                                           area=results$area[i], nationalEstimate=nationalEstimate, 
                                           signifs=signifs)
    }
    else {
      results[i, -1] <- region.time.HTDat(dataobj=dat_obj, svydesign=my.svydesign, 
                                           area=i, nationalEstimate=nationalEstimate, 
                                           signifs=signifs)
    }
  }
  
  if(getContrast) {
    # out = svyby(~y, by = ~urban, design = svydesign, svymean)
    glmob<-svyglm(y~urban,
                  design=my.svydesign,family=quasibinomial, maxit=50)
    
    # get contrast mean and variance
    est = glmob$coefficients[2]
    urbanVar = vcov(glmob)[2,2]
    
    # get confidence interval
    lowerSignifs = (1 - signifs)/2
    upperSignifs = 1 - (1 - signifs)/2
    lower = est + qnorm(lowerSignifs, sd=sqrt(urbanVar))
    upper = est + qnorm(upperSignifs, sd=sqrt(urbanVar))
    contrastStats = list(est=est, sd=sqrt(urbanVar), lower=lower, upper=upper, signifs=signifs)
    return(list(results=results, contrastStats=contrastStats))
  } else {
    return(results)
  }
  
}


# adapted from logitnorm package.  Calculates the mean of a distribution whose 
# logit is Gaussian. Each row of muSigmaMat is a mean and standard deviation 
# on the logit scale
logitNormMeanSimple = function(muSigmaMat, logisticApproximation=FALSE, ...) {
  if(length(muSigmaMat) > 2) {
    apply(muSigmaMat, 1, logitNormMeanSimple, logisticApproximation=logisticApproximation, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    if(sigma == 0)
      expit(mu)
    else {
      if(any(is.na(c(mu, sigma))))
        NA
      else if(!logisticApproximation) {
        # numerically calculate the mean
        fExp <- function(x) exp(plogis(x, log.p=TRUE) + dnorm(x, mean = mu, sd = sigma, log=TRUE))
        integrate(fExp, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
      } else {
        # use logistic approximation
        warning("logistic approximation is not always very accurate...")
        k = 16 * sqrt(3) / (15 * pi)
        expit(mu / sqrt(1 + k^2 * sigma^2))
      }
    }
  }
}

# Calculates the second moment of a distribution whose 
# logit is Gaussian. Each row of muSigmaMat is a mean and standard deviation 
# on the logit scale of the Gaussian that is squared.
logitNormSqMeanSimple = function(muSigmaMat, ...) 
{
  if (length(muSigmaMat) > 2) {
    apply(muSigmaMat, 1, logitNormSqMeanSimple, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    if (sigma == 0) 
      SUMMER::expit(mu)^2
    else {
      if (any(is.na(c(mu, sigma)))) 
        NA
      else {
        # plogis(x) = 1/(1 + e^-x)
        # fExp(x) = exp{2 * log(1/(1 + e^-x)) + log(phi(x; mu, sigma))}
        #         = (1/(1 + e^-x))^2 * phi(x; mu, sigma)
        #         = expit(x)^2 * phi(x; mu, sigma)
        fExp <- function(x) exp(stats::plogis(x, log.p = TRUE) * 2 + 
                                  stats::dnorm(x, mean = mu, sd = sigma, log = TRUE))
        stats::integrate(fExp, mu - 10 * sigma, mu + 
                           10 * sigma, abs.tol = 0, ...)$value
      }
    }
  }
}

# Calculates the second moment of a distribution whose 
# logit is Gaussian. Each row of muSigmaMat is a mean and standard deviation 
# on the logit scale of the Gaussian that is squared.
# NOTE: could use logitNormSqMeanSimple - logitNormMeanSimple^2, but that can be 
#       negative, whereas this function cannot.
logitNormVarSimple = function(muSigmaMat, ...) 
{
  if (length(muSigmaMat) > 2) {
    apply(muSigmaMat, 1, logitNormVarSimple, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    if (sigma == 0) 
      SUMMER::expit(mu)^2
    else {
      if (any(is.na(c(mu, sigma)))) 
        NA
      else {
        # plogis(x) = 1/(1 + e^-x) = expit(x)
        # fExp(x) = (expit(x) - probMean)^2 * phi(x; mu, sigma)
        probMean = logitNormMeanSimple(muSigmaMat)
        fExp <- function(x) exp(log((plogis(x) - probMean)^2) + stats::dnorm(x, mean = mu, sd = sigma, log = TRUE))
        stats::integrate(fExp, mu - 10 * sigma, mu + 
                           10 * sigma, abs.tol = 0, ...)$value
      }
    }
  }
}
