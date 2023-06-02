# library(haven)
library(survey)

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
getDirectEsts = function(clustDat, divideWeight=TRUE) {
  clustIDs = clustDat$clustID
  indivDat = extendDataset(clustDat, clustIDs=clustIDs, divideWeight=divideWeight)
  
  indivDat$regionRural <- with(indivDat, interaction(area, urbanRural), drop=TRUE)
  
  res = defineSurvey(indivDat, 
                     stratVar=indivDat$regionRural,
                     useSamplingWeights = TRUE)
  
  res
}

# clustDatDHS: a subset of edVal
# clustDatMICS: a subset of edMICSval
getCombinedDirectEsts = function(clustDatDHS, clustDatMICS, divideWeight=TRUE) {
  clustDatDHS$area = NULL
  # clustDatDHS$subarea = NULL only remove this because it could be confusing
  names(clustDatDHS)[grepl("Stratum", names(clustDatDHS))] = "area"
  names(clustDatDHS)[grepl("clusterID", names(clustDatDHS))] = "clustID"
  
  names(clustDatMICS)[grepl("ys", names(clustDatDHS))] = "y"
  names(clustDatMICS)[grepl("ns", names(clustDatDHS))] = "n"
  
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
    extendDataRow(data.frame(clustDat[i,]), clustID=clustIDs[i], divideWeight=divideWeight)
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
  # admin1 = rep(clustDatRow$admin1, n)
  
  res = merge(data.frame(y, clustID, hhID, urbanRural), tmp, by="clustID")
  
  if(divideWeight)
    res$samplingWeight = res$samplingWeight / n
  
  return(res)
}

# - a function that reads in a glm or svyglm - #
# - object and returns the estimate and SE - #
# - specifics in the supplementary materials - #
## This function takes care of the delta method
## to calculate the variance of the estimates.
get.est<-function(glm.ob){
  
  beta<-summary(glm.ob)$coef[,1]
  
  est <-expit(beta)
  var.est <- vcov(glm.ob)[1,1]
  
  # compute 80% CI intervals
  lower <- logit(est)+qnorm(c(0.9))*sqrt(var.est)
  upper <- logit(est)+qnorm(c(0.1))*sqrt(var.est)
  return(c(est,lower, upper,logit(est),var.est))
}

# -- a function to subset the design based on a region and time period -- #
# -- and then run the svyglm function and return the get.est() results -- #

## First line in function allows you to subset your data and ALSO the specified
## svydesign object into area (usually v024 variable in DHS) 
## and time (per5 is a variable we construct for the 5-year periods in the Stata step)
## Second line fits the survey-weighted glm

region.time.HT<-function(dataobj, svydesign, area){
  
  tmp<-subset(svydesign, (admin1==area))
  
  tt2 <- tryCatch(glmob<-svyglm(y.x~1,
                                design=tmp,family=quasibinomial, maxit=50), 
                  error=function(e) e, warning=function(w) w)
  
  if(is(tt2, "warning")){
    if(grepl("agegroups", tt2)){
      res <- get.est(glmob)
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
    res <- get.est(glmob)
    res = c(res, 0)
    return(res)
  }
}

region.time.HTDat<-function(dataobj, svydesign, area, nationalEstimate){
  
  if(!nationalEstimate) {
    
    tmp<-subset(svydesign, (admin1==area))
    
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
      res <- get.est(glmob)
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
    res <- get.est(glmob)
    res = c(res, 0)
    return(res)
  }
}

defineSurvey <- function(dat_obj, stratVar, useSamplingWeights=TRUE, nationalEstimate=FALSE, 
                            getContrast=nationalEstimate){
  
  options(survey.lonely.psu="adjust")
  
  # --- setting up a place to store results --- #
  regions <- sort(unique(dat_obj$admin1))
  regions_num  <- 1:length(regions)
  
  if(!nationalEstimate) {
    results<-data.frame(admin1=rep(regions,each=1))
    results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
    results$converge <- NA
  }
  else {
    results<-data.frame(urban=c(TRUE, FALSE))
    results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
    results$converge <- NA
  }
  
  if(useSamplingWeights){
    dat_obj$wt <- dat_obj$samplingWeight
  } else {
    dat_obj$wt <- NULL
  }
  
  if(is.null(stratVar)){
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata =NULL,
                              weights=NULL, data=dat_obj)
  } else {
    ## not in all surveys does v022 contain the correct sampling strata
    ## Thus, the correct vector has to be provided externally
    dat_obj$strat <- stratVar
    
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata=~strat, nest=T, 
                              weights=~wt, data=dat_obj)
  }
  
  for(i in 1:nrow(results)){
    if(!nationalEstimate) {
      results[i, 2:7] <- region.time.HTDat(dataobj=dat_obj, svydesign=my.svydesign, 
                                           area=results$admin1[i], nationalEstimate=nationalEstimate)
    }
    else {
      results[i, 2:7] <- region.time.HTDat(dataobj=dat_obj, svydesign=my.svydesign, 
                                           area=i, nationalEstimate=nationalEstimate)
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
    lower = est + qnorm(0.025, sd=sqrt(urbanVar))
    upper = est + qnorm(0.975, sd=sqrt(urbanVar))
    contrastStats = list(est=est, sd=sqrt(urbanVar), lower95=lower, upper95=upper)
    return(list(results=results, contrastStats=contrastStats))
  } else {
    return(results)
  }
  
}

# naive glm not accounting for survey weights within admin1 areas
# Set dat_obj$admin1 to be something else for different kinds of aggregations
run_naive <- function(dat_obj){
  regions <- sort(unique(dat_obj$admin1))
  regions_num  <- 1:length(regions)
  
  results<-data.frame(admin1=rep(regions,each=1))
  results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
  results$converge <- NA
  
  for(i in 1:nrow(results)){
    my.glm <- glm(y.x~1, family=binomial, 
                  data=dat_obj, 
                  subset = admin1 == results$admin1[i] ) 
    # newdat = dat_obj[dat_obj$admin1==results$admin1[i], ]
    # my.glm2 <- glm(y.x~1, family=binomial, 
    #               data=newdat) 
    
    results[i, 2:7] <- c(get.est(my.glm),0)
  }
  return(results)
}

