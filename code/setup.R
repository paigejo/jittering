# setup script for the jittering project
# install required packages if necessary
if(FALSE) {
  install.packages("Matrix")
  install.packages("spam")
  install.packages("fields")
  install.packages("LatticeKrig")
  install.packages("invgamma")
  install.packages("latex2exp")
  install.packages("xtable")
  install.packages("profvis")
  install.packages("geosphere")
  install.packages("viridis")
  install.packages("sp")
  install.packages("raster")
  install.packages("MCMCpack")
  install.packages("numDeriv")
  install.packages("INLA")
  install.packages("edfun") # for drawing from empirical distributions quickly
  install.packages("data.table")
  install.packages("sampling")
  install.packages("haven")
  install.packages("survey")
  install.packages("abind")
  install.packages("devtools")
  install.packages("ggplot2")
  install.packages("rgeos")
  install.packages("WeightedCluster")
  install.packages("colorspace")
  install.packages("TMB")
  install.packages("plyr")
  
  library(devtools)
  if(FALSE) {
    install_github("https://github.com/richardli/SUMMER/")
    install_github("https://github.com/paigejo/SUMMER/")
    install_local("~/git/SUMMER/")
    document("~/git/SUMMER/")
    load_all("~/git/SUMMER/")
  }
}

# load required packages and R scripts
library(Matrix)
library(spam)
library(fields)
library(LatticeKrig)
library(invgamma)
library(latex2exp)
library(xtable)
library(profvis)
library(geosphere)
library(viridis)
library(sp)
library(raster)
library(MCMCpack)
library(numDeriv)
library(INLA)
library(edfun) # for drawing from empirical distributions quickly
library(data.table)
library(sampling)
library(haven)
library(survey)
library(abind)
# install_github("https://github.com/richardli/SUMMER/tree/dev")
library(SUMMER)
# library(Rcpp)
library(ggplot2)
library(rgeos)
library(WeightedCluster)
library(shapefiles)
library(devtools)
library(colorspace)
library(TMB)
library(plyr)

codeDirectory <<- "~/git/jittering/code/"
figDirectory <<- "~/git/jittering/figures/"
dataDirectory <<- "~/git/jittering/data/"
outputDirectory <<- "~/git/jittering/savedOutput/"
globalDirectory <<- "~/git/jittering/savedOutput/global/"
elkDirectory <<- "~/git/LK-INLA/"

inf = sessionInfo()
if(inf$platform == "x86_64-apple-darwin17.0 (64-bit)") {
  setwd("~/git/jittering/")
  options(error=recover)
} else if(inf$platform != "x86_64-apple-darwin15.6.0 (64-bit)" && inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  INLA:::inla.dynload.workaround()
  # avoid setting too many threads and thereby using too much memory
  inla.setOption(num.threads=1)
  options(error=traceback)
  setwd("~/git/jittering/")
} else if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  setwd("~/git/jittering/")
  options(error=recover)
} else {
  setwd("~/git/jittering/")
  inla.setOption(num.threads=1) # consider raising
  options(error=recover)
}
options(warn=1)

setwd("~/git/jittering")
source("code/genericSpatialPlottingFunctions.R")
source('code/makeIntegrationPoints.R')
source('code/plotGenerator.R')
source('code/modSPDEJitter.R')
source('code/modBYM2.R')
source('code/utilityFuns.R')
source('code/covariates.R')
source('code/modAgg.R')
source('code/validation.R')
source("code/scores.R")
source("code/getDirectEsts.R")
\source("code/ed2AtRes.R")

## load in global variables made from the following script: 
if(FALSE) {
  source('code/preprocess.R')
}

out=load("savedOutput/global/NigeriaMapData.RData")
# load("savedOutput/global/datMICS.RData")
out=load("savedOutput/global/edMICS.RData")
out=load("savedOutput/global/ed.RData")
out=load("savedOutput/global/covariates.RData")
out=load("savedOutput/global/urbProps.RData")
out=load("savedOutput/global/poppaNGA.RData")
out=load("savedOutput/global/poppsubNGA.RData")
out=load("savedOutput/global/poppsubNGAThresh.RData")
out=load("savedOutput/global/popMatNGA.RData")
out=load("savedOutput/global/popMatNGAThresh.RData")
lonLimNGA = adm0@bbox[1,]
latLimNGA = adm0@bbox[2,]
adm0FullProjNGA = projNigeriaArea(adm0Full)
eastLimNGA = adm0FullProjNGA@bbox[1,]
northLimNGA = adm0FullProjNGA@bbox[2,]

# determine version of PROJ
ver = terra::gdal(lib="proj")
PROJ6 <- as.numeric(substr(ver, 1, 1)) >= 6




