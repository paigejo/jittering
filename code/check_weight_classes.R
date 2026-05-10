#!/usr/bin/env Rscript
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_DMSep.R")

inp <- makeInputsMDM(ed, edMICS,
                     KMICS=100,
                     KDHSurb=16, JInnerUrban=4,
                     KDHSrur=21, JInnerRural=4, JOuterRural=1,
                     admMICS=admFinal, adm2DHS=adm2Full,
                     adm2AsCovariate=FALSE)

cat("class intPtsMICS$wUrban:", class(inp$intPtsMICS$wUrban), "\n")
cat("dim intPtsMICS$wUrban:", paste(dim(inp$intPtsMICS$wUrban), collapse="x"), "\n")
cat("class intPtsMICS$wRural:", class(inp$intPtsMICS$wRural), "\n")
cat("dim intPtsMICS$wRural:", paste(dim(inp$intPtsMICS$wRural), collapse="x"), "\n")
cat("class intPtsDHS$wUrban:", class(inp$intPtsDHS$wUrban), "\n")
cat("dim intPtsDHS$wUrban:", paste(dim(inp$intPtsDHS$wUrban), collapse="x"), "\n")
cat("class intPtsDHS$wRural:", class(inp$intPtsDHS$wRural), "\n")
cat("dim intPtsDHS$wRural:", paste(dim(inp$intPtsDHS$wRural), collapse="x"), "\n")
