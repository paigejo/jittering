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

cat("MICS XUrb:", paste(dim(inp$intPtsMICS$XUrb), collapse="x"), "\n")
cat("MICS XRur:", paste(dim(inp$intPtsMICS$XRur), collapse="x"), "\n")
cat("DHS covsUrb:", paste(dim(inp$intPtsDHS$covsUrb), collapse="x"), "\n")
cat("DHS covsRur:", paste(dim(inp$intPtsDHS$covsRur), collapse="x"), "\n")
cat("MICS cols:", paste(colnames(inp$intPtsMICS$XUrb), collapse=","), "\n")
cat("DHS cols:", paste(colnames(inp$intPtsDHS$covsUrb), collapse=","), "\n")
