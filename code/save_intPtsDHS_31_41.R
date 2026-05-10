#!/usr/bin/env Rscript
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
cat("Generating DHS integration points (KDHSurb=31, KDHSrur=41)...\n")

intPtsDHS <- makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban,
                                         areaNames=ed$subarea, popPrior=TRUE,
                                         numPointsUrban=31, numPointsRural=41,
                                         JInnerUrban=4, JInnerRural=4, JOuterRural=1,
                                         adminMap=adm2Full,
                                         outFile="savedOutput/global/intPtsDHS_31_41.RData",
                                         saveOutput=TRUE)

cat("Saved intPtsDHS to savedOutput/global/intPtsDHS_31_41.RData\n")
