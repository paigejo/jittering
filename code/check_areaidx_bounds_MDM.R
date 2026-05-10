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
load("savedOutput/global/admFinalMat.RData")
q <- prepareBYM2argumentsForTMB(admFinalMat, u=0.5, alpha=2/3,
                                constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
nAreas <- ncol(q$Q)

chk <- function(v,name){
  cat(name,": min=", min(v), " max=", max(v), " nAreas=", nAreas, "\n", sep="")
}
chk(inp$areaidxlocUrbanMICS, "areaidxlocUrbanMICS")
chk(inp$areaidxlocRuralMICS, "areaidxlocRuralMICS")
chk(inp$areaidxlocUrbanDHS, "areaidxlocUrbanDHS")
chk(inp$areaidxlocRuralDHS, "areaidxlocRuralDHS")
