#!/usr/bin/env Rscript
# Run MICS-only FE+nug model on simulated survey 1
library(TMB)
setwd("c:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_FEnug.R")

simIndex <- 1
KMICS <- 100
KDHSu <- 16
KDHSr <- 21
Qgh <- 10
covs <- c("urban","access","elev","distRiversLakes","normPop")

# load simulated surveys
load('savedOutput/simStudy1/simPopsSurveys_BYM2.RData')
datMICS <- surveysMICS[[simIndex]]
datDHS <- surveysDHS[[simIndex]]
# ensure N/Z naming
nameTab <- rbind(c("N","n"), c("N","ns"), c("Z","y"), c("Z","ys"))
for(i in 1:nrow(nameTab)){
  fr <- nameTab[i,1]; to <- nameTab[i,2]
  if(!(to %in% names(datMICS))) datMICS[[to]] <- datMICS[[fr]]
  if(!(to %in% names(datDHS))) datDHS[[to]] <- datDHS[[fr]]
}
if(!("Stratum" %in% names(datMICS))) datMICS$Stratum <- adm2ToStratumMICS(datMICS$subarea)

# load precomputed inputs
inputsFile <- sprintf("savedOutput/simStudy1/inputs_sim%d_KMICS%03d_KDHS%02d_%02d.RData", simIndex, KMICS, KDHSu, KDHSr)
if(!file.exists(inputsFile)) stop('Inputs file not found: ', inputsFile)
load(inputsFile) # loads inputsMDM

cat('Running fitMFEM (MICS-only)...\n')
res_mics <- fitMFEM(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM,
                    covariates=covs, KMICS=KMICS, Qgh=Qgh,
                    fixedEffectsOnly=TRUE, getSDs=TRUE, verbose=FALSE)

if(is.null(res_mics$TMBsd)) {
  cat('fitMFEM: sdreport failed or returned NULL.\n')
} else {
  peM <- summary(res_mics$TMBsd, 'fixed')
  cat('\nMICS fixed parameter estimates:\n')
  print(peM)
  mics_alpha <- peM['alpha','Estimate']; mics_alpha_se <- peM['alpha','Std. Error']
  mics_beta <- peM[grep('^beta', rownames(peM)),'Estimate']; mics_beta_se <- peM[grep('^beta', rownames(peM)),'Std. Error']
  cat(sprintf('\nMICS timing: opt=%.2f s sd=%.2f s\n', res_mics$totalTime, res_mics$sdTime))
  cat(sprintf('NLL=%.4f\n', res_mics$opt$objective))
}

cat('\nDone.\n')
