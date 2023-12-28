# srun --partition=CPUQ --time=02:00:00 --mem-per-cpu=10000 --pty bash

source("setup.R")
options(error=traceback)
index = as.numeric(commandArgs(trailingOnly = TRUE)) # test with index == something
# index=53 # i=4, j=11
# index=43
# index=41
# i: 1-4 corresponding to M_d, M_D, M_dm, and M_DM
# j: 1-11 for M_d and M_D and 1-20 for M_dm and M_DM
jobInds = getJobIndices(index, maxJ=20, rev=TRUE)
i = jobInds[1]
j = jobInds[2]

model = c("Md", "MD", "Mdm", "MDM")[i]
fold = j



# load("savedOutput/simStudyResults/spde_prevRiskSimStudyCommandArgs.RData")
# argList = spde_prevRiskSimStudyCommandArgs[[i]]
# argList$j = j

# Rprof("savedOutput/simStudyResults/tempFiles/data.Rprof", interval = 0.01, line.profiling = TRUE,
#       gc.profiling = TRUE, memory.profiling = TRUE)

# p = profvis({
# randomBeta = ifelse(model=="MDM", TRUE, FALSE)
randomBeta = TRUE
randomAlpha=TRUE
system.time(out <- getValidationFit(fold, model, regenModFit=TRUE, sep=TRUE, 
                                    randomBeta=randomBeta, randomAlpha=randomAlpha, 
                                    varClust=FALSE))

# system.time(out <- getValidationFit(fold, model, regenModFit=FALSE, regenPreds=TRUE, sep=TRUE, 
#                                     randomBeta=randomBeta, randomAlpha=randomAlpha, 
#                                     varClust=TRUE))
# })
# save(p, file="savedOutput/simStudyResults/tempFiles/profFile.RData")

# Rprof(NULL)
# profvis(prof_input = "data.Rprof")