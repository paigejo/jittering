# srun --partition=CPUQ --time=02:00:00 --mem-per-cpu=20000 --pty bash

source("setup.R")
index = as.numeric(commandArgs(trailingOnly = TRUE)) # test with index == something
# index=53
# i: 1-4 corresponding to M_d, M_D, M_dm, and M_DM
# j: 1-11 for M_d and M_D and 1-20 for M_dm and M_DM
jobInds = getJobIndices(index, maxJ=c(11, 11, 20, 20), rev=TRUE)
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
system.time(out <- getValidationFit(fold, model, regenModFit=FALSE))
# })
# save(p, file="savedOutput/simStudyResults/tempFiles/profFile.RData")

# Rprof(NULL)
# profvis(prof_input = "data.Rprof")