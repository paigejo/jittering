# srun --partition=CPUQ --time=02:00:00 --mem-per-cpu=10000 --pty bash

source("setup.R")
options(error=traceback)
index = as.numeric(commandArgs(trailingOnly = TRUE)) # test with index == something

# i: 1-3 corresponding to M_D, M_M, and M_DM
# j: 1-100 corresponding to population/survey(s)
jobInds = getJobIndices(index, maxJ=100, rev=TRUE)
i = jobInds[1]
j = jobInds[2]

model = c("MD", "MM", "MDM")[i]
fold = j

# load the survey data for this population
out = load("savedOutput/simStudy1/simPopsSurveys.RData")
datDHS = surveysDHS[[j]]
datMICS = surveysMICS[[j]]

# load the MDM model inputs
out = load(paste0("savedOutput/simStudy1/inputsMDM_", j, ".RData"))

# Rprof("savedOutput/simStudyResults/tempFiles/data.Rprof", interval = 0.01, line.profiling = TRUE,
#       gc.profiling = TRUE, memory.profiling = TRUE)

# p = profvis({

if(model == "MD") {
  fit = fitMD(datDHS=datDHS, inputsMDM=inputsMDM)
} else if(model == "MM") {
  fit = fitMD(datMICS=datMICS, inputsMDM=inputsMDM)
} else if(model == "MDM") {
  fit = fitMD(datDHS=datDHS, datMICS=datMICS, inputsMDM=inputsMDM)
} else {
  stop(paste0("Model ", model, " not supported for simulation study 1"))
}



# })
# save(p, file="savedOutput/simStudyResults/tempFiles/profFile.RData")

# Rprof(NULL)
# profvis(prof_input = "data.Rprof")