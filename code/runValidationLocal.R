# run validation locally

# source("setup.R")
# options(error=recover)
# index = as.numeric(commandArgs(trailingOnly = TRUE)) # test with index == something
# index=53 # i=4, j=11
# index=43
# index=41
# i: 1-4 corresponding to M_d, M_D, M_dm, and M_DM
# j: 1-11 for M_d and M_D and 1-20 for M_dm and M_DM

for(ind in 1:80) {
  index = ind
  
  print(paste0("ind: ", index))
  
  jobInds = getJobIndices(index, maxJ=20, rev=TRUE)
  i = jobInds[1]
  j = jobInds[2]
  
  model = c("Md", "MD", "Mdm", "MDM")[i]
  fold = j
  
  randomBeta = TRUE
  system.time(out <- getValidationFit(fold, model, regenModFit=TRUE, randomBeta=randomBeta))
}
