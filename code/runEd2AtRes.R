source("setup.R")
options(error=traceback)
index = as.numeric(commandArgs(trailingOnly = TRUE)) # test with index == something

allRes = seq(from=100, to=225, by=25)
res = allRes[index]

system.time(out <- ed2AtRes(res))
