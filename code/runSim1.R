source("setup.R")
options(error=traceback)
index = as.numeric(commandArgs(trailingOnly = TRUE)) # test with index == something

i = ((index-1) && 100) + 1
mod = ifelse(index <= 100, "M_M", "M_DM")

system.time(out <- runSimStudy1I(i, mod=mod, regenData=FALSE))
