# Run FE-only validation locally (laptop sanity check).
# Invokes runValidationFE.R logic in-process for a chosen subset of indices.
# Default: all valid (model, fold) combinations. Set indicesToRun to a smaller subset
# (e.g. 1:4 for one fold per model) for a quick smoke test.
#
# Index decoding (uniform maxJ = 20):
#   i: 1..4 = M_D, M_M, M_DM, M_DM_comm
#   j: 1..10 = DHS held-out fold;  11..20 = MICS held-out fold

source("setup.R")
options(error=traceback)
source("code/modFEMD_comm.R")
source("code/modFED.R")
source("code/modFEMD.R")
source("code/validation.R")
source("code/validationCommon.R")

# Override commandArgs so runValidationFE.R picks up our index value.
runOne = function(index) {
    cat("\n========== local index =", index, "==========\n")
    assign("commandArgs",
           function(trailingOnly=FALSE) as.character(index),
           envir=globalenv())
    sys.source("code/runValidationFE.R", envir=new.env(parent=globalenv()))
}

indicesToRun = 1:80
# For a quick smoke test, override e.g.:
# indicesToRun = c(1, 21, 51, 71)   # one DHS fold per model

for(ind in indicesToRun) {
    job = decodeFEjob(ind, maxJ=20)
    if(!job$isValid) next
    cat("\n--- ind =", ind, " : ", job$model, " / ", job$tag, " ---\n", sep="")
    runOne(ind)
}

cat("\nAll local validation jobs complete.\n")
