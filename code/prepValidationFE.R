# One-shot prep before submitting the runValidationFE.sbatch array.
# Run this ONCE on the cluster (login node, or a single-task sbatch with
# --dependency afterok used to gate the array) so the parallel array workers
# don't race to (a) build fullMDMInputs_FE.RData or (b) compile TMB templates.

source("setup.R")
options(error=traceback)
source("code/modFEMD_comm.R")
source("code/modFED.R")
source("code/modFEMD.R")
source("code/validationCommon.R")

# 1. Pre-build the canonical full-MDM inputs (with fold columns).
fullCache = "savedOutput/validation/fullMDMInputs_FE.RData"
if(!dir.exists("savedOutput/validation")) {
    dir.create("savedOutput/validation", recursive=TRUE)
}
if(file.exists(fullCache)) {
    cat("[prep] cache already exists:", fullCache, "\n")
} else {
    cat("[prep] building", fullCache, "...\n")
    fullInp = buildFullMDMInputs(KMICS=100, KDHSurb=16, KDHSrur=21)
    save(fullInp, file=fullCache)
    cat("[prep] saved.\n")
}

# 2. Pre-compile the four TMB templates each fitter uses.
templates = c("modM_FEnug_GH",         # used by fitFEM
              "modD_BYM2_GH_v2",       # used by fitFED
              "modMDM_BYM2_GH_v2",     # used by fitFEMD
              "modMDM_BYM2_GH_comm")   # used by fitFEMD_comm
for(tpl in templates) {
    cppFile = paste0("code/", tpl, ".cpp")
    objFiles = paste0("code/", tpl, c(".o", ".so", ".dll"))
    if(any(file.exists(objFiles))) {
        cat("[prep] compiled artifact already present for", tpl, "\n")
    } else {
        cat("[prep] compiling", cppFile, "...\n")
        TMB::compile(cppFile, framework="TMBad", safebounds=FALSE)
    }
}

# 3. Make sure the per-fold output directory exists.
foldDir = "savedOutput/validation/folds"
if(!dir.exists(foldDir)) dir.create(foldDir, recursive=TRUE)

cat("\n[prep] all set. Submit the array with:  bash runValidationFE.sh\n")
