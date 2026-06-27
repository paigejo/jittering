# One-shot prep before submitting the runValidationFE.sbatch array.
# Run this ONCE on the cluster (login node, or a single-task sbatch with
# --dependency afterok used to gate the array) so the parallel array workers
# don't race to (a) build fullMDMInputs_FE.RData or (b) compile TMB templates.

source("setup.R")
options(error=traceback)
source("code/modFED.R")
source("code/modFEM.R")
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

# 2. Pre-compile the three TMB templates the fitters use.
templates = c("modM_FEnug_GH",         # used by fitFEM
              "modD_BYM2_GH_v2",       # used by fitFED
              "modMDM_BYM2_GH_v2")     # used by fitFEMD
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

# 3. Make sure the output directories exist (cluster-level + areal).
for(d in c("savedOutput/validation/folds", "savedOutput/validation/areal")) {
    if(!dir.exists(d)) dir.create(d, recursive=TRUE)
}

# 4. Verify the predGrid inputs the AREAL jobs need are loaded by setup.R (areal
#    jobs call predGrid + predArea on the admin1 shapefile / population surface,
#    exactly as applicationReal.R does). Fail loudly here, not in 147 workers.
needed = c("popMatNGAThresh", "adm1")
missing = needed[!sapply(needed, exists)]
if(length(missing)) {
    stop("[prep] predGrid inputs missing after source('setup.R'): ",
         paste(missing, collapse=", "),
         " — areal validation jobs would fail. Ship these objects / their .RData.")
} else {
    cat("[prep] predGrid inputs present:", paste(needed, collapse=", "), "\n")
}

cat("\n[prep] all set. Submit the array with:  bash runValidationFE.sh\n")
cat("[prep] array = 1-105%16 :  100 cluster jobs + 5 areal jobs over the 5 VALMODELS.\n")
