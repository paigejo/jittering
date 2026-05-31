# Compile all TMB templates used by the sim-study fitters, sequentially, so
# the parallel SPDE / BYM2 scoring runs can dyn.load existing .dlls without
# racing on `compile()`.
source("setup.R")
library(TMB)

tpls <- c(
    "modM_BYM2_GH_v2",
    "modD_BYM2_GH_v2",
    "modMDM_BYM2_GH_v2",
    "modM_DSepRepar",     # fitMd (repar)
    "modM_DSep",          # fitMd (non-repar fallback)
    "modFEM",
    "modFED",
    "modFEMD",
    "modM_FEnug"
)

for(t in tpls) {
    dllPath <- paste0("code/", t, ".dll")
    if(file.exists(dllPath)) {
        cat("[skip] ", t, " (", dllPath, " exists)\n", sep="")
        next
    }
    srcPath <- paste0("code/", t, ".cpp")
    if(!file.exists(srcPath)) {
        cat("[miss] ", t, " (", srcPath, " not found, skipping)\n", sep="")
        next
    }
    cat("[compile] ", t, " ...\n", sep="")
    t0 <- proc.time()[3]
    compile(srcPath, framework = "TMBad", safebounds = FALSE)
    cat(sprintf("[done]    %s in %.1f s\n", t, proc.time()[3] - t0))
}
cat("\nAll templates ready.\n")
