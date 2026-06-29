# Compile all TMB templates used by the sim-study AND validation fitters,
# sequentially, so the parallel SPDE / BYM2 / validation runs can dyn.load
# existing shared objects without racing on `compile()`.
#
# Run from the repo root:   Rscript code/precompileTemplates.R
#
# Tip: to force a clean rebuild first, delete the compiled artifacts:
#   rm -f code/*.o code/*.so code/*.dll      # (Linux/macOS cluster)
# then run this script.
source("code/setup.R")
library(TMB)

# The union of TMB .cpp templates the fitters actually use:
#   modM_BYM2_GH_v2 / modD_BYM2_GH_v2 / modMDM_BYM2_GH_v2 — fitMM / fitMD / fitMDM
#     (and fitFED/fitFEM/fitFEMD reuse modD_/modM_/modMDM_BYM2_GH_v2 with the
#      BYM2 field mapped out), so the FE family needs these too.
#   modM_FEnug_GH                              — FE nugget template (fitFEM/FE path)
#   modM_DSepRepar / modM_DSep                 — fitMd (observed coords; repar + fallback)
tpls <- c(
    "modM_BYM2_GH_v2",
    "modD_BYM2_GH_v2",
    "modMDM_BYM2_GH_v2",
    "modM_FEnug_GH",
    "modM_DSepRepar",
    "modM_DSep"
)

for(t in tpls) {
    srcPath <- paste0("code/", t, ".cpp")
    if(!file.exists(srcPath)) {
        cat("[miss] ", t, " (", srcPath, " not found, skipping)\n", sep="")
        next
    }
    # Cross-platform skip: an existing .o/.so/.dll means it's already built.
    artifacts <- paste0("code/", t, c(".o", ".so", ".dll"))
    if(any(file.exists(artifacts))) {
        cat("[skip] ", t, " (already compiled)\n", sep="")
        next
    }
    cat("[compile] ", t, " ...\n", sep="")
    t0 <- proc.time()[3]
    compile(srcPath, framework = "TMBad", safebounds = FALSE)
    cat(sprintf("[done]    %s in %.1f s\n", t, proc.time()[3] - t0))
}

# Sanity check the one that just changed (M_D fix).
soOrDll <- file.exists("code/modD_BYM2_GH_v2.so") || file.exists("code/modD_BYM2_GH_v2.dll")
if(!soOrDll) stop("modD_BYM2_GH_v2 did not produce a .so/.dll — check compile output above")
cat("\nAll templates ready (modD_BYM2_GH_v2 verified).\n")
