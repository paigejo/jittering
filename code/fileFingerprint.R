# Print size + md5 for every RData file the scoring pipeline reads, in a
# stable order so two machines can be diff'd line-for-line. Run on both
# laptop and cluster; pipe to a text file and `diff` the two.

library(tools)

want <- c(
    # Globals from setup.R
    "savedOutput/global/NigeriaMapData.RData",
    "savedOutput/global/edMICS.RData",
    "savedOutput/global/ed.RData",
    "savedOutput/global/covariates.RData",
    "savedOutput/global/urbProps.RData",
    "savedOutput/global/poppaNGA.RData",
    "savedOutput/global/poppsubNGA.RData",
    "savedOutput/global/poppsubNGAThresh.RData",
    "savedOutput/global/popMatNGA.RData",
    "savedOutput/global/popMatNGAThresh.RData",
    "savedOutput/global/easpaNGA.RData",
    "savedOutput/global/easpaNGAMICS.RData",
    "savedOutput/global/easpaNGAed.RData",
    "savedOutput/global/easpaNGAedMICS.RData",
    "savedOutput/global/popMatNGAed.RData",
    "savedOutput/global/popMatNGAedThresh.RData",
    "savedOutput/global/poppStratMICS.RData",
    "savedOutput/global/poppStratMICSThresh.RData",
    # MICS integration points
    "savedOutput/global/intPtsMICS_100.RData",
    # Per-model populations + surveys
    "savedOutput/simStudy1/simPopsSurveys_BYM2.RData",
    "savedOutput/simStudy1/simPopsSurveys_SPDE.RData"
)

# Per-sim DHS int pts (first 10 sims of each generative model — enough to
# diagnose; full 100 fingerprint is in a follow-up if needed).
for(mdl in c("BYM2","SPDE"))
    for(i in 1:10)
        want <- c(want,
                  sprintf("savedOutput/simStudy1/intPtsDHS_simStudy1_%s_%d_K16_21.RData",
                          mdl, i))

cat(sprintf("# Fingerprint generated %s on %s\n", Sys.time(),
            Sys.info()[["nodename"]]))
cat(sprintf("# %-72s  %12s  %s\n", "path", "bytes", "md5"))
for(p in want) {
    if(file.exists(p)) {
        fi <- file.info(p)
        h  <- as.character(md5sum(p))
        cat(sprintf("  %-72s  %12d  %s\n", p, fi$size, h))
    } else {
        cat(sprintf("  %-72s  %12s  %s\n", p, "MISSING", "-"))
    }
}
