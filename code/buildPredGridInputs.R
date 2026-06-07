# Precompute the inputs predGrid() pulls from rasters at every fit. They're
# pure functions of popMatNGAThresh's pixel coordinates, so we evaluate them
# ONCE on the laptop (which has a working raster/GDAL stack) and ship the
# resulting RData with the rest of savedOutput/global. The cluster then
# never has to call terra::extract, eliminating:
#   - per-fit raster I/O (small but adds up across 800+ fits)
#   - GDAL/proj cluster-vs-laptop version sensitivity
#   - covariates.RData cross-machine drift
#
# Saves: savedOutput/global/predGridInputs.RData with
#   predGridXmat       — getDesignMat(LLcoords, normalized=TRUE) result
#   predGridXmatRaw    — getDesignMat(LLcoords, normalized=FALSE) result (for
#                        callers that pass normalized=FALSE — rare)
#   predGridPopValsN   — (log1p(pop) - popMean) / popSD  (normalized log-pop)
#   predGridLLcoords   — the cbind(lon, lat) we evaluated on (for sanity check)
#   predGridPopMatNrow — nrow(popMatNGAThresh) (for sanity check on load)

source("setup.R")

cat(sprintf("\n=== %s | precomputing predGrid inputs ===\n", format(Sys.time())))

LLcoords <- cbind(popMatNGAThresh$lon, popMatNGAThresh$lat)
cat(sprintf("nrow(popMatNGAThresh) = %d\n", nrow(popMatNGAThresh)))

t0 <- proc.time()[3]
predGridXmat    <- getDesignMat(LLcoords, normalized = TRUE)
cat(sprintf("getDesignMat(normalized=TRUE)  : %.1f sec\n", proc.time()[3] - t0))

t1 <- proc.time()[3]
predGridXmatRaw <- getDesignMat(LLcoords, normalized = FALSE)
cat(sprintf("getDesignMat(normalized=FALSE) : %.1f sec\n", proc.time()[3] - t1))

# Match the popValsNorm computation inside predGrid (line 1028-1031 of modBYM2.R)
out <- load("savedOutput/global/popMeanSDCal.RData")
predGridPopValsN <- (log1p(predGridXmat[, 2]) - popMeanCalThresh) * (1 / popSDCalThresh)

predGridLLcoords   <- LLcoords
predGridPopMatNrow <- nrow(popMatNGAThresh)

outFile <- "savedOutput/global/predGridInputs.RData"
save(predGridXmat, predGridXmatRaw, predGridPopValsN,
     predGridLLcoords, predGridPopMatNrow,
     file = outFile)
cat(sprintf("Saved %s  (%.1f MB)\n", outFile,
            file.info(outFile)$size / 1024 / 1024))
cat(sprintf("\n=== %s | done ===\n", format(Sys.time())))
