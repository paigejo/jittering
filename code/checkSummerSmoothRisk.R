# Which smoothRisk definition does the INSTALLED SUMMER use at runtime?
# Print the lines of SUMMER's pixel-pop machinery that compute pSmoothRisk.
suppressMessages(library(SUMMER))

cat("SUMMER version:", as.character(packageVersion("SUMMER")), "\n\n")

showSmoothRiskLines <- function(fn, name) {
    txt <- deparse(fn)
    hits <- grep("smoothRisk|SmoothRisk|logitNormMean", txt, ignore.case = FALSE)
    if(length(hits) == 0) { cat("[", name, "] no smoothRisk lines\n"); return(invisible()) }
    cat("=====", name, "=====\n")
    for(h in hits) {
        lo <- max(1, h - 2); hi <- min(length(txt), h + 2)
        cat(sprintf("  L%d-%d:\n", lo, hi))
        cat(paste0("    ", txt[lo:hi]), sep = "\n")
        cat("\n")
    }
}

showSmoothRiskLines(SUMMER::simPopCustom,            "SUMMER::simPopCustom")
if(exists("simPopInternal", where = getNamespace("SUMMER")))
    showSmoothRiskLines(SUMMER:::simPopInternal,     "SUMMER:::simPopInternal")

# The pixel-level generator is usually an internal; scan all namespace fns
fns <- ls(getNamespace("SUMMER"))
for(f in fns) {
    obj <- get(f, envir = getNamespace("SUMMER"))
    if(!is.function(obj)) next
    txt <- deparse(obj)
    if(any(grepl("pSmoothRisk", txt))) {
        cat("\n### namespace fn with pSmoothRisk:", f, "\n")
        showSmoothRiskLines(obj, paste0("SUMMER:::", f))
    }
}
