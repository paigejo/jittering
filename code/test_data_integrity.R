# Quick data-integrity check for the Nigeria (NGA) .RData files setup.R loads.
# Focused on:
#   (1) the ordering of area / subarea / stratum identifiers,
#   (2) a handful of spot-check numerical values per object.
# Output is intentionally short and line-by-line stable so cluster vs local can be `diff`ed.
#
# Run from code/:
#   Rscript test_data_integrity.R > ../savedOutput/test/data_integrity_<host>.log 2>&1
# Then on the other machine, diff the two logs.

source("setup.R")

ngaFiles <- c(
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
    "savedOutput/global/poppStratMICSThresh.RData"
)
extraFiles <- c(
    "savedOutput/global/ed.RData",
    "savedOutput/global/edMICS.RData"
)

# Pick the column most likely to identify an area / subarea / stratum,
# in priority order (returns NULL if none of these names exist).
nameCol <- function(df) {
    candidates <- c("subarea", "area", "Stratum", "stratum", "strat",
                    "NAME_FINAL", "NAME_2", "NAME_1")
    hit <- intersect(candidates, names(df))
    if(length(hit) > 0) hit[1] else NULL
}

spotValues <- function(x, ks=c(1, 5, 25, 50, 100)) {
    ks <- ks[ks <= length(x)]
    if(length(ks) == 0) return("")
    paste(sprintf("[%d]=%.4g", ks, x[ks]), collapse="  ")
}

inspect <- function(path) {
    cat("\n----- ", path, " -----\n", sep="")
    if(!file.exists(path)) { cat("  MISSING\n"); return(invisible(NULL)) }
    e <- new.env(); nm <- load(path, envir=e)

    for(varName in nm) {
        obj <- get(varName, envir=e)
        cat("  $ ", varName, "   class=", paste(class(obj), collapse="/"), "\n", sep="")

        if(is.data.frame(obj)) {
            cat("    nrow=", nrow(obj), "   cols=", paste(names(obj), collapse=","), "\n", sep="")
            nc <- nameCol(obj)
            if(!is.null(nc)) {
                ids <- as.character(obj[[nc]])
                cat("    id col: ", nc, "\n", sep="")
                cat("    id head: ", paste(head(ids, 5), collapse=" | "), "\n", sep="")
                cat("    id tail: ", paste(tail(ids, 3), collapse=" | "), "\n", sep="")
                cat("    id at 1,5,25,50,nrow: ",
                    paste(ids[c(1, 5, 25, 50, length(ids))[c(1,5,25,50,length(ids)) <= length(ids)]],
                          collapse=" | "), "\n", sep="")
            }
            # spot-check the first numeric column
            firstNum <- names(obj)[sapply(obj, is.numeric)][1]
            if(!is.na(firstNum)) {
                v <- obj[[firstNum]]
                cat(sprintf("    num col [%s]:  sum=%.6g   %s\n",
                            firstNum, sum(v, na.rm=TRUE), spotValues(v)))
            }
        } else if(is.matrix(obj) || is.array(obj)) {
            d <- dim(obj)
            cat("    dim=", paste(d, collapse=" x "), "  sum=", sprintf("%.6g", sum(obj, na.rm=TRUE)), "\n", sep="")
            # print a handful of named-row entries if rownames exist
            rn <- rownames(obj)
            if(!is.null(rn)) {
                pick <- c(1, 5, 25, 50, nrow(obj))
                pick <- unique(pick[pick <= nrow(obj)])
                cat("    rownames at pick: ", paste(rn[pick], collapse=" | "), "\n", sep="")
            }
            cat(sprintf("    cell[1,1]=%.6g  cell[5,1]=%.6g  cell[nrow,1]=%.6g\n",
                        obj[1,1], obj[min(5,nrow(obj)),1], obj[nrow(obj),1]))
        } else if(is.numeric(obj)) {
            cat(sprintf("    length=%d  sum=%.6g   %s\n",
                        length(obj), sum(obj, na.rm=TRUE), spotValues(obj)))
        } else {
            cat("    length=", length(obj), "\n", sep="")
        }
    }
}

cat("##########################################\n")
cat("## Nigeria (NGA) data files\n")
cat("##########################################\n")
for(f in ngaFiles) inspect(f)

cat("\n##########################################\n")
cat("## Survey cluster data (ed, edMICS)\n")
cat("##########################################\n")
for(f in extraFiles) inspect(f)

# ---- area ordering checks across the in-memory objects ----
cat("\n##########################################\n")
cat("## Area ordering cross-check (in-memory after setup.R)\n")
cat("##########################################\n")
checkOrder <- function(label, value) {
    if(is.null(value)) { cat("\n>> ", label, ": missing\n", sep=""); return(invisible(NULL)) }
    x <- as.character(value)
    cat("\n>> ", label, " (n=", length(x), "):\n", sep="")
    cat("   head: ", paste(head(x, 5),  collapse=" | "), "\n", sep="")
    cat("   tail: ", paste(tail(x, 3),  collapse=" | "), "\n", sep="")
    invisible(NULL)
}
safeGet <- function(expr) tryCatch(eval(parse(text=expr)), error=function(e) NULL)

checkOrder("admFinal$NAME_FINAL",  safeGet("admFinal@data$NAME_FINAL"))
checkOrder("adm1$NAME_1",          safeGet("adm1@data$NAME_1"))
checkOrder("adm2$NAME_2",          safeGet("adm2@data$NAME_2"))
checkOrder("poppsubNGA$subarea",       safeGet("poppsubNGA$subarea"))
checkOrder("poppsubNGAThresh$subarea", safeGet("poppsubNGAThresh$subarea"))
checkOrder("easpaNGA$area",            safeGet("easpaNGA$area"))
checkOrder("easpaNGAMICS$area",        safeGet("easpaNGAMICS$area"))
checkOrder("easpaNGAed$area",          safeGet("easpaNGAed$area"))
checkOrder("easpaNGAedMICS$area",      safeGet("easpaNGAedMICS$area"))
checkOrder("poppStratMICS$strat",         safeGet("poppStratMICS$strat"))
checkOrder("poppStratMICSThresh$strat",   safeGet("poppStratMICSThresh$strat"))

cat("\nDone.\n")
