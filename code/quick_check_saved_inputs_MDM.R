#!/usr/bin/env Rscript
setwd("c:/Users/jpaige/git/jittering")
cat("Quick check of saved M_DM inputs (no TMB required)\n")

tryLoad <- function(f) {
  if(!file.exists(f)) { cat("File not found:", f, "\n"); return(NULL) }
  cat("Attempting to load:", f, "\n")
  res <- tryCatch({ out <- load(f); list(ok=TRUE, names=out) }, error=function(e) { list(ok=FALSE, err=conditionMessage(e)) })
  if(!res$ok) {
    cat("Failed to load:", f, " — ", res$err, "\n")
    return(NULL)
  }
  cat("Loaded:", paste(res$names, collapse=", "), "\n")
  return(res$names)
}

# Prefer loading component files that are likely to be compatible with base R
filesToTry <- c("savedOutput/global/edMICS.RData", "savedOutput/global/ed.RData",
                "savedOutput/global/intPtsMICS_25.RData", "savedOutput/global/intPtsDHS.RData",
                "savedOutput/global/admFinalMat.RData")

loadedNames <- list()
for(f in filesToTry) {
  nms <- tryLoad(f)
  if(!is.null(nms)) loadedNames[[f]] <- nms
}

if(length(loadedNames) == 0) stop("No compatible savedInput files could be loaded; please run makeInputsMDM or provide compatible RData files")

checkVar <- function(name) {
  if(exists(name)) {
    v <- get(name)
    if(is.null(dim(v))) {
      cat(name, ": length=", length(v), " min=", ifelse(length(v)>0, min(v), NA), " max=", ifelse(length(v)>0, max(v), NA), "\n")
    } else {
      cat(name, ": dim=", paste(dim(v), collapse="x"), "\n")
    }
  } else {
    cat(name, ": MISSING\n")
  }
}

vars <- c("areaidxlocUrbanMICS","areaidxlocRuralMICS","areaidxlocUrbanDHS","areaidxlocRuralDHS",
          "ysUrbMICS","nsUrbMICS","ysRurMICS","nsRurMICS",
          "ysUrbDHS","nsUrbDHS","ysRurDHS","nsRurDHS",
          "intPtsMICS","intPtsDHS")

for(v in vars) {
  checkVar(v)
}

# If intPtsMICS exists, inspect its internals
if(exists("intPtsMICS")) {
  ip <- get("intPtsMICS")
  cat("intPtsMICS contents:\n")
  print(names(ip))
  if(!is.null(ip$XUrb)) cat("XUrb dim:", paste(dim(as.matrix(ip$XUrb)), collapse="x"), "\n")
  if(!is.null(ip$XRur)) cat("XRur dim:", paste(dim(as.matrix(ip$XRur)), collapse="x"), "\n")
  if(!is.null(ip$wUrban)) cat("wUrban dim:", paste(dim(as.matrix(ip$wUrban)), collapse="x"), "\n")
  if(!is.null(ip$wRural)) cat("wRural dim:", paste(dim(as.matrix(ip$wRural)), collapse="x"), "\n")
}

if(exists("intPtsDHS")) {
  ipd <- get("intPtsDHS")
  cat("intPtsDHS contents:\n")
  print(names(ipd))
  if(!is.null(ipd$covsUrb)) cat("covsUrb dim:", paste(dim(as.matrix(ipd$covsUrb)), collapse="x"), "\n")
  if(!is.null(ipd$covsRur)) cat("covsRur dim:", paste(dim(as.matrix(ipd$covsRur)), collapse="x"), "\n")
  if(!is.null(ipd$wUrban)) cat("wUrban dim:", paste(dim(as.matrix(ipd$wUrban)), collapse="x"), "\n")
  if(!is.null(ipd$wRural)) cat("wRural dim:", paste(dim(as.matrix(ipd$wRural)), collapse="x"), "\n")
}

cat("Quick check completed.\n")
