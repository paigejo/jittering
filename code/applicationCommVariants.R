# Fit the two requested M_DM commensurate variants on the application data and
# compare to the standard (shared-intercept) M_DM:
#   - FULLY COMMENSURATE:     commIntercept=1, commSlope=1  (alpha_M~alpha, beta_M~beta)
#   - FULLY NON-COMMENSURATE: commIntercept=0, commSlope=0  (alpha_M, beta_M independent)
# Parameters reported are DHS-side (alpha, beta), matching the paper convention.
# Extracted from sdreport directly (fast; no full predGrid needed to compare params).
if(file.exists("code/setup.R")) { source("code/setup.R"); source("code/modFEMD_comm.R") } else { source("setup.R"); source("modFEMD_comm.R") }
options(error=traceback)
KMICS<-100; KDHSu<-16; KDHSr<-21; Qgh<-10
covs <- c("urban","access","elev","distRiversLakes","normPop")
outDir<-"savedOutput/application"; if(!dir.exists(outDir)) dir.create(outDir, recursive=TRUE)
rowNames <- c("(Int)","Urban","Healthcare Inaccessibility","Elevation","Dist. Rivers & Lakes","Population","sigma^2","phi","sigma_eps^2")

parTabFromSD <- function(res){
  SD<-res$TMBsd; rnd<-summary(SD,"random"); fx<-summary(SD,"fixed")
  g<-function(tb,nm) tb[rownames(tb)==nm,,drop=FALSE]
  aD<-g(rnd,"alpha"); bD<-g(rnd,"beta")
  est<-c(aD[1,1], bD[,1]); se<-c(aD[1,2], bD[,2])
  feEst<-est; feL<-est-1.96*se; feU<-est+1.96*se
  # hypers (variance scale) with delta-method-free transform of the CI endpoints
  lt<-fx["log_tau","Estimate"]; lts<-fx["log_tau","Std. Error"]
  lp<-fx["logit_phi","Estimate"]; lps<-fx["logit_phi","Std. Error"]
  le<-fx["log_tauEps","Estimate"]; les<-fx["log_tauEps","Std. Error"]
  sig2<-c(exp(-lt), exp(-(lt+1.96*lts)), exp(-(lt-1.96*lts)))
  phi <-c(plogis(lp), plogis(lp-1.96*lps), plogis(lp+1.96*lps))
  se2 <-c(exp(-le), exp(-(le+1.96*les)), exp(-(le-1.96*les)))
  tab<-rbind(cbind(feEst,feL,feU), sig2, phi, se2)
  rownames(tab)<-rowNames; colnames(tab)<-c("Est","L","U"); round(tab,3)
}
sigmaCommOf <- function(res){
  fx<-summary(res$TMBsd,"fixed")
  if(!("log_sigma_comm" %in% rownames(fx))) return("NA (mapped out)")
  l<-fx["log_sigma_comm","Estimate"]; s<-fx["log_sigma_comm","Std. Error"]
  sprintf("%.3f (%.3f, %.3f)", exp(l), exp(l-1.96*s), exp(l+1.96*s))
}

fitVariant <- function(tag, commI, commS){
  ff<-sprintf("%s/fit_%s.RData", outDir, tag)
  if(file.exists(ff)){ e<-new.env(); load(ff,envir=e); return(e$res) }
  cat(sprintf("[%s] fitting %s (commIntercept=%d, commSlope=%d) ...\n", format(Sys.time()), tag, commI, commS)); flush.console()
  t0<-proc.time()[3]
  res<-fitFEMD_comm(datDHS=ed, datMICS=edMICS, KMICS=KMICS, KDHSurb=KDHSu, KDHSrur=KDHSr,
                    Qgh=Qgh, fixedEffectsOnly=FALSE, commIntercept=commI, commSlope=commS,
                    getSDs=TRUE, verbose=FALSE)
  cat(sprintf("  fit %.1f min\n",(proc.time()[3]-t0)/60)); flush.console()
  save(res, file=ff); res
}

resFull <- fitVariant("M_DM_commFull", 1L, 1L)   # fully commensurate
resNon  <- fitVariant("M_DM_nonComm",  0L, 0L)   # fully non-commensurate
tabFull <- parTabFromSD(resFull); tabNon <- parTabFromSD(resNon)
save(tabFull, tabNon, file=sprintf("%s/application_commVariants.RData", outDir))

ci<-function(m) sprintf("(%.2f, %.2f)", m[,"L"], m[,"U"])
cmp<-data.frame(commFull_Est=tabFull[,"Est"], commFull_CI=ci(tabFull),
                nonComm_Est=tabNon[,"Est"],   nonComm_CI=ci(tabNon))
stdF<-sprintf("%s/application_M_DM.RData", outDir)
if(file.exists(stdF)){ s<-new.env(); load(stdF,envir=s); cmp$std_Est<-s$parTab[,"Est"]; cmp$std_CI<-sprintf("(%.2f, %.2f)", s$parTab[,"L"], s$parTab[,"U"]) }
rownames(cmp)<-rowNames
cat("\n===== M_DM variants vs standard M_DM (DHS-side params) =====\n"); print(cmp)
cat(sprintf("\nsigma_comm  fully-commensurate = %s\n", sigmaCommOf(resFull)))
cat(sprintf("sigma_comm  fully-non-comm     = %s\n", sigmaCommOf(resNon)))
cat(sprintf("[%s] done.\n", format(Sys.time())))
