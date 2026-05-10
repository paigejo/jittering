library(TMB)
files = list(
  "IID intpts (DHS+MICS)" = c("savedOutput/testMM_BYM2sim_IIDnugget_progress.RData", "result_opt"),
  "IID intpts (MICS only)" = c("savedOutput/testMM_BYM2sim2_IIDnugget.RData", "result_iid"),
  "BYM2fixHyp sim1" = c("savedOutput/testMM_BYM2sim_BYM2fixedHyper_sim1.RData", "result_sim1"),
  "BYM2fixHyp sim2" = c("savedOutput/testMM_BYM2sim_BYM2fixedHyper_sim2.RData", "result_sim2"),
  "trueCoords IID" = c("savedOutput/testMM_BYM2sim_trueCoords_IID_sim1.RData", "result_iid"),
  "trueCoords FE" = c("savedOutput/testMM_BYM2sim_trueCoords_FE_sim1.RData", "result_fe"),
  "trueCoords BYM2" = c("savedOutput/testMM_BYM2sim_trueCoords_BYM2_sim1.RData", "result_bym2")
)

cat(sprintf("%-25s %10s %10s %10s %10s %10s %10s\n", "Model", "alpha", "SE", "urban", "SE", "normPop", "SE"))
cat(sprintf("%-25s %10s %10s %10s %10s %10s %10s\n", "-----", "-----", "--", "-----", "--", "-------", "--"))
cat(sprintf("%-25s %10.4f %10s %10.4f %10s %10.4f %10s\n", "TRUTH", -1.25, "", 1.00, "", 0.50, ""))

for(nm in names(files)) {
  f = files[[nm]][1]; obj_nm = files[[nm]][2]
  e = new.env(); load(f, envir=e)
  v = get(obj_nm, envir=e)
  SD0 = v$TMBsd
  sr = summary(SD0, select="random")
  rn = names(SD0$par.random)
  sf = summary(SD0, select="fixed")
  fn = names(SD0$par.fixed)
  if("alpha" %in% rn) {
    aE = sr[rn=="alpha",1]; aSE = sr[rn=="alpha",2]
    bE = sr[rn=="beta",1]; bSE = sr[rn=="beta",2]
  } else {
    aE = sf[fn=="alpha",1]; aSE = sf[fn=="alpha",2]
    bE = sf[fn=="beta",1]; bSE = sf[fn=="beta",2]
  }
  cat(sprintf("%-25s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", nm, aE, aSE, bE[1], bSE[1], bE[2], bSE[2]))
}
