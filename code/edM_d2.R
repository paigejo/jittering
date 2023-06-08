# script for women's secondary education in Nigeria application

# load datasets ----
out = load("savedOutput/global/ed.RData")
out = load("savedOutput/global/edMICS.RData")

# set parameters ----
# Umut settings: 
#   5 urban rings of 15 each (61 points total)
#   10 rural rings of 15 each (136 points total)
# KMICS=100
# KDHSurb = 31 # 4 rings of 10 each
# JInnerUrban = 4
# KDHSrur = 71 # 4 inner + 4 outer rings of 10 each
# JInnerRural = 4
# JOuterRural = 4
KMICS=25
KDHSurb = 11 # 3 rings of 5 each
JInnerUrban = 3
KDHSrur = 16 # 3 inner + 1 outer rings of 5 each
JInnerRural = 3
JOuterRural = 1

if(FALSE) {
  # do some precomputation ----
  
  # make integration points if necessary
  intPtsMICS = makeAllIntegrationPointsMICS(kmresFineStart=2.5, loadSavedIntPoints=FALSE, 
                                            numPtsRur=KMICS, numPtsUrb=KMICS)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE)
  # intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, popPrior=TRUE, 
  #                                         numPointsUrban=KDHSurb, numPointsRural=KDHSrur, 
  #                                         JInnerUrban=JInnerUrban, JInnerRural=JInnerRural, 
  #                                         JOuterRural=JOuterRural)
  intPtsDHS = makeAllIntegrationPointsDHS(cbind(ed$east, ed$north), ed$urban, 
                                          areaNames=ed$subarea, popPrior=TRUE, 
                                          numPointsUrban=KDHSurb, numPointsRural=KDHSrur, 
                                          JInnerUrban=JInnerUrban, JInnerRural=JInnerRural, 
                                          JOuterRural=JOuterRural, adminMap=adm2Full)
  
  load("savedOutput/global/intPtsDHS.RData")
  load("savedOutput/global/intPtsMICS.RData")
  
  # AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  # ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  AUrbDHS = makeApointToArea(rep(ed$subarea[ed$urban], times=KDHSurb), adm2$NAME_2) # 775 x 6259 nArea x nObsUrb
  ARurDHS = makeApointToArea(rep(ed$subarea[!ed$urban], times=KDHSrur), adm2$NAME_2) # 775 x 12960
  
  # modify the integration points to be in the correct format for TMB
  allNumPerStrat = aggregate(edMICS$Stratum, by=list(strat=edMICS$Stratum, urb=edMICS$urban), FUN=length, drop=FALSE)
  numPerStratUrb = allNumPerStrat[allNumPerStrat[,2], 3]
  numPerStratRur = allNumPerStrat[!allNumPerStrat[,2], 3]
  numPerStratRur[is.na(numPerStratRur)] = 0
  
  # first extract only the relevant covariates
  # XUrb = intPtsMICS$XUrb # XUrb is 1025 x 16 [K x nStrat] x nVar
  # AUrbMICS = makeApointToArea(edMICS$Stratum[edMICS$urban], admFinal$NAME_FINAL)
  # TODO: EXTEND AMICS TO BE LARGER, INCLUDE DIFFERENT ROW FOR EACH INTEGRATION POINT AND OBSERVATION
  
  # numPerStratUrb = table(edMICS$Stratum[edMICS$urban])
  # stratIndexUrb = unlist(mapply(rep, 1:nrow(AUrbMICS), each=numPerStratUrb * KMICS))
  # obsIndexUrb = rep(1:sum(numPerStratUrb), KMICS)
  # intPtIndexUrb = rep(1:sum(numPerStratUrb), each=KMICS)
  # actualIndexUrb = unlist(mapply(rep, 1:nrow(XUrb), each=rep(numPerStratUrb, times=KMICS)))
  # XUrb = XUrb[actualIndexUrb,] # now XUrb is [K * nObsUrb] x nVar
  # AUrbMICS = makeApointToArea(XUrb$subarea, adm2$NAME_2)
  # XUrb = XUrb[,names(XUrb) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  # XRur = intPtsMICS$XRur # XRur is 1025 x 16 [nStrat * K] x nVar
  # ARurMICS = makeApointToArea(edMICS$Stratum[!edMICS$urban], admFinal$NAME_FINAL)
  # numPerStratRur = table(edMICS$Stratum[!edMICS$urban])
  # stratIndexRur = unlist(mapply(rep, 1:nrow(ARurMICS), each=numPerStratRur * KMICS))
  # obsIndexRur = rep(1:sum(numPerStratRur), KMICS)
  # intPtIndexRur = rep(1:sum(numPerStratRur), each=KMICS)
  # actualIndexRur = unlist(mapply(rep, 1:nrow(XRur), each=rep(numPerStratRur, times=KMICS)))
  # XRur = XRur[actualIndexRur,] # now XRur is [K * nObsRur] x nVar
  # ARurMICS = makeApointToArea(XRur$subarea, adm2$NAME_2)
  # XRur = XRur[,names(XRur) %in% c("strat", "int", "urban", "access", "elev", "distRiversLakes", "normPop")]
  
  # w matrices are nStrata x K. They should be nObs x K
  # wUrban = intPtsMICS$wUrban
  # stratIndexUrbW = unlist(mapply(rep, 1:nrow(wUrban), each=numPerStratUrb))
  # wUrban = wUrban[stratIndexUrbW,]
  # 
  # wRural = intPtsMICS$wRural
  # stratIndexRurW = unlist(mapply(rep, 1:nrow(wRural), each=numPerStratRur))
  # wRural = wRural[stratIndexRurW,]
  
  # make sure the dataset aligns with this ordering, i.e. is sorted by stratum and urbanicity
  # stratIDs = match(edMICS$Stratum, admFinal$NAME_FINAL)
  # edMICS = edMICS[order(stratIDs),]
  
  # extract cluster information (in the correct order)
  # ysUrbMICS = edMICS[edMICS$urban,]$ys
  # nsUrbMICS = edMICS[edMICS$urban,]$ns
  # ysRurMICS = edMICS[!edMICS$urban,]$ys
  # nsRurMICS = edMICS[!edMICS$urban,]$ns
  
  ysUrbDHS = ed$y[ed$urban]
  ysRurDHS = ed$y[!ed$urban]
  nsUrbDHS = ed$n[ed$urban]
  nsRurDHS = ed$n[!ed$urban]
  
  # make sure A matrices are nArea x nObs, as TMB expects
  # AUrbMICS = t(AUrbMICS)
  # ARurMICS = t(ARurMICS)
  AUrbDHS = t(AUrbDHS)
  ARurDHS = t(ARurDHS)
  # mode(AUrbMICS) = "numeric"
  # mode(ARurMICS) = "numeric"
  mode(AUrbDHS) = "numeric"
  mode(ARurDHS) = "numeric"
  
  # save everything
  # intPtsMICS$XUrb = XUrb[,-(2:3)] # don't include strata or intercept
  # intPtsMICS$XRur = XRur[,-(2:3)]
  # intPtsMICS$XUrb = as.matrix(intPtsMICS$XUrb)
  # intPtsMICS$XRur = as.matrix(intPtsMICS$XRur)
  # intPtsMICS$wUrban = wUrban
  # intPtsMICS$wRural = wRural
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[,-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[,-1]
  
  # convert A matrices to sparse matrices
  # AUrbMICS = as(AUrbMICS, "sparseMatrix")
  # ARurMICS = as(ARurMICS, "sparseMatrix")
  AUrbDHS = as(AUrbDHS, "sparseMatrix")
  ARurDHS = as(ARurDHS, "sparseMatrix")
  # AUrbMICS = as.matrix(AUrbMICS)
  # ARurMICS = as.matrix(ARurMICS)
  # AUrbDHS = as.matrix(AUrbDHS)
  # ARurDHS = as.matrix(ARurDHS)
  
  save(AUrbDHS, ARurDHS, intPtsDHS, 
       ysUrbDHS, ysRurDHS, nsUrbDHS, nsRurDHS, 
       file="savedOutput/global/edM_D2Inputs.RData")
  
  # compile model ----
  # dyn.unload( dynlib("code/modBYM2JitterDHS2sparse"))
  # compile( "code/modBYM2JitterDHS2sparse.cpp", framework="TMBad", safebounds=FALSE)
  
  dyn.unload( dynlib("code/modBYM2JitterDHS2"))
  compile( "code/modBYM2JitterDHS2.cpp", 
           framework="TMBad", safebounds=FALSE)
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/TMB/include" -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/RcppEigen/include"  -DTMB_SAFEBOUNDS -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modBYM2JitterDHS2  -DTMB_LIB_INIT=R_init_modBYM2JitterDHS2  -DTMBAD_FRAMEWORK  -I/usr/local/include   -fPIC  -Wall -g -O2  -c code/modBYM2JitterDHS2.cpp -o code/modBYM2JitterDHS2.o
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o code/modBYM2JitterDHS2.so code/modBYM2JitterDHS2.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
  
  # on Idun:
  # g++ -std=gnu++14 -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/include" -DNDEBUG -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/TMB/include" -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/RcppEigen/include"   -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modBYM2JitterDHS2  -DTMB_LIB_INIT=R_init_modBYM2JitterDHS2  -DTMBAD_FRAMEWORK  -I/cluster/apps/eb/software/OpenSSL/1.1/include -I/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/include -I/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/include -I/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/include -I/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/include -I/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/include -I/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Java/11.0.2/include -I/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/include -I/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/include -I/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include/flexiblas   -fpic  -O2 -ftree-vectorize -march=native -fno-math-errno  -c code/modBYM2JitterDHS2.cpp -o code/modBYM2JitterDHS2.o
  # g++ -std=gnu++14 -shared -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -L/cluster/apps/eb/software/OpenSSL/1.1/lib64 -L/cluster/apps/eb/software/OpenSSL/1.1/lib -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib64 -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib64 -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Java/11.0.2/lib64 -L/cluster/apps/eb/software/Java/11.0.2/lib -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib64 -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib -L/cluster/apps/eb/software/GCCcore/11.3.0/lib64 -L/cluster/apps/eb/software/GCCcore/11.3.0/lib -o code/modBYM2JitterDHS2.so code/modBYM2JitterDHS2.o -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -lR
}

# load in TMB function inputs
out = load("savedOutput/global/edM_D2Inputs.RData")

# set priors ----
alpha_pri = c(0, 100^2)
beta_pri = c(0, 10^2)

out = load("savedOutput/global/adm2Mat.RData")
bym2ArgsTMB = prepareBYM2argumentsForTMB(adm2Mat, u=0.5, alpha=2/3, 
                                         constr=TRUE, scale.model=TRUE, matrixType="TsparseMatrix")
lambdaTau = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for bym2 precision
lambdaTauEps = getLambdaPCprec(u=1, alpha=.1) # get PC prior lambda for nugget precision

# Specify inputs for TMB ----

## specify random effects
rand_effs <- c('Epsilon_bym2', 'nuggetUrbDHS', 'nuggetRurDHS')

# collect input data

data_full = list(
  y_iUrbanDHS=ysUrbDHS, # same as above but for DHS survey
  y_iRuralDHS=ysRurDHS, # 
  n_iUrbanDHS=nsUrbDHS, # number binomial trials
  n_iRuralDHS=nsRurDHS, # 
  AprojUrbanDHS=AUrbDHS, # [nIntegrationPointsUrban * nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  AprojRuralDHS=ARurDHS, # 
  X_betaUrbanDHS=intPtsDHS$covsUrb, # [nIntegrationPointsUrban * nObsUrban] x nPar design matrix. Indexed mod numObsUrban
  X_betaRuralDHS=intPtsDHS$covsRur, # 
  wUrbanDHS=intPtsDHS$wUrban, # nObsUrban x nIntegrationPointsUrban weight matrix
  wRuralDHS=intPtsDHS$wRural, # 
  
  Q_bym2=bym2ArgsTMB$Q, # BYM2 unit scaled structure matrix
  V_bym2=bym2ArgsTMB$V, # eigenvectors of Q (i.e. Q = V Lambda V^T)
  alpha_pri=alpha_pri, # 2-vector with (Gaussian) prior mean and variance for intercept
  beta_pri=beta_pri, # 2-vector with (Gaussian) prior mean and variance for covariates
  tr=bym2ArgsTMB$tr, # precomputed for Q_bym2
  gammaTildesm1=bym2ArgsTMB$gammaTildesm1, # precomputed for Q_bym2
  lambdaPhi=bym2ArgsTMB$lambda, # precomputed for Q_bym2
  lambdaTau=lambdaTau, # determines PC prior for tau
  lambdaTauEps=lambdaTauEps, 
  options=0 # 1 for adreport of log tau and logit phi
)

anyna = function(x) {any(is.na(x))}
myDim = function(x) {
  dims = dim(x)
  if(is.null(dims)) {
    length(x)
  }
  else {
    dims
  }
}

if(FALSE) {
  # for testing purposes
  sapply(data_full, anyna)
  sapply(data_full, myDim)
  hist(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS, breaks=50)
  hist(data_full$y_iRuralDHS/data_full$n_iRuralDHS, breaks=50)
  hist(data_full$n_iUrbanDHS-data_full$y_iUrbanDHS, breaks=seq(-0.5, 17.5, by=1))
  hist(data_full$n_iRuralDHS-data_full$y_iRuralDHS, breaks=seq(-0.5, 24.5, by=1))
  hist(data_full$y_iRuralDHS, breaks=seq(-0.5, 16.5, by=1))
  mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 1)
  mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 0)
  mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 1) + mean(data_full$y_iUrbanDHS/data_full$n_iUrbanDHS == 0)
  mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 0) + mean(data_full$y_iRuralDHS/data_full$n_iRuralDHS == 1)
  
  hist(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS, breaks=50)
  hist(data_full$y_iRuralMICS/data_full$n_iRuralMICS, breaks=50)
  hist(data_full$n_iUrbanMICS-data_full$y_iUrbanMICS, breaks=seq(-0.5, 17.5, by=1))
  hist(data_full$y_iRuralMICS, breaks=seq(-0.5, 16.5, by=1))
  mean(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS == 1)
  mean(data_full$y_iRuralMICS/data_full$n_iRuralMICS == 0)
  mean(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS == 1) + mean(data_full$y_iUrbanMICS/data_full$n_iUrbanMICS == 0)
  mean(data_full$y_iRuralMICS/data_full$n_iRuralMICS == 0) + mean(data_full$y_iRuralMICS/data_full$n_iRuralMICS == 1)
}

# initial parameters
initUrbP = sum(c(data_full$y_iUrbanMICS, data_full$y_iUrbanDHS))/sum(c(data_full$n_iUrbanMICS, data_full$n_iUrbanDHS))
initRurP = sum(c(data_full$y_iRuralMICS, data_full$y_iRuralDHS))/sum(c(data_full$n_iRuralMICS, data_full$n_iRuralDHS))
initAlpha = logit(initRurP)
initBeta1 = logit(initUrbP) - initAlpha

tmb_params <- list(alpha = initAlpha, # intercept
                   beta = c(initBeta1, rep(0, ncol(intPtsDHS$covsUrb)-1)), 
                   log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   logit_phi = 0, # SPDE parameter related to the range
                   log_tauEps = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                   nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                   nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
)

# make TMB fun and grad ----
# dyn.load( dynlib("code/modBYM2JitterDHS2sparse"))
dyn.load( dynlib("code/modBYM2JitterDHS2"))
TMB::config(tmbad.sparse_hessian_compress = 1)
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modBYM2JitterDHS2')
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modBYM2JitterDHS2')
# objFull <- MakeADFun(data=data_full,
#                      parameters=tmb_params,
#                      hessian=TRUE,
#                      DLL='modBYM2JitterDHS2')

lower = rep(-10, length(obj[['par']]))
upper = rep( 10, length(obj[['par']]))

# make wrapper functions that print out parameters and function values
funWrapper = function(par, badParVal=1e10) {
  if(any(par < lower) || any(par > upper)) {
    return(badParVal)
  }
  print(par)
  objVal = testObj[['fn']](par)
  parNames = names(par)
  parVals = par
  parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
  parStr = paste(parStrs, collapse=", ")
  print(paste0("objective: ", objVal, " for parameters, ", parStr))
  objVal
}

grWrapper = function(par, badParVal=1e10) {
  if(any(par < lower) || any(par > upper)) {
    return(rep(badParVal, length(par)))
  }
  print(par)
  grVal = testObj[['gr']](par)
  parNames = names(par)
  parVals = par
  parStrs = sapply(1:length(par), function(ind) {paste(parNames[ind], ": ", parVals[ind], sep="")})
  parStr = paste(parStrs, collapse=", ")
  grStr = paste(grVal, collapse=", ")
  print(paste0("gradient: ", grStr, " for parameters, ", parStr))
  grVal
}

if(FALSE) {
  # testing before running optimization
  initAlphaBeta = rep(1, 6)
  initParFull = objFull$par
  initParFull[1:6] = initAlphaBeta
  objFull$fn(initParFull)
  testRep = obj$report(initParFull)
  
  system.time(test <- obj$fn(obj$par))
  # 363.236  for non-sparse
  # 355.211 for sparse
  
  system.time(test <- obj$gr(obj$par))
  # 54.5 for non-sparse
  # 55.3 for sparse
  
  testRep = obj$report()
  
  range(testRep$latentFieldUrbMICS)
  range(testRep$latentFieldRurMICS)
  range(testRep$latentFieldUrbDHS)
  range(testRep$latentFieldRurDHS)
  
  range(testRep$fe_iUrbanMICS)
  range(testRep$fe_iRuralMICS)
  range(testRep$fe_iUrbanDHS)
  range(testRep$fe_iRuralDHS)
  
  range(testRep$projepsilon_iUrbanMICS)
  range(testRep$projepsilon_iRuralMICS)
  range(testRep$projepsilon_iUrbanDHS)
  range(testRep$projepsilon_iRuralDHS)
  
  range(testRep$liksUrbMICS)
  range(testRep$liksRurMICS)
  range(testRep$liksUrbDHS)
  range(testRep$liksRurDHS)
  
  range(testRep$logDet)
  sapply(testRep, range)
  
  any(apply(data_full$wUrbanDHS, 1, function(x) {all(x == 0)}))
  any(apply(data_full$wRuralDHS, 1, function(x) {all(x == 0)}))
  any(apply(data_full$wUrbanMICS, 1, function(x) {all(x == 0)}))
  any(apply(data_full$wRuralMICS, 1, function(x) {all(x == 0)}))
}

# * Run TMB ----

{
  # tolSeq = c(1e-06, 1e-08, 1e-10, 1e-12, 1e-14)
  tolSeq = 1e-06
  testObj = obj
  optPar = testObj$par
  startTime = proc.time()[3]
  for(thisTol in tolSeq) {
    testObj = obj
    testObj$env$inner.control = list(maxit=1000, tol10=thisTol)
    testObj$env$tracepar = TRUE
    print(paste0("optimizing for tol = ", thisTol, "."))
    opt1 <- optim(par=optPar, fn=funWrapper, gr=grWrapper,
                  method = c("BFGS"), hessian = FALSE, control=list(reltol=thisTol))
    optPar = opt1$par
    if(!is.null(opt1$message)) {
      print(paste0("error for tol = ", thisTol, ". Message:"))
      print(opt1$message)
      next
    }
    else {
      print(paste0("completed optimization for tol = ", thisTol, ""))
      
      ## Get standard errors
      print("getting standard errors...")
      sdTime = system.time(
        SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                             bias.correct = TRUE,
                             bias.correct.control = list(sd = TRUE))
      )[3]
      # SD0
      print(sdTime/60)
      # 3.9447 minutes for intern=FALSE
      
      if(SD0$pdHess) {
        print("Optimization and PD hess calculation done!")
        break
      }
      else {
        print("Hessan not PD. Rerunning optimization with stricter tol...")
        
      }
    }
    
    
    
  }
  endTime = proc.time()[3]
  sdTime/60
  totalTime = endTime - startTime
  print(paste0("optimization took ", totalTime/60, " minutes"))
  # optimization took 174.831866666667 minutes (for intern=FALSE)
}

if(FALSE) {
  badPar = c(SD0$par.fixed, SD0$par.random)
  
  tempGrad = obj$env$f(badPar,order=1)
  tempHess = obj$env$spHess(badPar,random=TRUE)
  range(tempHess)
  
  tempEig = eigen(tempHess)
  range(tempEig$values)
}

if(!SD0$pdHess) {
  thisTMBpar = tmb_params
  thisTMBpar$log_tauEps = optPar[names(optPar) == "log_tauEps"]
  
  map = as.list(factor(c(alpha=1, beta=2:6, log_tau=7, logit_phi=8, log_tauEps=NA)))
  temp = unlist(map[grepl("beta", names(map))])
  map[grepl("beta", names(map))] = NULL
  map$beta = temp
  map=list(factor(c(log_tauEps=NA)))
  names(map) = "log_tauEps"
  objFixed <- MakeADFun(data=data_full,
                        parameters=tmb_params,
                        random=rand_effs,
                        map=map, 
                        hessian=TRUE,
                        DLL='modBYM2JitterFusionNugget')
  testObj = objFixed
  thisOptPar = optPar[-which(names(optPar) == "log_tauEps")]
  lower = lower[-which(names(optPar) == "log_tauEps")]
  upper = upper[-which(names(optPar) == "log_tauEps")]
  
  newStartTime = proc.time()[3]
  for(thisTol in tolSeq) {
    testObj = objFixed
    testObj$env$inner.control = list(maxit=1000, tol10=thisTol)
    testObj$env$tracepar = TRUE
    print(paste0("optimizing for tol = ", thisTol, "."))
    opt1 <- optim(par=thisOptPar, fn=funWrapper, gr=grWrapper,
                  method = c("BFGS"), hessian = FALSE, control=list(reltol=thisTol))
    # opt1 <- optim(par=optPar, fn=funWrapper, 
    #               hessian = FALSE, control=list(reltol=thisTol))
    thisOptPar = opt1$par
    if(!is.null(opt1$message)) {
      print(paste0("error for tol = ", thisTol, ". Message:"))
      print(opt1$message)
      next
    }
    else {
      print(paste0("completed optimization for tol = ", thisTol, ""))
      
      ## Get standard errors
      print("getting standard errors...")
      sdTime = system.time(
        SD0 <- TMB::sdreport(testObj, getJointPrecision=TRUE,
                             bias.correct = TRUE,
                             bias.correct.control = list(sd = TRUE))
      )[3]
      # SD0
      print(sdTime/60)
      # 3.9447 minutes for intern=FALSE
      
      if(SD0$pdHess) {
        print("Optimization and PD hess calculation done!")
        break
      }
      else {
        print("Hessan not PD. Rerunning optimization with stricter tol...")
        
      }
    }
    
    
    
  }
  thisEndTime = proc.time()[3]
  sdTime/60
  thisTotalTime = thisEndTime - newStartTime
  totalTime = thisEndTime - startTime
  print(paste0("second optimization took ", thisTotalTime/60, " minutes"))
  print(paste0("all optimization took ", totalTime/60, " minutes"))
}

# opt0 <- nlminb(start       =    obj[['par']],
#                objective   =    obj[['fn']],
#                gradient    =    obj[['gr']],
#                lower = rep(-10, length(obj[['par']])),
#                upper = rep( 10, length(obj[['par']])),
#                control     =    list(trace=1))

# * TMB Posterior Sampling ----

if(FALSE) {
  hessTime = system.time(testHess <- hessian(testObj$fn, optPar))[3]
  hessTime/60
  # 10.28792 minutes for intern=FALSE
  eig = eigen(testHess)
  eig
  
  for(i in 1:length(eig$values)) {
    thisVal = eig$values[i]
    thisVec = eig$vectors[,i]
    barplot(eig$vectors[,i], names.arg=names(SD0$par.fixed), 
            main=paste0("Eigenvector ", i, " with value ", thisVal))
  }
  
  testObj$fn(SD0$par.fixed)
  testObj$fn(SD0$par.fixed + eig$vectors[,9]*.05)
  testObj$fn(SD0$par.fixed - eig$vectors[,9]*.05)
}



## summary(SD0, 'report')
## summary(SD0, 'fixed')

save(SD0, obj, totalTime, sdTime, file="savedOutput/ed/fitM_D2.RData")
out = load("savedOutput/ed/fitM_D2.RData")

gridPreds = predGrid(SD0, obj, admLevel="adm2")
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & Est & Q0.025 & Q0.1 & Q0.9 & Q0.975 \\
# \hline
# (Int) & -2.00 & -2.22 & -2.15 & -1.86 & -1.79 \\
# urb & 0.28 & 0.07 & 0.13 & 0.42 & 0.49 \\
# access & -0.25 & -0.38 & -0.34 & -0.17 & -0.13 \\
# elev & 0.20 & 0.06 & 0.12 & 0.29 & 0.34 \\
# distRiversLakes & 0.07 & -0.05 & -0.01 & 0.16 & 0.20 \\
# popValsNorm & 0.90 & 0.66 & 0.74 & 1.07 & 1.16 \\
# sigmaSq & 0.65 & 0.49 & 0.53 & 0.77 & 0.85 \\
# phi & 0.13 & 0.10 & 0.11 & 0.16 & 0.17 \\
# sigmaEpsSq & 0.87 & 0.71 & 0.76 & 0.99 & 1.07 \\
# \hline
# \end{tabular}
# \end{table} 
save(gridPreds, file="savedOutput/ed/gridPredsM_D2.RData")
out = load("savedOutput/ed/gridPredsM_D2.RData")

stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
save(stratPreds, file="savedOutput/ed/stratPredsM_D2.RData")
save(admin1Preds, file="savedOutput/ed/admin1PredsM_D2.RData")
save(admin2Preds, file="savedOutput/ed/admin2PredsM_D2.RData")
out = load("savedOutput/ed/stratPredsM_D2.RData")
out = load("savedOutput/ed/admin1PredsM_D2.RData")
out = load("savedOutput/ed/admin2PredsM_D2.RData")

summaryTabBYM2(SD0, obj, popMat=popMatNGAThresh, 
               gridPreds=gridPreds)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrr}
# \hline
# & Est & Q0.025 & Q0.975 \\ 
# \hline
# (Int) & -2.01 & -2.22 & -2.14 & -1.86 & -1.80 \\
# urb & 0.28 & 0.05 & 0.14 & 0.42 & 0.48 \\
# access & -0.26 & -0.38 & -0.33 & -0.18 & -0.13 \\
# elev & 0.20 & 0.07 & 0.11 & 0.29 & 0.34 \\
# distRiversLakes & 0.07 & -0.06 & -0.01 & 0.15 & 0.19 \\
# popValsNorm & 0.90 & 0.66 & 0.74 & 1.06 & 1.13 \\
# sigmaSq & 0.65 & 0.49 & 0.53 & 0.76 & 0.83 \\
# phi & 0.13 & 0.10 & 0.11 & 0.16 & 0.17 \\
# sigmaEpsSq & 0.88 & 0.72 & 0.77 & 0.98 & 1.06 \\
# \hline
# \end{tabular}
# \end{table}
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=NULL, 
          plotNameRoot="edM_D2")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=stratPreds, 
          plotNameRoot="edM_D2", plotNameRootAreal="Strat")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin1Preds, 
          plotNameRoot="edM_D2", plotNameRootAreal="Admin1")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin2Preds, 
          plotNameRoot="edM_D2", plotNameRootAreal="Admin2")


