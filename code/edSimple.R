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
KMICS=300
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
  load(paste0("savedOutput/global/intPtsMICS_", KMICS, "_adm2Cov.RData"))
  
  # AUrbDHS = makeApointToArea(intPtsDHS$areasUrban, admFinal$NAME_FINAL) # 41 x 569 nStrat x nObsUrb
  # ARurDHS = makeApointToArea(intPtsDHS$areasRural, admFinal$NAME_FINAL) # 41 x 810
  
  AUrbDHS = makeApointToArea(ed$subarea[ed$urban], adm2$NAME_2) # 775 x 569 nArea x nObsUrb
  ARurDHS = makeApointToArea(ed$subarea[!ed$urban], adm2$NAME_2) # 775 x 810
  
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
  intPtsDHS$covsUrb = intPtsDHS$covsUrb[1:nrow(AUrbDHS),-1] # don't include intercepts
  intPtsDHS$covsRur = intPtsDHS$covsRur[1:nrow(ARurDHS),-1]
  
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
       file="savedOutput/global/edSimpleInputs.RData")
  
  # compile model ----
  # dyn.unload( dynlib("code/modBYM2JitterDHS2sparse"))
  # compile( "code/modBYM2JitterDHS2sparse.cpp", framework="TMBad", safebounds=FALSE)
  
  dyn.unload( dynlib("code/modBYM2simple"))
  compile( "code/modBYM2simple.cpp", 
           framework="TMBad", safebounds=FALSE)
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/TMB/include" -I"/Library/Frameworks/R.framework/Versions/4.2/Resources/library/RcppEigen/include"  -DTMB_SAFEBOUNDS -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modBYM2JitterDHS2  -DTMB_LIB_INIT=R_init_modBYM2JitterDHS2  -DTMBAD_FRAMEWORK  -I/usr/local/include   -fPIC  -Wall -g -O2  -c code/modBYM2JitterDHS2.cpp -o code/modBYM2JitterDHS2.o
  # clang++ -mmacosx-version-min=10.13 -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o code/modBYM2JitterDHS2.so code/modBYM2JitterDHS2.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
  
  # on Idun:
  # g++ -std=gnu++14 -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/include" -DNDEBUG -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/TMB/include" -I"/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/library/RcppEigen/include"   -DTMB_EIGEN_DISABLE_WARNINGS -DLIB_UNLOAD=R_unload_modBYM2JitterDHS2  -DTMB_LIB_INIT=R_init_modBYM2JitterDHS2  -DTMBAD_FRAMEWORK  -I/cluster/apps/eb/software/OpenSSL/1.1/include -I/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/include -I/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/include -I/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/include -I/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/include -I/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/include -I/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Java/11.0.2/include -I/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/include -I/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/include -I/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/include -I/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/include -I/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/include -I/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/include -I/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/include -I/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/include -I/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include -I/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/include/flexiblas   -fpic  -O2 -ftree-vectorize -march=native -fno-math-errno  -c code/modBYM2JitterDHS2.cpp -o code/modBYM2JitterDHS2.o
  # g++ -std=gnu++14 -shared -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -L/cluster/apps/eb/software/OpenSSL/1.1/lib64 -L/cluster/apps/eb/software/OpenSSL/1.1/lib -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libgit2/1.4.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/MPFR/4.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib64 -L/cluster/apps/eb/software/GDAL/3.5.0-foss-2022a/lib -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/nodejs/16.15.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GLPK/5.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ImageMagick/7.1.0-37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/GSL/2.7-GCC-11.3.0/lib -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/UDUNITS/2.2.28-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib64 -L/cluster/apps/eb/software/HDF5/1.12.2-gompi-2022a/lib -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ICU/71.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libsndfile/1.1.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FFTW/3.3.10-GCC-11.3.0/lib -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/NLopt/2.7.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/GMP/6.2.1-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libxml2/2.9.13-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cURL/7.83.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Tk/8.6.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Java/11.0.2/lib64 -L/cluster/apps/eb/software/Java/11.0.2/lib -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/LibTIFF/4.3.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libjpeg-turbo/2.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libpng/1.6.37-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/PCRE2/10.40-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/SQLite/3.38.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/zlib/1.2.12-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/XZ/5.2.5-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/bzip2/1.0.8-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/ncurses/6.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libreadline/8.1.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/cairo/1.17.4-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/libGLU/9.0.2-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Mesa/22.0.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/X11/20220504-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/Xvfb/21.1.3-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib64 -L/cluster/apps/eb/software/pkgconf/1.8.0-GCCcore-11.3.0/lib -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib64 -L/cluster/apps/eb/software/ScaLAPACK/2.2.0-gompi-2022a-fb/lib -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib64 -L/cluster/apps/eb/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib -L/cluster/apps/eb/software/GCCcore/11.3.0/lib64 -L/cluster/apps/eb/software/GCCcore/11.3.0/lib -o code/modBYM2JitterDHS2.so code/modBYM2JitterDHS2.o -L/cluster/apps/eb/software/R/4.2.1-foss-2022a/lib64/R/lib -lR
}

# load in TMB function inputs
out = load("savedOutput/global/edSimpleInputs.RData")

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
rand_effs <- c('Epsilon_bym2', 'nuggetUrbDHS', 'nuggetRurDHS', 'alpha', 'beta')
# rand_effs <- c('Epsilon_bym2', 'nuggetUrbDHS', 'nuggetRurDHS')

# collect input data

data_full = list(
  y_iUrbanDHS=ysUrbDHS, # same as above but for DHS survey
  y_iRuralDHS=ysRurDHS, # 
  n_iUrbanDHS=nsUrbDHS, # number binomial trials
  n_iRuralDHS=nsRurDHS, # 
  AprojUrbanDHS=AUrbDHS, # [nObsUrban] x nArea matrix with ij-th entry = 1 if cluster i associated with area j and 0 o.w.
  AprojRuralDHS=ARurDHS, # 
  X_betaUrbanDHS=intPtsDHS$covsUrb, # [nObsUrban] x nPar design matrix
  X_betaRuralDHS=intPtsDHS$covsRur, # 
  
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
  sapply(data_full, class)
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

if(FALSE) {
  tmb_params <- obj$env$last.par
  
  tmb_params <- list(alpha = SD0$par.random[1], # intercept
                     beta = SD0$par.random[2:6], 
                     log_tau = optPar[1], # Log tau (i.e. log spatial precision, Epsilon)
                     logit_phi = optPar[2], # SPDE parameter related to the range
                     log_tauEps = optPar[3], # Log tau (i.e. log spatial precision, Epsilon)
                     Epsilon_bym2 = rep(0, ncol(bym2ArgsTMB$Q)), # RE on mesh vertices
                     nuggetUrbDHS = rep(0, length(data_full$y_iUrbanDHS)), 
                     nuggetRurDHS = rep(0, length(data_full$y_iRuralDHS))
  )
}

# make TMB fun and grad ----
# dyn.load( dynlib("code/modBYM2JitterDHS2sparse"))
dyn.load( dynlib("code/modBYM2simple"))
TMB::config(tmbad.sparse_hessian_compress = 1)
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modBYM2simple')
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
  
  # test cases
  testPar = c(-0.0907214120323802, 0.409165119671409, -0.358399195083911)
  obj$fn(rep(0, 3))
  obj$fn(rep(0.5, 3))
  obj$fn(rep(1, 3))
  obj$fn(optParINLA)
  obj$fn(testPar)
  
  testRep = obj$report(obj$env$last.par)
  testRep2 = obj$report(obj$env$last.par)
  
  # compare densities
  testRep$bym2LogLik
  testRep2$bym2LogLik
  
  testRep$logPriPhi
  testRep2$logPriPhi
  
  -sum(testRep$liksUrbDHS)
  -sum(testRep2$liksUrbDHS)
  
  -sum(testRep$liksRurDHS)
  -sum(testRep2$liksRurDHS)
  
  testRep$logPriTau
  testRep2$logPriTau
  
  testRep$logPriTauEps
  testRep2$logPriTauEps
  
  testRep$bym2LogLik
  
  logTau = obj$env$last.par[names(obj$env$last.par) == "log_tau"]
  tau = exp(logTau)
  logitPhi = obj$env$last.par[names(obj$env$last.par) == "logit_phi"]
  phi = expit(logitPhi)
  testRep$bym2LogLik
  dBYM2naive(x=obj$env$last.par[grepl("Epsilon", names(obj$env$last.par))], 
             phi=phi, 
             tau=tau, 
             precomputedValues=bym2ArgsTMB)
  
  testRep$logPriTau
  logTau = obj$env$last.par[names(obj$env$last.par) == "log_tau"]
  tau = exp(logTau)
  lambdaTau = data_full$lambdaTau
  logTauEps = obj$env$last.par[names(obj$env$last.par) == "log_tauEps"]
  tauEps = exp(logTauEps)
  lambdaTauEps = data_full$lambdaTauEps
  
  dPCprec(tau=tau, lambda=lambdaTau, doLog=TRUE) + logTau
  testRep$logPriTau
  
  dPCprec(tau=tauEps, lambda=lambdaTauEps, doLog=TRUE) + logTauEps
  testRep$logPriTauEps
  
  logitPhi = obj$env$last.par[names(obj$env$last.par) == "logit_phi"]
  phi = expit(logitPhi)
  lambdaPhi = data_full$lambdaPhi
  dBYM2phiPC(phi=logitPhi, lambda=lambdaPhi, logitPhi=TRUE, Q=data_full$Q_bym2, 
             gammaTildesm1=data_full$gammaTildesm1, tr=data_full$tr, doLog=TRUE)
  testRep$logPriPhi
  
  nuggets = obj$env$last.par[grepl("nugget", names(obj$env$last.par))]
  logTauEps = obj$env$last.par[names(obj$env$last.par) == "log_tauEps"]
  tauEps = exp(logTauEps)
  sigmaEps = sqrt(1/tauEps)
  testRep$nuggetLogLik
  sum(dnorm(nuggets, sd=sigmaEps, log=TRUE))
  
  alpha = obj$env$last.par[grepl("alpha", names(obj$env$last.par))]
  betas = obj$env$last.par[grepl("beta", names(obj$env$last.par))]
  Epsilon = obj$env$last.par[grepl("Epsilon", names(obj$env$last.par))]
  nuggetsUrb = obj$env$last.par[grepl("nuggetUrbDHS", names(obj$env$last.par))]
  nuggetsRur = obj$env$last.par[grepl("nuggetRurDHS", names(obj$env$last.par))]
  
  # goodIndsUrb = seq(1, nrow(data_full$AprojUrbanDHS), by=11)
  # goodIndsRur = seq(1, nrow(data_full$AprojRuralDHS), by=16)
  goodIndsUrb = 1:length(nuggetsUrb)
  goodIndsRur = 1:length(nuggetsRur)
  logitPurb = data_full$AprojUrbanDHS[goodIndsUrb,] %*% Epsilon + data_full$X_betaUrbanDHS[goodIndsUrb,] %*% betas + nuggetsUrb + alpha
  logitPrur = data_full$AprojRuralDHS[goodIndsRur,] %*% Epsilon + data_full$X_betaRuralDHS[goodIndsRur,] %*% betas + nuggetsRur + alpha
  pUrb = c(expit(as.matrix(logitPurb)))
  pRur = c(expit(as.matrix(logitPrur)))
  sum(dbinom(data_full$y_iUrbanDHS, data_full$n_iUrbanDHS, prob=pUrb, log=FALSE))
  sum(dbinom(data_full$y_iRuralDHS, data_full$n_iRuralDHS, prob=pRur, log=FALSE))
  sum(testRep$liksUrbDHS)
  sum(testRep$liksRurDHS)
  
  system.time(test <- obj$fn(obj$par))
  # 363.236  for non-sparse
  # 355.211 for sparse
  
  system.time(test <- obj$gr(obj$par))
  # 54.5 for non-sparse
  # 55.3 for sparse
  
  testRep = obj$report()
  
  range(testRep$latentFieldUrbDHS)
  range(testRep$latentFieldRurDHS)
  
  range(testRep$fe_iUrbanDHS)
  range(testRep$fe_iRuralDHS)
  
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
# set start of optimization
# obj$par = optParINLA
# obj$par = optPar

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
  # optimization took 29.2183499999999 minutes (for intern=FALSE)
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
                        DLL='modBYM2simple')
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

save(SD0, obj, totalTime, sdTime, file="savedOutput/ed/fitSimple2.RData")
out = load("savedOutput/ed/fitSimple2.RData")

# gridPreds = predGrid(SD0, popMat=popMatNGAThresh, nsim=nsim, admLevel="adm2", 
#                      predAtArea=foldArea,
#                      quantiles=c(0.025, 0.1, 0.9, 0.975))
# preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
# preds$fixedMat = gridPreds$fixedMat
gridPreds = predGrid(SD0, popMat=popMatNGAThresh, nsim=1000, admLevel="adm2", 
                     quantiles=c(0.025, 0.1, 0.9, 0.975))
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & Est & Q0.025 & Q0.1 & Q0.9 & Q0.975 \\ 
# \hline
# (Int) & 0.01 & -0.19 & -0.11 & 0.12 & 0.19 \\ 
# urb & 3.06 & 2.84 & 2.92 & 3.21 & 3.28 \\ 
# access & 0.08 & -0.02 & 0.02 & 0.15 & 0.19 \\ 
# elev & -1.76 & -1.93 & -1.86 & -1.66 & -1.61 \\ 
# distRiversLakes & 0.30 & 0.14 & 0.19 & 0.40 & 0.45 \\ 
# popValsNorm & -0.19 & -0.40 & -0.32 & -0.05 & 0.02 \\ 
# sigmaSq & 0.85 & 0.63 & 0.70 & 1.01 & 1.11 \\ 
# phi & 0.50 & 0.13 & 0.22 & 0.79 & 0.89 \\ 
# sigmaEpsSq & 0.45 & 0.36 & 0.39 & 0.50 & 0.53 \\
# \hline
# \end{tabular}
# \end{table}
save(gridPreds, file="savedOutput/ed/gridPredsSimple2.RData")
out = load("savedOutput/ed/gridPredsSimple2.RData")

stratPreds = predArea(gridPreds, areaVarName="stratumMICS", orderedAreas=admFinal@data$NAME_FINAL)
admin1Preds = predArea(gridPreds, areaVarName="area", orderedAreas=adm1@data$NAME_1)
admin2Preds = predArea(gridPreds, areaVarName="subarea", orderedAreas=adm2@data$NAME_2)
save(stratPreds, file="savedOutput/ed/stratPredsSimple2.RData")
save(admin1Preds, file="savedOutput/ed/admin1PredsSimple2.RData")
save(admin2Preds, file="savedOutput/ed/admin2PredsSimple2.RData")
out = load("savedOutput/ed/stratPredsSimple2.RData")
out = load("savedOutput/ed/admin1PredsSimple2.RData")
out = load("savedOutput/ed/admin2PredsSimple2.RData")

summaryTabBYM2(SD0, obj, popMat=popMatNGAThresh, 
               gridPreds=gridPreds)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrr}
# \hline
# & Est & Q0.025 & Q0.1 & Q0.9 & Q0.975 \\
# \hline
# X.Int. & -1.76 & -1.91 & -1.61 \\ 
# beta & 0.30 & 0.14 & 0.45 \\ 
# beta.1 & -0.19 & -0.41 & 0.02 \\ 
# beta.2 & 0.18 & -0.10 & 0.49 \\ 
# beta.3 & 0.03 & -1.85 & 1.88 \\ 
# beta.4 & 0.81 & 0.61 & 1.01 \\ 
# sigmaSq & 0.99 & 0.81 & 1.21 \\ 
# phi & 0.95 & 0.94 & 0.96 \\ 
# sigmaEpsSq & 0.92 & 0.82 & 1.02 \\ 
# \hline
# \end{tabular}
# \end{table}
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=NULL, 
          plotNameRoot="edSimple2")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=stratPreds, 
          plotNameRoot="edSimple2", plotNameRootAreal="Strat")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin1Preds, 
          plotNameRoot="edSimple2", plotNameRootAreal="Admin1")
plotPreds(SD0, obj, popMat=popMatNGAThresh, 
          gridPreds=gridPreds, arealPreds=admin2Preds, 
          plotNameRoot="edSimple2", plotNameRootAreal="Admin2")

if(FALSE) {
  logTau = obj$env$last.par[names(obj$env$last.par) == "log_tau"]
  tau = exp(logTau)
  logitPhi = obj$env$last.par[names(obj$env$last.par) == "logit_phi"]
  phi = expit(logitPhi)
  plotBYM2priorCIwidth80(tau=tau, phi=phi, bym2ArgsTMB=bym2ArgsTMB, nsim=1000)
}

