# some tests


nUrban = length(ysUrban)
is = 1:nUrban
par(mfrow=c(2,2))
plot(coordsUrban[is,], pch=".")
plot(coordsUrban[is+nUrban,], pch=".")
plot(coordsUrban[is+nUrban*2,], pch=".")
# plot(coordsUrban[is+nUrban*3,], pch=".")

is = apply(AUrban!= 0, 1, which)
coords = mesh.s$loc[is,1:2]
plot(coords, pch=".", col="blue")

nRural = length(ysRural)
is = 1:nRural
par(mfrow=c(2,2))
plot(coordsRural[is,], pch=".")
plot(coordsRural[is+nRural,], pch=".")
plot(coordsRural[is+nRural*2,], pch=".")
# plot(coordsRural[is+nRural*3,], pch=".")

is = apply(ARural!= 0, 1, which)
coords = mesh.s$loc[is,1:2]
plot(coords, pch=".", col="blue")


wUrbanTest = matrix(c(rep(1, nrow(wUrban)), rep(0, nrow(wUrban) * (ncol(wUrban) - 1))), ncol=ncol(wUrban))
wRuralTest = matrix(c(rep(1, nrow(wRural)), rep(0, nrow(wRural) * (ncol(wRural) - 1))), ncol=ncol(wRural))
image.plot(wRural)
image.plot(wRuralTest)

data_full <- list(num_iUrban = length(ysUrban),  # Total number of urban observations
                  num_iRural = length(ysRural),  # Total number of rural observations
                  num_s = mesh.s[['n']], # num. of vertices in SPDE mesh
                  y_iUrban   = ysUrban, # num. of pos. urban obs in the cluster
                  y_iRural   = ysRural, # num. of pos. rural obs in the cluster
                  n_iUrban   = nsUrban,  # num. of urban exposures in the cluster
                  n_iRural   = nsRural,  # num. of rural exposures in the cluster
                  n_integrationPointsUrban = n_integrationPointsUrban, 
                  n_integrationPointsRural = n_integrationPointsRural, 
                  wUrban = wUrbanTest, 
                  wRural = wRuralTest, 
                  X_alphaUrban  = matrix(1, nrow = length(ysUrban), ncol = 1),# des.mat for urban observations
                  X_alphaRural  = matrix(1, nrow = length(ysRural), ncol = 1),# des.mat for rural observations
                  M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
                  M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
                  M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
                  AprojUrban = AUrban,             # Projection matrix (urban)
                  AprojRural = ARural,             # Projection matrix (rural)
                  options = c(1, ## if 1, use normalization trick
                              1), ## if 1, run adreport
                  # normalization flag.
                  flag = 1,
                  alpha_pri = alpha.pri, ## normal vector of mean, sd for intercept
                  matern_pri = matern.pri ## 
)

tmb_params <- list(alpha = 0.0, # intercept
                   log_tau = 0, # Log tau (i.e. log spatial precision, Epsilon)
                   log_kappa = 0, # SPDE parameter related to the range
                   Epsilon_s = rep(0, mesh.s[['n']]) # RE on mesh vertices
)

## make a list of things that are random effects
rand_effs <- c('Epsilon_s')

## make the autodiff generated likelihood func & gradient
obj <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='modSPDEJitter')

## we can normalize the GMRF outside of the nested optimization,
## avoiding unnecessary and expensive cholesky operations.
obj <- normalize(obj, flag="flag", value = 0)

obj[['par']]
objSimple[['par']]
obj[['fn']](obj[['par']])
objSimple[['fn']](objSimple[['par']])

# urbanInds = 1 + c(outer(which(dat$urban)-1, (0:(length(ysUrban)-1))*length(ysUrban), "+"))
A.projUrban = A.proj[dat$urban,]
testA.projUrban = AUrban[1:length(ysUrban),]
image.plot(as.matrix(A.projUrban))
image.plot(as.matrix(testA.projUrban))
all.equal(as.matrix(A.projUrban), as.matrix(testA.projUrban))
