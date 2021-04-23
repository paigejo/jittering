# Section 0: Load packages, source code ----
library(raster)
library(TMB)
library(SUMMER)
library(fields)
library(INLA)
library(viridisLite)

# Section 1: Set up the data ----

# Load data
load("~/git/continuousNugget/savedOutput/global/kenyaData.RData")
# head(mort)
#   clusterID  n y regionRural  region urban samplingWeight      lon
# 1         1 10 1         9.1 Nairobi  TRUE        5476381 36.75296
# 2        10  9 0         9.1 Nairobi  TRUE        3067129 36.80952
# 3       100  4 0         4.2 Central FALSE         735830 37.14386
# 4      1000 26 1         8.2  Nyanza FALSE         742949 34.91504
# 5      1001  8 0         8.2  Nyanza FALSE         751588 34.86289
# 6      1002  9 1         8.2  Nyanza FALSE         379868 34.87003
#         lat  admin1      east      north
# 1 -1.282723 Nairobi 249.95384 -141.87640
# 2 -1.315103 Nairobi 256.25440 -145.45222
# 3 -0.430589   Nyeri 293.42387  -47.61379
# 4 -0.640455 Nyamira  45.07855  -70.96453
# 5 -0.663048 Nyamira  39.26281  -73.47274
# 6 -0.732350 Nyamira  40.06611  -81.15141
#    admin2
# 1                                                                                    Westlands + Dagoretti North
# 2 Kamukunji + Starehe + Langata + Dagoretti South + Kibra + Mathare + Embakasi North + Embakasi South + Makadara
# 3                                                                                                        Mathira
# 4                                                                                                 West Mugirango
# 5                                                                                                  Kitutu Masaba
# 6                                                                                                  Kitutu Masaba

# sample a small number of observations for testing purposes
N = 100
set.seed(123)
inds = sample(1:nrow(mort), N, replace=FALSE)
dat = mort[inds,]

coords = cbind(dat$east, dat$north) # projected easting/northing in km
urbanVals = dat$urban # TRUE/FALSE if urban/rural
ys = dat$y # population numerator
ns = dat$n # population denominator

# Section 2: Fitting and prediction with TMB ----
# * Compile TMB function and load it in ----
if(FALSE) {
  compile( "code/modSPDEJitter.cpp")
}
dyn.load( dynlib("code/modSPDEJitter"))

# * Prep inputs to TMB ----

# construct to the integration points and their weights
# scaling factor is factor by which to increase jittering distance relative to DHS
intPointInfo = makeAllIntegrationPoints(coords, urbanVals, 
                                        numPointsUrban=11, numPointsRural=16, 
                                        scalingFactor=1, 
                                        JInnerUrban=3, JOuterUrban=0, 
                                        JInnerRural=3, JOuterRural=1)
wUrban = intPointInfo$wUrban
wRural = intPointInfo$wRural
n_integrationPointsUrban = ncol(wUrban)
n_integrationPointsRural = ncol(wRural)

# construct SPDE mesh and prior: 
# jitter the dataset coordinates a bit so triangulation vertices aren't all 
# put on observation locations
# prior for range is set so median is 1/5 times minimum of east and north range
# prior for spatial standard deviation is set so P(sigma > 1) = 0.05
jitterAmount = 7 / 4
jitteredCoords = locs=cbind(jitter(dat$east, amount=jitterAmount), jitter(dat$north, amount=jitterAmount))
mesh.s = getSPDEMeshKenya(jitteredCoords)
USpatial = 1
alphaSpatial = 0.05
spde = getSPDEPrior(mesh.s, U=USpatial, alpha=alphaSpatial)
range0 = min(c(diff(range(mesh.s$loc[, 1])), diff(range(mesh.s$loc[, 2])))) # 1444.772
range0 = range0 / 5

# construct projection matrices, and get other relevant info for TMB
out = makeJitterDataForTMB(intPointInfo, ys, urbanVals, ns, spdeMesh=mesh.s)
ysUrban = out$ysUrban
ysRural = out$ysRural
nsUrban = out$nsUrban
nsRural = out$nsRural
AUrban = out$AUrban
ARural = out$ARural

## set priors
# normal vector of mean, sd for intercept
alpha.pri = c(0, 100)

# parameters for the Matern joint PC prior
# range limit:    rho0
# range prob:     alpha_rho
# field sd limit: sigma0
# field sd prob:  alpha_sigma
matern.pri = c(range0, 0.5, USpatial, alphaSpatial)

# compile inputs for TMB
data_full <- list(num_iUrban = length(ysUrban),  # Total number of urban observations
                  num_iRural = length(ysRural),  # Total number of rural observations
                  num_s = mesh.s[['n']], # num. of vertices in SPDE mesh
                  y_iUrban   = ysUrban, # num. of pos. urban obs in the cluster
                  y_iRural   = ysRural, # num. of pos. rural obs in the cluster
                  n_iUrban   = nsUrban,  # num. of urban exposures in the cluster
                  n_iRural   = nsRural,  # num. of rural exposures in the cluster
                  n_integrationPointsUrban = n_integrationPointsUrban, 
                  n_integrationPointsRural = n_integrationPointsRural, 
                  wUrban = wUrban, 
                  wRural = wRural, 
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

## Specify starting values for TMB params
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

# * Run TMB ----
opt0 <- nlminb(start       =    obj[['par']],
               objective   =    obj[['fn']],
               gradient    =    obj[['gr']],
               lower = rep(-10, length(obj[['par']])),
               upper = rep( 10, length(obj[['par']])),
               control     =    list(trace=1))

# * TMB Posterior Sampling ----
## Get standard errors
SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                     bias.correct = TRUE,
                     bias.correct.control = list(sd = TRUE))
## summary(SD0, 'report')
## summary(SD0, 'fixed')

## take samples from fitted model
mu <- c(SD0$par.fixed,SD0$par.random)

## simulate draws
rmvnorm_prec <- function(mu, chol_prec, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- chol_prec #Cholesky(prec, super=TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  mu + z
}

L <- Cholesky(SD0[['jointPrecision']], super = T)
t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = 500)

## summarize the draws
parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
epsilon_tmb_draws  <- t.draws[parnames == 'Epsilon_s',]
alpha_tmb_draws    <- matrix(t.draws[parnames == 'alpha',], nrow = 1)

# project from mesh to raster, add intercept
data(kenyaPopulationData)
predCoords = cbind(popMatKenya$east, popMatKenya$north)
A.pred = inla.spde.make.A(mesh = mesh.s, loc = predCoords)
pred_tmb <- as.matrix(A.pred %*% epsilon_tmb_draws)
pred_tmb <- sweep(pred_tmb, 2, alpha_tmb_draws, '+')

# convert to probability scale
pred_tmb = expit(pred_tmb)

## find the median and sd across draws, as well as 90% intervals
summ_tmb <- cbind(median = (apply(pred_tmb, 1, median)),
                  sd     = (apply(pred_tmb, 1, sd)),
                  lower = (apply(pred_tmb, 1, quantile, .05)),
                  upper = (apply(pred_tmb, 1, quantile, .95)))

# * Plot TMB results ----
## make summary rasters
# kenyaExtent = extent(cbind(popMatKenya$east, popMatKenya$north))
# pop.rast <- raster(nrows=length(unique(popMatKenya$north)), ncols=length(length(unique(popMatKenya$north))),
#                   xmn=min(popMatKenya$east), xmx=max(length(unique(popMatKenya$east))), 
#                   ymn=min(popMatKenya$north), ymx=max(popMatKenya$north),
#                   vals=popMatKenya$pop, ext = kenyaExtent)
pop.rast = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=popMatKenya$pop), 
                         res=c(5, 5))
ras_med_tmb = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_tmb[, 1]), 
                            res=c(5, 5))
ras_sdv_tmb = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_tmb[, 2]), 
                            res=c(5, 5))
ras_lower_tmb = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_tmb[, 3]), 
                            res=c(5, 5))
ras_upper_tmb = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_tmb[, 4]), 
                            res=c(5, 5))

## plot truth, pixels falling within/without the 90% interval,
##  post. median, and post sd

# set the range for the truth and median
rast.zrange <- range(c(dat$y/dat$n, values(ras_med_tmb)), na.rm = T)

# plot tmb
par(mfrow = c(2, 2))
quilt.plot(dat$east, dat$north, dat$y/dat$n, main = 'Observations', col = (viridis(100)), 
           nx=30, ny=30)
# plot(ras_inInt_tmb, main = 'Pixels where 90% CIs did not cover Truth')
# points(dat$east, dat$north, pch=".")
plot(ras_med_tmb, main = 'TMB Posterior Median',
     col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_tmb, main='TMB Posterior Standard Deviation')
points(dat$east, dat$north, pch=".")

## plot TMB meds and stdevs
med.zrange <- range(values(ras_med_tmb), na.rm = T)
sdv.zrange <- range(values(ras_sdv_tmb), na.rm = T)

par(mfrow = c(1, 2))
plot(ras_med_tmb, main = 'TMB Posterior Median',
     zlim = med.zrange, col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_tmb, main = 'TMB Posterior Standard Deviation',
     zlim = sdv.zrange)
points(dat$east, dat$north, pch=".")



# Section 3: Fitting and prediction with INLA ----

# * Prep inputs for INLA ----
design_matrix <- data.frame(int = rep(1, nrow(dat)))
A.proj = inla.spde.make.A(mesh = mesh.s, loc = cbind(dat$east, dat$north))
stack.obs <- inla.stack(tag='est',
                        data=list(Y = dat$y, ## response
                                  N = dat$n), ## binom trials
                        A=list(A.proj, ## A.proj for space
                               1),     ## 1 for design.mat
                        effects=list(
                          space = 1:mesh.s[['n']],
                          design_matrix))

## define the INLA model
formula <- formula(Y ~ -1 + int + f(space, model = spde))

# * Run INLA ----
i.fit <- inla(formula,
              data = inla.stack.data(stack.obs),
              control.predictor = list(A = inla.stack.A(stack.obs),
                                       compute = FALSE),
              control.fixed = list(expand.factor.strategy = 'inla',
                                   prec = list(default = 1 / alpha.pri[2] ^ 2)),
              control.inla = list(strategy = 'simplified.laplace',
                                  int.strategy = 'ccd'),
              control.compute=list(config = TRUE),
              family = 'binomial',
              Ntrials = N,
              verbose = FALSE,
              keep = FALSE)

# * INLA Posterior Sampling ----
# i.draws <- inla.posterior.sample(n = 500, i.fit,
#                                  use.improved.mean = TRUE,
#                                  skew.corr = TRUE)
i.draws <- inla.posterior.sample(n = 500, i.fit,
                                 use.improved.mean = TRUE)

## summarize the draws
par_names <- rownames(i.draws[[1]][['latent']])
s_idx <- grep('^space.*', par_names)
a_idx <- which(!c(1:length(par_names)) %in%
                 grep('^space.*|Predictor|clust.id', par_names))

# project from mesh to raster, add intercept
pred_s <- sapply(i.draws, function (x) x[['latent']][s_idx])
pred_inla <- as.matrix(A.pred %*% pred_s)
alpha_inla_draws <- sapply(i.draws, function (x) x[['latent']][a_idx])
pred_inla <- sweep(pred_inla, 2, alpha_inla_draws, '+')

# convert to probability scale
pred_inla = expit(pred_inla)

## find the median and sd across draws, as well as 90% intervals
summ_inla <- cbind(median = (apply(pred_inla, 1, median)),
                   sd     = (apply(pred_inla, 1, sd)),
                   lower = (apply(pred_inla, 1, quantile, .05)),
                   upper = (apply(pred_inla, 1, quantile, .95)))

# * Plot INLA results ----
## make summary rasters
ras_med_inla = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_inla[, 1]), 
                            res=c(5, 5))
ras_sdv_inla = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_inla[, 2]), 
                            res=c(5, 5))
ras_lower_inla = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_inla[, 3]), 
                              res=c(5, 5))
ras_upper_inla = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_inla[, 4]), 
                              res=c(5, 5))

## plot truth, pixels falling within/without the 90% interval,
##  post. median, and post sd

# set the range for the truth and median
# rast.zrange <- range(c(dat$y/dat$n, values(ras_med_inla)), na.rm = T)

# plot
par(mfrow = c(2, 2))
quilt.plot(dat$east, dat$north, dat$y/dat$n, main = 'Observations', col = (viridis(100)), nx=30, ny=30)
# points(dat[, .(x, y)])
# plot(ras_inInt_inla, main = 'Pixels where 90% CIs did not cover Truth')
# points(dat$east, dat$north, pch=".")
plot(ras_med_inla, main = 'INLA Posterior Median',
     col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_inla, main = 'INLA Posterior Standard Deviation')
points(dat$east, dat$north, pch=".")

# Section 4: Compare TMB and INLA ----

## compare INLA and TMB meds and stdevs
med.zrange <- range(c(values(ras_med_tmb), values(ras_med_inla)), na.rm = T)
sdv.zrange <- c(0, max(c(values(ras_sdv_tmb), values(ras_sdv_inla)), na.rm = T))

par(mfrow = c(2, 2))
plot(ras_med_inla, main = 'INLA Posterior Median',
     zlim = med.zrange, col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_inla, main = 'INLA Posterior Standard Deviation',
     zlim = sdv.zrange)
points(dat$east, dat$north, pch=".")
plot(ras_med_tmb, main = 'TMB Posterior Median',
     zlim = med.zrange, col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_tmb, main = 'TMB Posterior Standard Deviation',
     zlim = sdv.zrange)
points(dat$east, dat$north, pch=".")

par(mfrow = c(1, 2))
plot(ras_med_tmb, main = 'TMB Posterior Median',
     col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_tmb, main = 'TMB Posterior Standard Deviation')
points(dat$east, dat$north, pch=".")

par(mfrow = c(1, 2))
plot(ras_med_inla, main = 'INLA Posterior Median',
     col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_inla, main = 'INLA Posterior Standard Deviation')
points(dat$east, dat$north, pch=".")

summary(SD0, 'report')
summary(i.fit)

# Section 5: Simple (non-jittered) Example ----
# * Compile TMB function and load it in ----
if(FALSE) {
  compile( "code/tmb_spde.cpp")
}
dyn.load( dynlib("code/tmb_spde") )

# * Prep inputs to TMB ----
data_full <- list(num_i = nrow(dat),  # Total number of observations
                  num_s = mesh.s[['n']], # num. of vertices in SPDE mesh
                  y_i   = dat$y, # num. of pos. obs in the cluster
                  n_i   = dat$n,  # num. of exposures in the cluster
                  X_alpha  = matrix(1, nrow = nrow(dat), ncol = 1),# des.mat
                  M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
                  M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
                  M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
                  Aproj = A.proj,             # Projection matrix
                  options = c(1, ## if 1, use normalization trick
                              1), ## if 1, run adreport
                  # normalization flag.
                  flag = 1,
                  alpha_pri = alpha.pri, ## normal
                  matern_pri = matern.pri
)

## Specify starting values for TMB params
tmb_params <- list(alpha = 0.0, # intercept
                   log_tau = 0, # Log inverse of tau (Epsilon)
                   log_kappa = 0, # Matern range parameter
                   Epsilon_s = rep(0, mesh.s[['n']]) # RE on mesh vertices
)

## make a list of things that are random effects
rand_effs <- c('Epsilon_s')

## make the autodiff generated liklihood func & gradient
objSimple <- MakeADFun(data=data_full,
                 parameters=tmb_params,
                 random=rand_effs,
                 hessian=TRUE,
                 DLL='code/tmb_spde')

## we can normalize the GMRF outside of the nested optimization,
## avoiding unnecessary and expensive cholesky operations.
objSimple <- normalize(objSimple, flag="flag", value = 0)

# * Run TMB ----
opt0 <- nlminb(start       =    objSimple[['par']],
               objective   =    objSimple[['fn']],
               gradient    =    objSimple[['gr']],
               lower = rep(-10, length(objSimple[['par']])),
               upper = rep( 10, length(objSimple[['par']])),
               control     =    list(trace=1))

# * TMB Posterior Sampling ----
## Get standard errors
SD0 <- TMB::sdreport(objSimple, getJointPrecision=TRUE,
                     bias.correct = TRUE,
                     bias.correct.control = list(sd = TRUE))
## summary(SD0, 'report')
## summary(SD0, 'fixed')

## take samples from fitted model
mu <- c(SD0$par.fixed,SD0$par.random)

## simulate draws
rmvnorm_prec <- function(mu, chol_prec, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- chol_prec #Cholesky(prec, super=TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  mu + z
}

L <- Cholesky(SD0[['jointPrecision']], super = T)
t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = 500)

## summarize the draws
parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
epsilon_simple_draws  <- t.draws[parnames == 'Epsilon_s',]
alpha_simple_draws    <- matrix(t.draws[parnames == 'alpha',], nrow = 1)

# project from mesh to raster, add intercept
pred_simple <- as.matrix(A.pred %*% epsilon_simple_draws)
pred_simple <- sweep(pred_simple, 2, alpha_simple_draws, '+')

# convert to probability scale
pred_simple = expit(pred_simple)

## find the median and sd across draws, as well as 90% intervals
summ_simple <- cbind(median = (apply(pred_simple, 1, median)),
                  sd     = (apply(pred_simple, 1, sd)),
                  lower = (apply(pred_simple, 1, quantile, .05)),
                  upper = (apply(pred_simple, 1, quantile, .95)))

# * Plot TMB results ----
## make summary rasters
ras_med_simple = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_simple[, 1]), 
                            res=c(5, 5))
ras_sdv_simple = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_simple[, 2]), 
                            res=c(5, 5))
ras_lower_simple = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_simple[, 3]), 
                              res=c(5, 5))
ras_upper_simple = rasterFromXYZ(data.frame(x=popMatKenya$east, y=popMatKenya$north, z=summ_simple[, 4]), 
                              res=c(5, 5))

## plot truth, pixels falling within/without the 90% interval,
##  post. median, and post sd

# set the range for the truth and median
rast.zrange <- range(c(values(ras_med_simple)), na.rm = T)

# plot tmb
par(mfrow = c(1, 2))
plot(ras_med_simple, main = 'Simple (non-jittered) TMB Posterior Median',
     col = (viridis(100)))
points(dat$east, dat$north, pch=".")
plot(ras_sdv_simple, main='Simple (non-jittered) TMB Posterior Standard Deviation')
points(dat$east, dat$north, pch=".")

