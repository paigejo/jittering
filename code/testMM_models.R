# Test script for M_M models (MICS-only) using simStudy1, simulation i=1
# Compares: SepRepar (original), Marg (marginalized), Constr (constrained 2n-2)
# Runs standard TMB optimization AND TMBstan MCMC

# ---- Tracking execution status ----
executionLog = list(
  startTime = Sys.time(),
  optimizationStatus = list(),
  mcmcStatus = list(),
  errors = list()
)

# ---- Setup ----
Sys.setenv(HOME = "C:/Users/jpaige")
setwd("C:/Users/jpaige/git/jittering")
source("code/setup.R")
source("code/modM_MSepMarg.R")
source("code/modM_MSepConstr.R")

# Load TMBstan if available
if(!require("tmbstan", quietly=TRUE)) {
  cat("Warning: tmbstan package not available - MCMC results will be skipped\n")
  doMCMC = FALSE
  executionLog$mcmcStatus$available = FALSE
} else {
  doMCMC = TRUE
  executionLog$mcmcStatus$available = TRUE
}

# ---- Which models to run ----
runSepRepar = TRUE
runMarg     = TRUE
runConstr   = TRUE

# ---- Regenerate all results? (set TRUE to ignore progress file and rerun all models) ----
rerun = FALSE  # Set to TRUE to regenerate all results from scratch

# If rerun=TRUE, force all models to run
if(rerun) {
  runSepRepar = TRUE
  runMarg     = TRUE
  runConstr   = TRUE
}

# ---- Get standard errors? (sdreport is slow; can skip and just compare estimates) ----
getSDs = FALSE  # Set to TRUE to get full results including random effects

# ---- True parameters from simData1() defaults ----
# beta0=-1.25, gamma=1, betaRest=c(0, 0, 0, .5), margVar=.5, sigmaEpsilon=sqrt(1.5)
# In the model: alpha = beta0, beta = c(gamma, betaRest)
# sigmaSq = margVar = 0.5, sigmaEpsSq = sigmaEpsilon^2 = 1.5
# phi is not directly set in sim (sim uses SPDE, not BYM2), but truePar from runSimStudy1I:
truePar = c(-1.25, 1, 0, 0, 0, 0.5)
trueParNames = c("(Int)", "Urban", "Healthcare inaccessibility", 
                 "Elevation", "Dist. rivers & lakes", "Population")
trueSigmaEpsSq = 1.5
trueSigmaSq = 0.5

cat("True parameters:\n")
print(data.frame(name=trueParNames, value=truePar))
cat("True sigmaSq:", trueSigmaSq, "\n")
cat("True sigmaEpsSq:", trueSigmaEpsSq, "\n\n")

# ---- Load simulated data with validation ----
cat("Loading simulated data...\n")

# Check and load main data file
dataFile = "savedOutput/simStudy1/simPopsSurveys.RData"
if(!file.exists(dataFile)) {
  stop("ERROR: Data file not found: ", dataFile)
}
out = load(dataFile)
cat("Loaded from", dataFile, ":", paste(out, collapse=", "), "\n")

# Validate loaded objects
if(!exists("surveysDHS") || !exists("surveysMICS")) {
  stop("ERROR: surveysDHS or surveysMICS not found in loaded data")
}
if(length(surveysDHS) == 0 || length(surveysMICS) == 0) {
  stop("ERROR: surveysDHS or surveysMICS are empty")
}

i = 1
thisDHS = surveysDHS[[i]]
thisMICS = surveysMICS[[i]]
cat("Simulation index i =", i, "\n")
cat("  DHS survey: n =", nrow(thisDHS), "\n")
cat("  MICS survey: n =", nrow(thisMICS), "\n")

# Load DHS integration points (needed by makeInputsMDM even for M_M)
dhsIntFile = paste0("savedOutput/simStudy1/intPtsDHS_simStudy1_", i, ".RData")
if(!file.exists(dhsIntFile)) {
  stop("ERROR: DHS integration points file not found: ", dhsIntFile)
}
out = load(dhsIntFile)
cat("Loaded DHS integration points from:", dhsIntFile, "\n")

# Validate integration points
if(!exists("intPtsDHS")) {
  stop("ERROR: intPtsDHS not found in loaded data")
}
cat("DHS integration points: n =", nrow(intPtsDHS), "\n\n")

# ---- Progress tracking: load or initialize ----
progressFile = "savedOutput/testMM_progress.RData"
if(rerun) {
  cat("Rerun mode enabled - deleting existing progress and starting fresh\n")
  if(file.exists(progressFile)) {
    file.remove(progressFile)
    cat("Deleted:", progressFile, "\n")
  }
  allResults = list()
  cat("Starting fresh run\n\n")
} else if(file.exists(progressFile)) {
  cat("Loading progress from:", progressFile, "\n")
  load(progressFile)
  cat("Restored results for models:", names(allResults), "\n\n")
} else {
  allResults = list()
  cat("Starting fresh (no progress file found)\n\n")
}

# ---- Helper to summarize results ----
summarizeResults = function(SD0, obj, modelName) {
  cat("\n========================================\n")
  cat("Results for:", modelName, "\n")
  cat("========================================\n")
  
  # Check if SD0 is an sdreport object or a report list
  isSDReport = !is.null(SD0) && inherits(SD0, "sdreport")
  isReportList = !is.null(SD0) && is.list(SD0) && !inherits(SD0, "sdreport")
  
  # Get fixed parameters
  if(isSDReport) {
    parFixed = SD0$par.fixed
    if(!SD0$pdHess) {
      cat("WARNING: Hessian is not positive definite!\n")
    }
  } else if(isReportList) {
    # Extract fixed parameters from report list and obj$par
    cat("(Report object - using point estimates only, no standard errors)\n")
    # Construct parFixed from obj$par (fixed effects) and SD0 (hyperparameters from REPORT)
    parFixed = c(
      obj$par,  # alpha, beta, log_tau, logit_phi, log_tauEps
      log_tau = as.numeric(SD0$log_tau),
      logit_phi = as.numeric(SD0$logit_phi),
      log_tauEps = as.numeric(SD0$log_tauEps)
    )
    # Remove duplicates, keeping the reported values
    parFixed = parFixed[!duplicated(names(parFixed), fromLast=TRUE)]
  } else {
    cat("(SD0 is NULL - using point estimates only, no standard errors)\n")
    # Extract only numeric parameters from obj$par
    parFixed = obj$par
    # Ensure it's numeric
    parFixed = as.numeric(parFixed)
    # Try to get parameter names from obj if available
    if(!is.null(names(obj$par))) {
      names(parFixed) = names(obj$par)
    }
  }
  
  # Validate parFixed is numeric
  if(!is.numeric(parFixed)) {
    parFixed = as.numeric(parFixed)
  }
  
  parNames = names(parFixed)
  cat("\nFixed parameters:\n")
  tryCatch({
    print(round(parFixed, 4))
  }, error = function(e) {
    cat("Warning: Could not display parameters nicely:", conditionMessage(e), "\n")
    print(parFixed)
  })
  
  # Get random effects summary
  hasRandomEffects = FALSE
  if(isSDReport) {
    parRandom = SD0$par.random
    randomNames = names(parRandom)
    hasRandomEffects = TRUE
    cat("\nRandom effect counts:\n")
    print(table(randomNames))
  } else {
    cat("\n(Skipping random effect summaries - standard errors not computed)\n")
  }
  
  # Transform to interpretable scale
  log_tau = tryCatch(as.numeric(parFixed[parNames == "log_tau"]), error=function(e) NA)
  logit_phi = tryCatch(as.numeric(parFixed[parNames == "logit_phi"]), error=function(e) NA)
  log_tauEps = tryCatch(as.numeric(parFixed[parNames == "log_tauEps"]), error=function(e) NA)
  
  sigmaSq_est = if(!is.na(log_tau)) 1/exp(log_tau) else NA
  phi_est = if(!is.na(logit_phi)) plogis(logit_phi) else NA
  sigmaEpsSq_est = if(!is.na(log_tauEps)) 1/exp(log_tauEps) else NA
  
  cat("\nTransformed hyperparameters:\n")
  if(!is.na(sigmaSq_est)) {
    cat("  sigmaSq (1/tau):", round(sigmaSq_est, 4), "(true:", trueSigmaSq, ")\n")
  } else {
    cat("  sigmaSq: NOT AVAILABLE\n")
  }
  if(!is.na(phi_est)) {
    cat("  phi:", round(phi_est, 4), "\n")
  } else {
    cat("  phi: NOT AVAILABLE\n")
  }
  if(!is.na(sigmaEpsSq_est)) {
    cat("  sigmaEpsSq (1/tauEps):", round(sigmaEpsSq_est, 4), "(true:", trueSigmaEpsSq, ")\n")
  } else {
    cat("  sigmaEpsSq: NOT AVAILABLE\n")
  }
  
  # Check if model has alpha parameter (only if we have random effects)
  if(hasRandomEffects) {
    hasAlpha = "alpha" %in% names(parRandom)
    if(hasAlpha) {
      alpha_matches = which(randomNames == "alpha")
      if(length(alpha_matches) > 0) {
        alpha_est = parRandom[alpha_matches[1]]
        tryCatch({
          cat("  alpha (intercept):", round(as.numeric(alpha_est), 4), "(true:", truePar[1], ")\n")
        }, error = function(e) {
          cat("  alpha (intercept): extraction failed\n")
        })
      }
    }
    
    # Beta estimates
    beta_matches = which(randomNames == "beta")
    if(length(beta_matches) > 0) {
      beta_est = parRandom[beta_matches]
      cat("\nBeta estimates vs truth:\n")
      tryCatch({
        betaComp = data.frame(
          name = trueParNames[2:length(trueParNames)],
          truth = truePar[2:length(truePar)],
          estimate = round(as.numeric(beta_est), 4)
        )
        print(betaComp)
      }, error = function(e) {
        cat("Warning: Could not format beta comparison:", conditionMessage(e), "\n")
        print(beta_est)
      })
    } else {
      cat("\nBeta estimates: not found in random effects\n")
      beta_est = NA
    }
    
    invisible(list(sigmaSq=sigmaSq_est, phi=phi_est, sigmaEpsSq=sigmaEpsSq_est,
                   beta=beta_est, parFixed=parFixed))
  } else {
    cat("(Random effects not available - set getSDs=TRUE for full estimates)\n")
    invisible(list(sigmaSq=sigmaSq_est, phi=phi_est, sigmaEpsSq=sigmaEpsSq_est,
                   beta=NA, parFixed=parFixed))
  }
}


# ---- Run SepRepar model (original) ----
if(runSepRepar && (!('SepRepar' %in% names(allResults)) || rerun)) {
  cat("\n############################\n")
  cat("Running M_M SepRepar model\n")
  cat("############################\n")
  
  tryCatch({
    outSepRepar = fitMM(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, 
                         repar=TRUE, doMCMC=FALSE, MdInit=FALSE)
    
    resSepRepar = summarizeResults(outSepRepar$TMBsd, outSepRepar$TMBobj, "M_M SepRepar")
    cat("\nTotal time:", round(outSepRepar$totalTime/60, 2), "minutes\n")
    
    allResults$SepRepar = resSepRepar
    save(allResults, file=progressFile)
    cat("Progress saved to:", progressFile, "\n")
  }, error = function(e) {
    cat("\nERROR in SepRepar model:\n")
    cat(conditionMessage(e), "\n")
    cat("Skipping SepRepar model and continuing with others\n")
  })
} else if(runSepRepar && 'SepRepar' %in% names(allResults)) {
  cat("\nSkipping SepRepar model (already completed)\n")
}


# ---- Run Marg model (marginalized) ----
if(runMarg && (!('Marg' %in% names(allResults)) || rerun)) {
  cat("\n############################\n")
  cat("Running M_M Marg model\n")
  cat("############################\n")
  
  tryCatch({
    outMarg = fitMMMarg(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, 
                        doMCMC=FALSE, getSDs=getSDs)
    
    resMarg = summarizeResults(outMarg$TMBsd, outMarg$TMBobj, "M_M Marg")
    cat("\nTotal time:", round(outMarg$totalTime/60, 2), "minutes\n")
    
    allResults$Marg = resMarg
    save(allResults, file=progressFile)
    cat("Progress saved to:", progressFile, "\n")
  }, error = function(e) {
    cat("\nERROR in Marg model:\n")
    cat(conditionMessage(e), "\n")
    cat("Skipping Marg model and continuing with others\n")
  })
} else if(runMarg && 'Marg' %in% names(allResults)) {
  cat("\nSkipping Marg model (already completed)\n")
}


# ---- Run Constr model (constrained 2n-2) ----
if(runConstr && (!('Constr' %in% names(allResults)) || rerun)) {
  cat("\n############################\n")
  cat("Running M_M Constr model\n")
  cat("############################\n")
  
  tryCatch({
    outConstr = fitMMConstr(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, 
                            doMCMC=FALSE, getSDs=getSDs)
    
    resConstr = summarizeResults(outConstr$TMBsd, outConstr$TMBobj, "M_M Constr")
    cat("\nTotal time:", round(outConstr$totalTime/60, 2), "minutes\n")
    
    allResults$Constr = resConstr
    save(allResults, file=progressFile)
    cat("Progress saved to:", progressFile, "\n")
  }, error = function(e) {
    cat("\nERROR in Constr model:\n")
    cat(conditionMessage(e), "\n")
    cat("Skipping Constr model and continuing with others\n")
  })
} else if(runConstr && 'Constr' %in% names(allResults)) {
  cat("\nSkipping Constr model (already completed)\n")
}


# ---- Initialize MCMC results tracking ----
mcmcProgressFile = "savedOutput/testMM_mcmc_progress.RData"
if(rerun) {
  mcmcResults = list()
  cat("MCMC rerun mode - starting fresh\n\n")
} else if(file.exists(mcmcProgressFile)) {
  cat("Loading MCMC progress from:", mcmcProgressFile, "\n")
  load(mcmcProgressFile)
  cat("Restored MCMC results for models:", names(mcmcResults), "\n\n")
} else {
  mcmcResults = list()
  cat("Starting fresh MCMC runs (no progress file found)\n\n")
}


# ----  Run Marg MCMC model (marginalized) ----
if(doMCMC && runMarg && (!('MargMCMC' %in% names(mcmcResults)) || rerun)) {
  cat("\n############################\n")
  cat("Running M_M Marg MCMC model\n")
  cat("############################\n")
  
  tryCatch({
    outMargMCMC = fitMMMarg(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, 
                            doMCMC=TRUE, getSDs=FALSE)
    
    cat("\nTotal MCMC time:", round(outMargMCMC$totalTime/60, 2), "minutes\n")
    
    # Store the fit object
    mcmcResults$MargMCMC = outMargMCMC
    save(mcmcResults, file=mcmcProgressFile)
    cat("MCMC progress saved to:", mcmcProgressFile, "\n")
  }, error = function(e) {
    cat("\nERROR in Marg MCMC model:\n")
    cat(conditionMessage(e), "\n")
    cat("Skipping Marg MCMC model and continuing with others\n")
  })
} else if(doMCMC && runMarg && 'MargMCMC' %in% names(mcmcResults)) {
  cat("\nSkipping Marg MCMC model (already completed)\n")
} else if(!doMCMC) {
  cat("\nSkipping MCMC models (tmbstan not available)\n")
}


# ---- Run Constr MCMC model (constrained 2n-2) ----
if(doMCMC && runConstr && (!('ConstrMCMC' %in% names(mcmcResults)) || rerun)) {
  cat("\n############################\n")
  cat("Running M_M Constr MCMC model\n")
  cat("############################\n")
  
  tryCatch({
    outConstrMCMC = fitMMConstr(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, 
                               doMCMC=TRUE, getSDs=FALSE)
    
    cat("\nTotal MCMC time:", round(outConstrMCMC$totalTime/60, 2), "minutes\n")
    
    # Store the fit object
    mcmcResults$ConstrMCMC = outConstrMCMC
    save(mcmcResults, file=mcmcProgressFile)
    cat("MCMC progress saved to:", mcmcProgressFile, "\n")
  }, error = function(e) {
    cat("\nERROR in Constr MCMC model:\n")
    cat(conditionMessage(e), "\n")
    cat("Skipping Constr MCMC model and continuing with others\n")
  })
} else if(doMCMC && runConstr && 'ConstrMCMC' %in% names(mcmcResults)) {
  cat("\nSkipping Constr MCMC model (already completed)\n")
}


# ---- Run SepRepar MCMC model (original) ----
if(doMCMC && runSepRepar && (!('SepReparMCMC' %in% names(mcmcResults)) || rerun)) {
  cat("\n############################\n")
  cat("Running M_M SepRepar MCMC model\n")
  cat("############################\n")
  
  tryCatch({
    outSepReparMCMC = fitMM(datDHS=thisDHS, datMICS=thisMICS, intPtsDHS=intPtsDHS, 
                            repar=TRUE, doMCMC=TRUE, MdInit=FALSE)
    
    cat("\nTotal MCMC time:", round(outSepReparMCMC$totalTime/60, 2), "minutes\n")
    
    # Store the fit object
    mcmcResults$SepReparMCMC = outSepReparMCMC
    save(mcmcResults, file=mcmcProgressFile)
    cat("MCMC progress saved to:", mcmcProgressFile, "\n")
  }, error = function(e) {
    cat("\nERROR in SepRepar MCMC model:\n")
    cat(conditionMessage(e), "\n")
    cat("Skipping SepRepar MCMC model and continuing with others\n")
  })
} else if(doMCMC && runSepRepar && 'SepReparMCMC' %in% names(mcmcResults)) {
  cat("\nSkipping SepRepar MCMC model (already completed)\n")
}


# ---- Reload all results before comparison ----
if(file.exists(progressFile)) {
  cat("Reloading all results from:", progressFile, "\n")
  load(progressFile)
  cat("Loaded models:", names(allResults), "\n")
} else {
  cat("No progress file found to reload\n")
}


# ---- Compare results across models ----
cat("\n\n========================================\n")
cat("COMPARISON SUMMARY\n")
cat("========================================\n")

if(length(allResults) > 0) {
  cat("\nCompleted models:", names(allResults), "\n\n")
  cat("\nHyperparameter estimates:\n")
  compTab = data.frame(
    parameter = c("sigmaSq", "phi", "sigmaEpsSq"),
    truth = c(trueSigmaSq, NA, trueSigmaEpsSq)
  )
  for(nm in names(allResults)) {
    r = allResults[[nm]]
    compTab[[nm]] = round(c(r$sigmaSq, r$phi, r$sigmaEpsSq), 4)
  }
  print(compTab)
  
  cat("\nBeta estimates:\n")
  betaTab = data.frame(
    name = trueParNames[2:length(trueParNames)],
    truth = truePar[2:length(truePar)]
  )
  for(nm in names(allResults)) {
    betaTab[[nm]] = round(allResults[[nm]]$beta, 4)
  }
  print(betaTab)
  
  # ---- Final comprehensive report ----
  reportFile = "savedOutput/testMM_report.txt"
  sink(reportFile)
  
  cat("========================================\n")
  cat("M_M MODEL COMPARISON REPORT\n")
  cat("========================================\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Simulation:", i, "\n\n")
  
  cat("True parameters:\n")
  print(data.frame(name=trueParNames, value=truePar))
  cat("True sigmaSq:", trueSigmaSq, "\n")
  cat("True sigmaEpsSq:", trueSigmaEpsSq, "\n\n")
  
  cat("Completed models:", paste(names(allResults), collapse=", "), "\n\n")
  
  cat("Hyperparameter estimates:\n")
  print(compTab)
  
  cat("\n\nBeta estimates vs truth:\n")
  print(betaTab)
  
  cat("\n\nDetailed results by model:\n")
  for(nm in names(allResults)) {
    cat("\n---", nm, "---\n")
    print(allResults[[nm]])
  }
  
  sink()
  cat("\nReport saved to:", reportFile, "\n")
  
} else {
  cat("No results to compare (no models completed)\n")
}

cat("\n========================================\n")
cat("TEST COMPLETED\n")
cat("========================================\n")
cat("Script finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")


# ---- Reload and visualize MCMC results ----
if(file.exists(mcmcProgressFile)) {
  cat("\n\nReloading MCMC results from:", mcmcProgressFile, "\n")
  load(mcmcProgressFile)
  cat("MCMC models available:", paste(names(mcmcResults), collapse=", "), "\n\n")
  
  # Check if bayesplot is available for visualization
  if(!require("bayesplot", quietly=TRUE)) {
    cat("Warning: bayesplot package not available - MCMC visualizations will be skipped\n")
  } else {
    # Create output directory
    figDir = "Figures/testMCMC"
    if(!dir.exists(figDir)) {
      dir.create(figDir, recursive=TRUE)
      cat("Created directory:", figDir, "\n\n")
    }
    
    # Generate pairplots for each MCMC result
    cat("Generating MCMC pairplots...\n\n")
    
    # Helper function to extract available parameters from fit object
    extract_available_pars = function(fit_obj, requested_pars) {
      tryCatch({
        available = fit_obj$model@par_dims
        available_names = names(available)
        # Check which requested parameters exist
        found_pars = intersect(requested_pars, available_names)
        if(length(found_pars) == 0) {
          # Try with [1] notation
          found_pars = c()
          for(par in requested_pars) {
            if(grepl("\\[", par)) {
              base_name = sub("\\[[0-9]+\\]", "", par)
              if(base_name %in% available_names) {
                found_pars = c(found_pars, par)
              }
            } else if(par %in% available_names) {
              found_pars = c(found_pars, par)
            }
          }
        }
        return(found_pars)
      }, error = function(e) {
        cat("Warning: Could not extract parameters:", conditionMessage(e), "\n")
        return(requested_pars)  # Return requested as fallback
      })
    }
    
    # ---- MargMCMC pairplot ----
    if('MargMCMC' %in% names(mcmcResults)) {
      cat("Creating pairplot for Marg MCMC...\n")
      tryCatch({
        fit_marg = mcmcResults$MargMCMC$TMBstan_fit
        
        # Parameters for Marginalized model
        pars_marg_requested = c("alpha", "beta[5]", "logit_phi", "log_tau", "log_tauEps", "lp__", "Epsilon_bym2[1]")
        pars_marg = extract_available_pars(fit_marg, pars_marg_requested)
        
        if(length(pars_marg) > 0) {
          png(file.path(figDir, "pairplot_Marg_MCMC.png"), width=1200, height=1000, res=300)
          print(mcmc_pairs(fit_marg, pars=pars_marg, off_diag_args=list(size=1)))
          dev.off()
          cat("Saved: pairplot_Marg_MCMC.png (parameters:", paste(pars_marg, collapse=", "), ")\n")
        } else {
          cat("Warning: No requested parameters found in Marg MCMC fit\n")
        }
      }, error = function(e) {
        cat("ERROR creating Marg MCMC pairplot:", conditionMessage(e), "\n")
      })
    }
    
    # ---- ConstrMCMC pairplot ----
    if('ConstrMCMC' %in% names(mcmcResults)) {
      cat("Creating pairplot for Constr MCMC...\n")
      tryCatch({
        fit_constr = mcmcResults$ConstrMCMC$TMBstan_fit
        
        # Parameters for Constrained model
        pars_constr_requested = c("alpha", "beta[5]", "logit_phi", "log_tau", "log_tauEps", "lp__", "u_bym2[1]", "w_bym2[1]")
        pars_constr = extract_available_pars(fit_constr, pars_constr_requested)
        
        if(length(pars_constr) > 0) {
          png(file.path(figDir, "pairplot_Constr_MCMC.png"), width=1200, height=1000, res=300)
          print(mcmc_pairs(fit_constr, pars=pars_constr, off_diag_args=list(size=1)))
          dev.off()
          cat("Saved: pairplot_Constr_MCMC.png (parameters:", paste(pars_constr, collapse=", "), ")\n")
        } else {
          cat("Warning: No requested parameters found in Constr MCMC fit\n")
        }
      }, error = function(e) {
        cat("ERROR creating Constr MCMC pairplot:", conditionMessage(e), "\n")
      })
    }
    
    # ---- SepReparMCMC pairplot ----
    if('SepReparMCMC' %in% names(mcmcResults)) {
      cat("Creating pairplot for SepRepar MCMC...\n")
      tryCatch({
        fit_seprepar = mcmcResults$SepReparMCMC$TMBstan_fit
        
        # Parameters for SepRepar model
        pars_seprepar_requested = c("alpha", "beta[5]", "logit_phi", "log_tau", "log_tauEps", "lp__", "w[1]", "u[1]")
        pars_seprepar = extract_available_pars(fit_seprepar, pars_seprepar_requested)
        
        if(length(pars_seprepar) > 0) {
          png(file.path(figDir, "pairplot_SepRepar_MCMC.png"), width=1200, height=1000, res=300)
          print(mcmc_pairs(fit_seprepar, pars=pars_seprepar, off_diag_args=list(size=1)))
          dev.off()
          cat("Saved: pairplot_SepRepar_MCMC.png (parameters:", paste(pars_seprepar, collapse=", "), ")\n")
        } else {
          cat("Warning: No requested parameters found in SepRepar MCMC fit\n")
        }
      }, error = function(e) {
        cat("ERROR creating SepRepar MCMC pairplot:", conditionMessage(e), "\n")
      })
    }
    
    cat("\nAll MCMC pairplots saved to:", normalizePath(figDir), "\n")
  }
} else {
  cat("\nNo MCMC results to visualize\n")
}

# ---- Final execution summary ----
executionLog$endTime = Sys.time()
executionLog$totalTime = difftime(executionLog$endTime, executionLog$startTime, units="mins")
executionLog$optimizationCount = length(allResults)
executionLog$mcmcCount = length(mcmcResults)

cat("\n\n========================================\n")
cat("EXECUTION SUMMARY\n")
cat("========================================\n")
cat("Start time:", format(executionLog$startTime, "%Y-%m-%d %H:%M:%S"), "\n")
cat("End time:", format(executionLog$endTime, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total duration:", round(executionLog$totalTime, 2), "minutes\n")
cat("Optimization models completed:", executionLog$optimizationCount, "/3\n")
cat("MCMC models completed:", executionLog$mcmcCount, "/3\n")
cat("MCMC available:", executionLog$mcmcStatus$available, "\n")
cat("\n✓ Script completed successfully\n")
cat("========================================\n")
