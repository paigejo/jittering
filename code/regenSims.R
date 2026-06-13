source("code/setup.R")

# simData1BYM2 now builds bym2ArgsTMB internally by default (same constructor
# the fits use), so no need to pass it explicitly here.
cat(sprintf("\n--- %s | BYM2 regenerate ---\n", format(Sys.time())))
simulateSurveys("bym2", nsim = 100, seed = 123, regenerate = TRUE)
cat(sprintf("\n--- %s | SPDE regenerate ---\n", format(Sys.time())))
simulateSurveys("spde", nsim = 100, seed = 123, regenerate = TRUE)
cat(sprintf("\n--- %s | All done ---\n", format(Sys.time())))
