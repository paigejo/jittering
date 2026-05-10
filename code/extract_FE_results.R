#!/usr/bin/env Rscript
# Extract and print FE model test results in a readable format

load("savedOutput/global/FE_models_test_results.RData")

cat("\n============================================================\n")
cat("FE MODEL TEST RESULTS SUMMARY\n")
cat("============================================================\n\n")

for (model_name in c("fitFEM", "fitFED", "fitFEMD")) {
  cat("\n", strrep("=", 60), "\n")
  cat("MODEL:", model_name, "\n")
  cat(strrep("=", 60), "\n")
  
  res <- out[[model_name]]
  
  cat("Convergence:", res$opt$convergence, "\n")
  cat("Negative Log-Likelihood (NLL):", sprintf("%.4f", res$opt$objective), "\n")
  cat("Total Time (seconds):", sprintf("%.2f", res$totalTime), "\n")
  cat("SD Report Time (seconds):", sprintf("%.2f", res$sdTime), "\n")
  
  if (!is.null(res$TMBsd)) {
    cat("\nFixed Effects Summary:\n")
    cat(strrep("-", 60), "\n")
    print(summary(res$TMBsd, "fixed"))
  }
  cat("\n")
}

cat("\n============================================================\n")
cat("END OF RESULTS\n")
cat("============================================================\n\n")

# Also save as CSV for easy importing
results_df <- data.frame()

for (model_name in c("fitFEM", "fitFED", "fitFEMD")) {
  res <- out[[model_name]]
  if (!is.null(res$TMBsd)) {
    summary_table <- summary(res$TMBsd, "fixed")
    if (is.matrix(summary_table)) {
      df_temp <- data.frame(
        Model = model_name,
        Parameter = rownames(summary_table),
        Estimate = summary_table[, "Estimate"],
        Std.Error = summary_table[, "Std. Error"]
      )
      results_df <- rbind(results_df, df_temp)
    }
  }
}

write.csv(results_df, "savedOutput/global/FE_models_results_summary.csv", row.names = FALSE)
cat("Results also saved to: savedOutput/global/FE_models_results_summary.csv\n")
