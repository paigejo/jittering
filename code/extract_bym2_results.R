load("savedOutput/global/BYM2_models_test_results.RData")

for(nm in names(out)) {
  res <- out[[nm]]
  nll <- res$opt$objective
  conv <- res$opt$convergence
  s <- summary(res$TMBsd, "fixed")
  fe_rows <- !grepl("^w_bym2|^u_bym2", rownames(s))
  fe <- s[fe_rows, , drop=FALSE]
  cat(sprintf("\n== %s | CONV: %d | NLL: %.4f ==\n", nm, conv, nll))
  print(round(fe, 4))
}
