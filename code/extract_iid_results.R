load("savedOutput/global/IID_models_test_results.RData")

for(nm in names(out)) {
  res <- out[[nm]]
  nll <- res$opt$value
  conv <- res$opt$convergence
  s <- summary(res$TMBsd, "fixed")
  fe_rows <- !grepl("^v_iid", rownames(s))
  fe <- s[fe_rows, , drop=FALSE]
  cat(sprintf("\n== %s | CONV: %d | NLL: %.4f ==\n", nm, conv, nll))
  print(round(fe, 4))
}
