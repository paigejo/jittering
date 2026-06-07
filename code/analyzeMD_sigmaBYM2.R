# Pull sigmaBYM2 estimates for M_D_BYM2 across all available sims, on both
# cluster (100 sims) and laptop (10 sims). This is independent of predGrid
# inputs — pure fit-side parameter recovery — so any disagreement is a
# property of the FIT, not the prediction.
#
# Reports:
#   - Per-sim {est, sd, truth, bias} table on cluster
#   - Same on laptop
#   - Distribution summary (quantiles) for each side
#   - Side-by-side overlap on sims 1-10
#
# Also does the same for M_DM_BYM2 for comparison — we expect M_DM_BYM2 to
# behave well (BYM2 truth = BYM2 model), while M_D_BYM2 is "MICS-only BYM2"
# which is known to be data-poor and may genuinely have poor recovery.

clustBase <- "c:/Users/jpaige/OneDrive - Norsk Regnesentral/Projects/jittering/jitterScores/scores_full_cluster"
localBase <- "c:/Users/jpaige/git/jittering/savedOutput/simStudy1/scores"

extractHyper <- function(baseDir, tag, model, nMax) {
    rows <- list()
    for(i in seq_len(nMax)) {
        p <- sprintf("%s/%s/scores_%s_sim%d.RData", baseDir, tag, model, i)
        if(!file.exists(p)) next
        e <- new.env(); load(p, envir = e)
        if(!exists("scoresHyper", envir = e) || is.null(e$scoresHyper)) next
        sh <- e$scoresHyper
        if(!"sigmaBYM2" %in% rownames(sh)) next
        rows[[length(rows) + 1]] <- data.frame(
            sim   = i,
            est   = as.numeric(sh["sigmaBYM2", "est"]),
            sd    = as.numeric(sh["sigmaBYM2", "sd"]),
            truth = as.numeric(sh["sigmaBYM2", "truth"]),
            bias  = as.numeric(sh["sigmaBYM2", "est"]) - as.numeric(sh["sigmaBYM2", "truth"]),
            stringsAsFactors = FALSE
        )
    }
    if(length(rows) == 0) return(NULL)
    do.call(rbind, rows)
}

summarizeDist <- function(d, label) {
    cat(sprintf("\n[%s]  n=%d\n", label, nrow(d)))
    cat(sprintf("  est:   median=%.4f  mean=%.4f  min=%.4f  max=%.4f  IQR=[%.4f, %.4f]\n",
                median(d$est),  mean(d$est),  min(d$est),  max(d$est),
                quantile(d$est, 0.25), quantile(d$est, 0.75)))
    cat(sprintf("  truth: median=%.4f  mean=%.4f  (constant across sims if DGP fixed)\n",
                median(d$truth), mean(d$truth)))
    cat(sprintf("  bias:  median=%.4f  mean=%.4f  IQR=[%.4f, %.4f]\n",
                median(d$bias), mean(d$bias),
                quantile(d$bias, 0.25), quantile(d$bias, 0.75)))
    # How many sims have outrageous estimates (>5 or >10)?
    cat(sprintf("  # sims with est > 2: %d   > 5: %d   > 10: %d   > 50: %d\n",
                sum(d$est > 2), sum(d$est > 5), sum(d$est > 10), sum(d$est > 50)))
}

for(tag in c("BYM2", "SPDE")) {
    for(model in c("M_D_BYM2", "M_DM_BYM2")) {
        cat(sprintf("\n=========== generative=%s | model=%s ===========\n", tag, model))
        clust <- extractHyper(clustBase, tag, model, 100)
        local <- extractHyper(localBase, tag, model, 10)

        if(is.null(clust)) cat("(no cluster files)\n")
        else { summarizeDist(clust, "cluster (100 sims)") }
        if(is.null(local)) cat("(no local files)\n")
        else { summarizeDist(local, "local (10 sims)") }

        # Side-by-side for sims 1-10
        if(!is.null(clust) && !is.null(local)) {
            cat("\n  Side-by-side sims 1-10:\n")
            m <- merge(clust[clust$sim <= 10, c("sim","est","sd","bias")],
                       local[, c("sim","est","sd","bias")],
                       by = "sim", suffixes = c("_clust","_local"))
            print(round(m, 4), row.names = FALSE)
        }

        # Always show the cluster's worst 5 sims (helps identify outliers)
        if(!is.null(clust)) {
            cat("\n  Cluster worst 5 (largest |est - truth|):\n")
            worst <- clust[order(-abs(clust$bias)), ][1:min(5, nrow(clust)), ]
            print(round(worst, 4), row.names = FALSE)
        }
    }
}
