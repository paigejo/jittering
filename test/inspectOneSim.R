# Drill into a single (sim, model) pair on both machines.
clustBase <- "c:/Users/jpaige/OneDrive - Norsk Regnesentral/Projects/jittering/jitterScores/scores_full_cluster"
localBase <- "c:/Users/jpaige/git/jittering/savedOutput/simStudy1/scores"

pickPair <- function(tag, model, sim) {
    cPath <- sprintf("%s/%s/scores_%s_sim%d.RData", clustBase, tag, model, sim)
    lPath <- sprintf("%s/%s/scores_%s_sim%d.RData", localBase, tag, model, sim)
    cE <- new.env(); load(cPath, envir = cE)
    lE <- new.env(); load(lPath, envir = lE)
    cat(sprintf("\n#### %s | %s | sim%d ####\n", tag, model, sim))
    cat("cluster: ", paste(ls(cE), collapse=", "), "\n")
    cat("local:   ", paste(ls(lE), collapse=", "), "\n")

    for(field in c("scoresFE", "scoresHyper", "scoresArea")) {
        cat(sprintf("\n--- %s ---\n", field))
        cHas <- exists(field, envir=cE); lHas <- exists(field, envir=lE)
        if(!cHas) { cat("(missing in cluster)\n") }
        if(!lHas) { cat("(missing in local)\n") }
        if(cHas && lHas) {
            cM <- as.matrix(get(field, envir=cE))
            lM <- as.matrix(get(field, envir=lE))
            cat("cluster dim: ", paste(dim(cM), collapse=" x "),
                "  local dim: ", paste(dim(lM), collapse=" x "), "\n")
            if(all(dim(cM) == dim(lM))) {
                colMax <- apply(abs(cM - lM), 2, max, na.rm = TRUE)
                cat("max|cluster - local| per column:\n")
                print(round(colMax, 4))
                cat("\ncluster (head):\n");  print(round(cM, 4))
                cat("\nlocal (head):\n");    print(round(lM, 4))
            } else {
                cat("DIM MISMATCH\n")
                cat("cluster colnames: ", paste(colnames(cM), collapse=","), "\n")
                cat("local   colnames: ", paste(colnames(lM), collapse=","), "\n")
            }
        }
    }
}

pickPair("BYM2", "M_DM_BYM2", 1)
pickPair("BYM2", "M_D_BYM2",  1)
pickPair("BYM2", "Md_FE",     1)
