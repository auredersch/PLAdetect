# Scoring of the platelet and the immune signature
# signatures: named list, e.g. list(platelet = genes, immune = immune_genes)
# returns a data.frame with one column per signature, rows = cells
score_cells <- function(seu, signatures,
                        method = c("AUCell", "UCell", "AddModuleScore"),
                        assay = "RNA", layer = "data") {
  method <- match.arg(method)

  # --- SCORING ---
  if (method == "AUCell") {
    expression_matrix <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
    rankings <- AUCell::AUCell_buildRankings(expression_matrix, plotStats = FALSE, verbose = FALSE)
    auc <- AUCell::AUCell_calcAUC(signatures, rankings, verbose = FALSE)
    scores <- t(AUCell::getAUC(auc))

  } else if (method == "UCell") {
    if (!requireNamespace("UCell", quietly = TRUE)) {
      stop("method = 'UCell' needs the package UCell: BiocManager::install('UCell')")
    }
    expression_matrix <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
    scores <- UCell::ScoreSignatures_UCell(expression_matrix, features = signatures, name = "")
    colnames(scores) <- sub("_UCell$", "", colnames(scores))  

  } else {
    # AddModuleScore appends one column per feature set: PLAenrichedAMS1, PLAenrichedAMS2
    tmp <- Seurat::AddModuleScore(seu, features = unname(signatures), assay = assay,
                                  name = "PLAenrichedAMS")
    scores <- as.matrix(tmp@meta.data[, paste0("PLAenrichedAMS", seq_along(signatures)), drop = FALSE])
    colnames(scores) <- names(signatures)
  }

  scores <- as.data.frame(scores[colnames(seu), names(signatures), drop = FALSE])
  if (anyNA(scores)) stop("Scoring produced NA values (method = ", method, ").")
  scores
}