#' Detect platelet-leukocyte aggregates (PLAs) in scRNA-seq data
#'
#' Scores every cell for a platelet signature and for a leukocyte activation
#' signature and calls PLAs with a global dual-gate GMM.
#'
#' @param seu Seurat (v5) object with log-normalized RNA data and gene symbols.
#' @param platelet_signature Built-in name (see [list_signatures()]), path to a
#'   .csv/.gmt file, or a character vector of gene symbols.
#' @param immune_signature Same options as `platelet_signature`; used for the
#'   second gate. Default is the built-in `GOBP_LEUKOCYTE`.
#' @param method `"AUCell"` (default, benchmarked), `"UCell"` or
#'   `"AddModuleScore"` (not benchmarked).
#' @param assay,layer Assay and layer with the log-normalized data.
#' @param min_genes,min_overlap Stop if fewer genes (absolute) or a smaller
#'   fraction of a signature is found in the data.
#' @param max_plausible_frac Warn if the share of PLA-labeled cells is higher.
#' @param seed Seed set at the start of the run (mclust starts from a random
#'   subset in large datasets). `NULL` leaves the random number generator untouched.
#' @return The Seurat object with new meta data columns `PLA_status`,
#'   `PLA_score`, `PLA_score_z`, `PLA_immune_score`, `PLA_platelet_prob`,
#'   `PLA_immune_prob`, `PLA_platelet_high` and `PLA_immune_high`. Thresholds
#'   and settings are stored in `seu@misc$PLAdetect`.
#' @export
detect_pla <- function(seu,
                       platelet_signature = "REACTOME",
                       immune_signature = "GOBP_LEUKOCYTE",
                       method = c("AUCell", "UCell", "AddModuleScore"),
                       assay = "RNA", layer = "data",
                       min_genes = 3, min_overlap = 0.5,
                       max_plausible_frac = 0.30,
                       seed = 42) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  # --- INPUT CHECKS ---
  if (!inherits(seu, "Seurat")) stop("`seu` must be a Seurat object.")
  if (!assay %in% SeuratObject::Assays(seu)) stop("Assay '", assay, "' not found.")
  available_layers <- SeuratObject::Layers(seu[[assay]])
  if (!layer %in% available_layers) {
    stop("Layer '", layer, "' not found in assay '", assay, "' (available: ",
         paste(available_layers, collapse = ", "), "). Run NormalizeData(), ",
         "or JoinLayers() if the object is split.")
  }
  expression_matrix <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)

  # integer values only -> probably raw counts, not log-normalized
  if (methods::is(expression_matrix, "dgCMatrix")) {
    test_values <- expression_matrix@x[seq_len(min(1e5, length(expression_matrix@x)))]
    if (all(test_values == round(test_values))) {
      warning("Layer '", layer, "' only contains integers, it looks like raw counts ",
              "and not log-normalized data.")
    }
  }

  # --- SIGNATURES ---
  genes <- .check_overlap(read_signature(platelet_signature), rownames(expression_matrix),
                          "platelet", min_genes, min_overlap)
  immune_genes <- .check_overlap(read_signature(immune_signature), rownames(expression_matrix),
                                 "immune", min_genes, min_overlap)

  # --- SCORING ---
  message("--- Calculating scores using ", method, " ---")
  scores <- score_cells(seu, list(platelet = genes, immune = immune_genes),
                        method = method, assay = assay, layer = layer)
  if (sd(scores$platelet) == 0 || sd(scores$immune) == 0) {
    stop("A score has zero variance, the GMM cannot be fitted.")
  }

  # --- GATING ---
  gate <- gate_dual_gmm(scores$platelet, scores$immune)

  # --- OUTPUT ---
  meta <- data.frame(
    PLA_status        = ifelse(gate$PLA, "PLA", "platelet-free"),
    PLA_score         = scores$platelet,
    PLA_score_z       = as.vector(scale(scores$platelet)),
    PLA_immune_score  = scores$immune,
    PLA_platelet_prob = gate$Platelet_Prob,
    PLA_immune_prob   = gate$Immune_Prob,
    PLA_platelet_high = gate$Platelet_High,
    PLA_immune_high   = gate$Immune_High,
    row.names = colnames(seu)
  )
  seu <- SeuratObject::AddMetaData(seu, meta)

  pla_frac <- mean(gate$PLA)
  message(sprintf("%d of %d cells labeled as PLA (%.1f%%).",
                  sum(gate$PLA), length(gate$PLA), 100 * pla_frac))
  if (!is.null(max_plausible_frac) && pla_frac > max_plausible_frac) {
    warning(sprintf("%.1f%% of the cells are labeled as PLA, more than the plausible maximum of %.0f%%. ",
                    100 * pla_frac, 100 * max_plausible_frac),
            "The method fails mostly by labeling too many cells, check signature and score distribution.")
  }

  seu@misc[["PLAdetect"]] <- list(
    method = method, assay = assay, layer = layer,
    n_platelet_genes = length(genes), n_immune_genes = length(immune_genes),
    threshold_platelet = gate$Threshold_Platelet,
    threshold_immune = gate$Threshold_Immune,
    version = as.character(utils::packageVersion("PLAdetect")),
    seed = seed
  )
  seu
}