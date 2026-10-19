# Built-in signatures shipped in inst/extdata/signatures/.
.builtin_signatures <- function() {
  c(
    # curated lists
    REACTOME     = "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION.v2025.1.Hs.csv",
    GOBP         = "GOBP_REGULATION_OF_PLATELET_ACTIVATION.v2025.1.Hs.csv",
    MANNE_DN     = "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs.csv",
    MANNE_UP     = "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_UP.v2025.1.Hs.csv",
    GNATENKO     = "GNATENKO_PLATELET_SIGNATURE.v2025.1.Hs.csv",
    HP           = "HP_ABNORMAL_PLATELET_MEMBRANE_PROTEIN_EXPRESSION.v2025.1.Hs.csv",
    WP           = "WP_PLATELETMEDIATED_INTERACTIONS_WITH_VASCULAR_AND_CIRCULATING_CELLS.v2025.1.Hs.csv",
    OVERLAP      = "OVERLAP_GOBP_MANNEDN_REACTOME.csv",
    OVERLAP_MIN2 = "OVERLAP_MIN2.csv",
    UNION_ALL    = "UNION_ALL.csv",
    # derived from cell-cell communication analyses (Sepsis, Immune Aging, Skin)
    LIANAPLUS     = "LIANAPLUS_top20_sepsis_immune_skin.csv",
    MULTINICHENET = "MULTINICHENET_top20_sepsis_immune_skin.csv",
    SCDIFFCOM     = "scDIFFCOM_top20_sepsis_immune_skin.csv",
    # leukocyte activation list, default for the immune gate
    GOBP_LEUKOCYTE = "GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs.csv"
  )
}

#' Names of the built-in gene signatures
#' @return Character vector of names usable as `platelet_signature` / `immune_signature`.
#' @export
list_signatures <- function() names(.builtin_signatures())

#' Read a gene signature
#'
#' @param x One of: a character vector of gene symbols (length > 1), a path to a
#'   `.csv` file (column `geneName`) or `.gmt` file, or the name of a built-in
#'   signature (see [list_signatures()]).
#' @param set For GMT files with several gene sets: name of the set to use.
#' @return Character vector of gene symbols.
#' @export
read_signature <- function(x, set = NULL) {
  if (!is.character(x) || length(x) < 1) {
    stop("`x` must be a gene vector, a file path, or a built-in signature name.")
  }
  if (length(x) > 1) return(.clean_genes(x))

  path <- NULL
  if (file.exists(x)) {
    path <- x
  } else if (x %in% names(.builtin_signatures())) {
    path <- system.file("extdata", "signatures", .builtin_signatures()[[x]],
                        package = "PLAdetect")
    if (!nzchar(path)) {
      stop("Built-in signature '", x, "' is registered but its file is missing in the package.")
    }
  } else {
    stop("'", x, "' is neither an existing file nor a built-in signature. ",
         "Built-in: ", paste(list_signatures(), collapse = ", "), ".")
  }

  if (grepl("\\.gmt$", path, ignore.case = TRUE)) {
    genes <- .read_gmt(path, set)
  } else {
    df <- utils::read.csv(path, stringsAsFactors = FALSE)
    if ("geneName" %in% colnames(df)) {
      genes <- df$geneName
    } else if (ncol(df) == 1) {
      genes <- df[[1]]
    } else {
      stop("CSV needs a column named 'geneName' (found: ",
           paste(colnames(df), collapse = ", "), ").")
    }
  }
  .clean_genes(genes)
}

.clean_genes <- function(genes) {
  genes <- trimws(as.character(genes))
  unique(genes[!is.na(genes) & nzchar(genes)])
}

.read_gmt <- function(path, set = NULL) {
  parts <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  parts <- parts[lengths(parts) >= 3]
  if (length(parts) == 0) stop("No gene sets found in GMT file: ", path)
  names(parts) <- vapply(parts, `[`, character(1), 1)
  if (length(parts) > 1) {
    if (is.null(set)) {
      stop("GMT file has several sets (", paste(names(parts), collapse = ", "),
           "). Choose one with `set =`.")
    }
    if (!set %in% names(parts)) stop("Set '", set, "' not found in GMT file.")
    parts <- parts[set]
  }
  parts[[1]][-(1:2)]
}

.check_overlap <- function(genes, available, label, min_genes = 3, min_overlap = 0.5) {
  overlap <- intersect(genes, available)
  frac <- length(overlap) / length(genes)
  message(sprintf("%s signature: %d of %d genes found in the data (%.1f%%).",
                  label, length(overlap), length(genes), 100 * frac))
  if (length(overlap) < min_genes) {
    stop("Only ", length(overlap), " ", label, " genes found (minimum ", min_genes, "). ",
         "Data genes e.g.: ", paste(utils::head(available, 3), collapse = ", "),
         "; signature genes e.g.: ", paste(utils::head(genes, 3), collapse = ", "), ".")
  }
  if (frac < min_overlap) {
    stop("Only ", round(100 * frac, 1), "% of the ", label, " signature found in the data. ",
         "This usually means a gene identifier mismatch (symbols vs. Ensembl IDs).")
  }
  if (length(overlap) < 10) {
    warning("Only ", length(overlap), " ", label, " genes found; the score will be noisy.")
  }
  overlap
}