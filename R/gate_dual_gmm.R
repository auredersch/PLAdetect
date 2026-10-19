# fits a 2-component GMM, returns NULL if mclust fails
fit_gmm_pair <- function(values) {
  tryCatch(mclust::Mclust(values, G = 2, verbose = FALSE), error = function(e) NULL)
}

# global dual-gate GMM; gate 1: platelet score of all cells, gate 2: immune score
gate_dual_gmm <- function(raw_platelet, raw_immune, min_cells = 100, min_gate2_cells = 20) {
  n_cells <- length(raw_platelet)
  if (n_cells < min_cells) {
    stop("Only ", n_cells, " cells, the GMM fit is not reliable below ", min_cells, " cells.")
  }

  # --- GATE 1: PLATELET ---
  fit_plat <- fit_gmm_pair(raw_platelet)
  if (is.null(fit_plat)) stop("Mclust (Platelet) failed!")
  plat_high <- which.max(fit_plat$parameters$mean)
  Platelet_High <- fit_plat$classification == plat_high
  Platelet_Prob <- fit_plat$z[, plat_high]

  # --- GATE 2: IMMUNE (only platelet-high cells) ---
  idx <- which(Platelet_High)
  if (length(idx) < min_gate2_cells) {
    stop("Only ", length(idx), " platelet-high cells after gate 1, cannot fit gate 2.")
  }
  fit_imm <- fit_gmm_pair(raw_immune[idx])
  if (is.null(fit_imm)) stop("Mclust (Immune) failed!")
  imm_high <- which.max(fit_imm$parameters$mean)

  Immune_High <- rep(FALSE, n_cells)
  Immune_Prob <- rep(NA_real_, n_cells)    # NA for cells that did not pass gate 1
  Immune_High[idx] <- fit_imm$classification == imm_high
  Immune_Prob[idx] <- fit_imm$z[, imm_high]

  list(
    PLA = Platelet_High & Immune_High,
    Platelet_High = Platelet_High, Immune_High = Immune_High,
    Platelet_Prob = Platelet_Prob, Immune_Prob = Immune_Prob,
    Threshold_Platelet = min(raw_platelet[Platelet_High]),
    Threshold_Immune = if (any(Immune_High)) min(raw_immune[Immune_High]) else NA_real_
  )
}