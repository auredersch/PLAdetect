# validation: PLAenriched vs. benchmark
# vaccine | GOBP | AUCell | gmm_dist_dual | global | raw | GT biologist

library(Seurat)
devtools::load_all("/nfs/home/students/a.dersch/PLAenriched")

project_root <- "/nfs/home/students/a.dersch/FoPra_PLAs"
BENCH_FILE   <- "metrics_vaccine_AUCell_GOBP_gmm_dist_dual_Scopeglobal_raw_GT_biologist.csv"   

pbmc <- readRDS(file.path(project_root, "data/datasets_automated/vaccine_raw_automated_gating.rds"))
DefaultAssay(pbmc) <- "RNA"

# --- PLAenriched ---
pbmc <- detect_pla(pbmc, platelet_signature = "GOBP", method = "AUCell")

# --- CONFUSION MATRIX AGAINST THE BIOLOGIST ---
GT_COLUMN <- if ("pla_status" %in% colnames(pbmc@meta.data)) "pla_status" else "pla.status"
gt_values <- as.character(pbmc@meta.data[[GT_COLUMN]])
gt_values[is.na(gt_values)] <- "platelet-free"

pred_pla <- pbmc$PLA_status == "PLA"
tp <- sum(pred_pla & gt_values == "PLA")
fp <- sum(pred_pla & gt_values != "PLA")
fn <- sum(!pred_pla & gt_values == "PLA")
tn <- sum(!pred_pla & gt_values != "PLA")

prec <- if ((tp + fp) > 0) tp / (tp + fp) else 0
rec  <- if ((tp + fn) > 0) tp / (tp + fn) else 0
f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0

# --- COMPARISON ---
bench <- read.csv(file.path(project_root, "results/benchmarking/metrics", BENCH_FILE))
vergleich <- rbind(
  PLAenriched = data.frame(TP = tp, FP = fp, FN = fn, TN = tn, F1 = round(f1, 4)),
  Benchmark = data.frame(TP = bench$TP, FP = bench$FP, FN = bench$FN, TN = bench$TN,
                         F1 = round(bench$F1, 4))
)
print(vergleich)
print(all(vergleich["PLAenriched", ] == vergleich["Benchmark", ]))