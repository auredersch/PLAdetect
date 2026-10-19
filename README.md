
# PLAenriched

RNA-based detection of platelet-leukocyte aggregates (PLAs) in scRNA-seq data.

PLAs are leukocytes with adherent platelets. In single-cell data they show up as cells that carry a platelet RNA signature on top of their own lineage. PLAenriched scores every cell for a platelet signature and for a leukocyte activation signature and calls PLAs with a global dual-gate Gaussian mixture model. It needs no ADT data, no manual gating and no cell type annotation.

The package was called PLAdetect in the FoPra report. The method and its benchmark across six datasets (Vaccine, Sepsis, Heart, Skin, Immune Aging, Impact) are described in that report "Beyond Enrichment: Computational Detection and Characterization of Platelet-Leukocyte Aggregates in Single-Cell Data" (TUM, List Lab).


## Installation

You can install the development version of PLAenriched from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("auredersch/PLAenriched")
```

## Example

The input is a Seurat v5 object with log-normalized RNA data and gene symbols as feature names.

``` r
library(PLAenriched)
library(Seurat)

# built-in gene lists
list_signatures()

# run with defaults: AUCell, REACTOME + GOBP_LEUKOCYTE
seu <- detect_pla(seu)

table(seu$PLA_status)
DimPlot(seu, group.by = "PLA_status")
FeaturePlot(seu, features = "PLA_score")
```

A custom platelet signature can be given as a gene vector or as a path to a .csv or .gmt file:

```r
seu <- detect_pla(seu, platelet_signature = c("PF4", "PPBP", "GP9", "ITGA2B", "TUBB1"))
```

## How it works

1. Each cell is scored for the platelet signature and the immune activation signature -> AUCell, top 5% of the gene ranking
2. Gate 1: a two-component GMM on the platelet score of all cells separates platelet-high from platelet-low cells
3. Gate 2: a two-component GMM on the immune activation score of the platelet-high cells only separates immune-high from immune-low cells
4. A cell is called PLA if it passes both gates

The thresholds are fitted once per dataset (global scope). Per-cell-type thresholds were tested during the benchmark and performed worse.

## Output

detect_pla()` returns the Seurat object with these metadata columns:

| Column | Content |
|---|---|
| `PLA_status` | PLA label per cell |
| `PLA_score` | continuous platelet score |
| `PLA_score_z` | z-scored platelet score |
| `PLA_immune_score` | immune activation score |
| `PLA_platelet_prob` | posterior of the platelet gate |
| `PLA_immune_prob` | posterior of the immune gate |
| `PLA_platelet_high` | passed gate 1 |
| `PLA_immune_high` | passed gate 2 |

The settings of the run (method, assay, layer, signature gene counts, thresholds, package version, seed) are stored in `seu@misc$PLAenriched`.

## Arguments

| Argument | Default | Meaning |
|---|---|---|
| `platelet_signature` | `"REACTOME"` | built-in name, `.csv`/`.gmt` path or gene vector |
| `immune_signature` | `"GOBP_LEUKOCYTE"` | signature for the second gate |
| `method` | `"AUCell"` | `"AUCell"` (benchmarked), `"UCell"`, `"AddModuleScore"` (not benchmarked) |
| `assay`, `layer` | `"RNA"`, `"data"` | location of the log-normalized data |
| `min_genes` | `3` | minimum number of signature genes found in the data |
| `min_overlap` | `0.5` | minimum fraction of the signature found in the data |
| `max_plausible_frac` | `0.30` | PLA share above which a warning is issued |
| `seed` | `42` | random seed, `NULL` leaves the RNG untouched |

## Checks

- The function stops if the assay or layer is missing, if a signature has fewer than `min_genes` genes or less than `min_overlap` of its genes in the data (usually a gene identifier mismatch), or if a score has no variance.
- It warns if the data look like raw counts. Run `NormalizeData()` first.
- It reports the share of PLA-labeled cells and warns if it exceeds `max_plausible_frac`.

## Limitations

- The GMM does not know the true PLA prevalence. In cohorts with weak signal it can label a large share of cells as PLA. Always check the reported share and the score distribution
- The platelet signature score is influenced by lineage identity and by platelet RNA that is not bound to a cell. A PLA label is a hypothesis that should be confirmed, for example with ADT 
- Defaults are based on the benchmark in the FoPra report. UCell and AddModuleScore were not benchmarked

## Citation

Dersch A, Kaciran I, Mathis F, Raka L (supervisor). Beyond Enrichment: Computational Detection and Characterization of Platelet-Leukocyte Aggregates in Single-Cell Data. FoPra report, Technical University of Munich, List Lab, 2026.

