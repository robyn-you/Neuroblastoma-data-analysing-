## 01b_import_qc.R
## Import, alignment and missingness QC for the DepMap dependency data.
## This formalises the checks in the student's own 01_import_check.R so
## that the QC evidence reported in Section 4.1 is reproducible.
##
## Three checks: (1) how many classified models have dependency data,
## (2) whether metadata and dependency rows are aligned, (3) how much
## dependency data is missing, at gene and cell-line level.

meta <- read.delim("data_raw/depmap_meta.txt", stringsAsFactors = FALSE, check.names = FALSE)
gd   <- read.delim("data_raw/depmap_GD.txt", row.names = 1,
                   stringsAsFactors = FALSE, check.names = FALSE)

cat("Metadata rows:", nrow(meta), "| Dependency matrix:", dim(gd), "\n\n")
cat("Cell-state classification, all metadata (stringent):\n")
print(table(meta$AM_class_stringent, useNA = "ifany"))

## ---- 1. Models that are both classified and measured -----------------------
matched_meta <- meta[!is.na(meta$AM_class_stringent) &
                       meta$condition %in% rownames(gd), ]
gd_classified <- gd[matched_meta$condition, , drop = FALSE]

cat("\nClassified models with dependency data:", nrow(matched_meta), "\n")
print(table(matched_meta$AM_class_stringent))

## ---- 2. Row alignment ------------------------------------------------------
cat("\nRow order identical to metadata:",
    identical(rownames(gd_classified), matched_meta$condition), "\n")

## ---- 3. Missingness --------------------------------------------------------
cat("Overall missing entries:", sum(is.na(gd_classified)),
    sprintf("(%.2f%%)\n", mean(is.na(gd_classified)) * 100))

missing_per_gene <- colSums(is.na(gd_classified))
gene_missing_summary <- data.frame(
  category = c("Complete", "Partly missing", "Entirely missing"),
  number_of_genes = c(sum(missing_per_gene == 0),
                      sum(missing_per_gene > 0 & missing_per_gene < nrow(gd_classified)),
                      sum(missing_per_gene == nrow(gd_classified))))
cat("\nGene-level missingness:\n"); print(gene_missing_summary)

cell_missing_summary <- data.frame(
  cell_line = rownames(gd_classified),
  cell_state = matched_meta$AM_class_stringent,
  missing_genes = rowSums(is.na(gd_classified)))
cat("\nCell-line missingness (range):",
    min(cell_missing_summary$missing_genes), "-",
    max(cell_missing_summary$missing_genes), "genes\n")

## ---- 4. Drop entirely-missing genes, then apply an 80% coverage rule -------
gd_available <- gd_classified[, missing_per_gene < nrow(gd_classified), drop = FALSE]
cat("\nGenes retained after dropping entirely-missing:", ncol(gd_available), "\n")

primary <- matched_meta$AM_class_stringent %in% c("ADRN", "MES")
sub  <- gd_available[matched_meta$condition[primary], , drop = FALSE]
cls  <- matched_meta$AM_class_stringent[primary]
n_adrn <- colSums(!is.na(sub[cls == "ADRN", , drop = FALSE]))
n_mes  <- colSums(!is.na(sub[cls == "MES",  , drop = FALSE]))

cat("\nGene coverage across the primary comparison groups:\n")
print(table(ADRN_observed = n_adrn, MES_observed = n_mes))

keep80 <- n_adrn >= ceiling(0.8 * sum(cls == "ADRN")) &
          n_mes  >= ceiling(0.8 * sum(cls == "MES"))
cat("\nGenes passing the 80% per-group coverage rule:", sum(keep80),
    "of", ncol(sub), "\n")

saveRDS(list(matched_meta = matched_meta, gd_available = gd_available,
             keep80 = names(keep80)[keep80]), "outputs/depmap_qc.rds")
cat("\nSaved: outputs/depmap_qc.rds\n")
