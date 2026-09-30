## 02b_aim1_sensitivity_nonstringent.R
## Sensitivity analysis for Aim 1 (proposal Section 6.4 risk mitigation:
## "sensitivity analysis will be run with alternative classification
## thresholds to confirm robustness of findings").
##
## Same limma differential-dependency test as 02_aim1_differential_dependency.R,
## but using the non-stringent AM_class instead of AM_class_stringent. The
## stringent call reassigns borderline lines to HYBRID and drops them; the
## non-stringent call forces every line into ADRN or MES, which should
## recover a larger (if noisier) MES group.

suppressMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

dir.create("outputs", showWarnings = FALSE)

read_tab <- function(path, id_col = "condition") {
  df <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, quote = "\"")
  names(df)[1] <- id_col
  df
}

# 1. Load data (identical sources to script 02) 
gd   <- read_tab("data_raw/depmap_GD.txt")
meta <- readRDS("outputs/depmap_meta.rds")
tf_panel <- read_tab("data_raw/depmap_TF_Activity.txt")
tf_genes <- setdiff(names(tf_panel), "condition")

meta_sub <- meta[, c("condition", "AM_class", "AM_class_stringent")]
gd_annot <- merge(meta_sub, gd, by = "condition")

cat("Class breakdown, NON-stringent (AM_class):\n")
print(table(gd_annot$AM_class, useNA = "ifany"))
cat("\nClass breakdown, stringent (for comparison):\n")
print(table(gd_annot$AM_class_stringent, useNA = "ifany"))

##  2. Subset to ADRN/MES using non-stringent AM_class 
keep <- gd_annot$AM_class %in% c("ADRN", "MES")
gd_use <- gd_annot[keep, ]
cat("\nCell lines used (non-stringent):", nrow(gd_use), "\n")
print(table(gd_use$AM_class))

tf_genes_present <- intersect(tf_genes, names(gd_use))
expr_mat <- t(as.matrix(gd_use[, tf_genes_present]))
colnames(expr_mat) <- gd_use$condition
expr_mat <- expr_mat[rowSums(is.na(expr_mat)) == 0, ]
cat("TF genes with complete data:", nrow(expr_mat), "\n")

##  3. limma, same design as before 
group <- factor(gd_use$AM_class, levels = c("ADRN", "MES"))
design <- model.matrix(~group)
fit <- eBayes(lmFit(expr_mat, design))
res_ns <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
res_ns$Gene <- rownames(res_ns)
res_ns <- res_ns[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
names(res_ns)[names(res_ns) == "logFC"] <- "diff_MES_minus_ADRN"

write.csv(res_ns, "outputs/aim1_TF_dependency_limma_NONSTRINGENT.csv", row.names = FALSE)

cat("\n== Top 15, non-stringent classification ==\n")
print(head(res_ns, 15))

##  4. Side-by-side comparison with the stringent (primary) result 
res_strict <- read.csv("outputs/aim1_TF_dependency_limma_results.csv")

compare <- merge(
  res_strict[, c("Gene", "diff_MES_minus_ADRN", "P.Value", "adj.P.Val")],
  res_ns[,     c("Gene", "diff_MES_minus_ADRN", "P.Value", "adj.P.Val")],
  by = "Gene", suffixes = c("_stringent", "_nonstringent")
)

known_TFs <- c("PHOX2B", "GATA3", "HAND2", "ISL1", "TBX2", "ASCL1",
               "WWTR1", "FOSL2", "TEAD4", "PRRX1", "RUNX1", "RUNX2")

cat("\n== Literature CRC TFs: stringent vs non-stringent side by side ==\n")
print(compare[compare$Gene %in% known_TFs, ], row.names = FALSE)

cat("\nGenes newly significant at FDR<0.05 in non-stringent but not stringent:\n")
newly_sig <- compare[compare$adj.P.Val_nonstringent < 0.05 &
                        (is.na(compare$adj.P.Val_stringent) | compare$adj.P.Val_stringent >= 0.05), ]
print(newly_sig[order(newly_sig$adj.P.Val_nonstringent), ], row.names = FALSE)

cat("\nSaved: outputs/aim1_TF_dependency_limma_NONSTRINGENT.csv\n")
