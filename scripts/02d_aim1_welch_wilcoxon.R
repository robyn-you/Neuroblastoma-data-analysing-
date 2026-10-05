## 02d_aim1_welch_wilcoxon.R
## Same Aim 1 differential dependency question as 02_aim1_differential_dependency.R,
## using the identical data, grouping, and TF panel -- only the statistical
## test differs. Structured to mirror 02's six steps exactly, so the two
## scripts can be read side by side. Requested by supervisor (Daniel) to
## check whether limma's empirical-Bayes moderation (suited to raw,
## high-throughput expression data) is appropriate for the already-processed
## Chronos dependency scores, or whether a plain Welch t-test / Wilcoxon
## rank-sum test is more suitable.

suppressMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("outputs", showWarnings = FALSE)

read_tab <- function(path, id_col = "condition") {
  df <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, quote = "\"")
  names(df)[1] <- id_col
  df
}

##  1. data load 
gd   <- read_tab("data_raw/depmap_GD.txt")
meta <- readRDS("outputs/depmap_meta.rds")
tf_panel <- read_tab("data_raw/depmap_TF_Activity.txt")

tf_genes <- setdiff(names(tf_panel), "condition")
cat("TF panel size:", length(tf_genes), "\n")

##  2. Match cell lines to classification 
meta_sub <- meta[, c("condition", "AM_class", "AM_class_stringent")]
gd_annot <- merge(meta_sub, gd, by = "condition")

keep <- gd_annot$AM_class_stringent %in% c("ADRN", "MES")
gd_use <- gd_annot[keep, ]
cat("Cell lines used for ADRN vs MES comparison:", nrow(gd_use), "\n")
print(table(gd_use$AM_class_stringent))

##  3. Build TF-restricted dependency matrix 
tf_genes_present <- intersect(tf_genes, names(gd_use))
expr_mat <- t(as.matrix(gd_use[, tf_genes_present]))   # genes x samples
colnames(expr_mat) <- gd_use$condition
complete_genes <- rowSums(is.na(expr_mat)) == 0
expr_mat <- expr_mat[complete_genes, ]
cat("TF genes with complete data across all used lines:", nrow(expr_mat), "\n")

group <- factor(gd_use$AM_class_stringent, levels = c("ADRN", "MES"))
adrn_cols <- colnames(expr_mat)[group == "ADRN"]
mes_cols  <- colnames(expr_mat)[group == "MES"]

##  4. Per-gene Welch t-test and Wilcoxon rank-sum test
## Same comparison direction as limma's coefficient (MES - ADRN): a positive
## difference means the gene is more essential in ADRN, negative = more
## essential in MES (Chronos scores are negative = essential).

run_welch <- function(x) {
  a <- x[adrn_cols]; m <- x[mes_cols]
  tt <- t.test(a, m)
  c(diff_MES_minus_ADRN = unname(tt$estimate[2] - tt$estimate[1]), P.Value = tt$p.value)
}
run_wilcox <- function(x) {
  a <- x[adrn_cols]; m <- x[mes_cols]
  wt <- suppressWarnings(wilcox.test(a, m, conf.int = TRUE))
  ## Hodges-Lehmann pseudo-median difference, matching the "MES - ADRN" sign
  c(diff_MES_minus_ADRN = unname(-wt$estimate), P.Value = wt$p.value)
}

welch_raw  <- t(apply(expr_mat, 1, run_welch))
wilcox_raw <- t(apply(expr_mat, 1, run_wilcox))

welch_res  <- data.frame(Gene = rownames(expr_mat), welch_raw, row.names = NULL)
wilcox_res <- data.frame(Gene = rownames(expr_mat), wilcox_raw, row.names = NULL)

welch_res$adj.P.Val  <- p.adjust(welch_res$P.Value,  method = "BH")
wilcox_res$adj.P.Val <- p.adjust(wilcox_res$P.Value, method = "BH")

welch_res  <- welch_res[order(welch_res$P.Value), ]
wilcox_res <- wilcox_res[order(wilcox_res$P.Value), ]

write.csv(welch_res,  "outputs/aim1_TF_dependency_welch_results.csv",  row.names = FALSE)
write.csv(wilcox_res, "outputs/aim1_TF_dependency_wilcoxon_results.csv", row.names = FALSE)

cat("\n== Top 15 differentially dependent TFs, Welch t-test ==\n")
print(head(welch_res, 15), digits = 3)
cat("\n== Top 15 differentially dependent TFs, Wilcoxon rank-sum ==\n")
print(head(wilcox_res, 15), digits = 3)

##  5. Sanity check against literature-known CRC TFs 
known_ADRN <- c("PHOX2B", "GATA3", "HAND2", "ISL1", "TBX2", "ASCL1")
known_MES  <- c("WWTR1", "FOSL2", "TEAD4", "PRRX1", "RUNX1", "RUNX2")

for (res_name in c("welch_res", "wilcox_res")) {
  res <- get(res_name)
  cat(sprintf("\n== %s: known ADRN-CRC TFs (expect positive diff) ==\n", res_name))
  print(res[res$Gene %in% known_ADRN, ], digits = 3)
  cat(sprintf("\n== %s: candidate MES-CRC TFs (expect negative diff) ==\n", res_name))
  print(res[res$Gene %in% known_MES, ], digits = 3)
}

##  6. Volcano plots 
make_volcano <- function(res, label, fname) {
  res$sig <- ifelse(res$adj.P.Val < 0.05, "FDR < 0.05", "n.s.")
  res$lbl <- ifelse(res$Gene %in% c(known_ADRN, known_MES), res$Gene, NA)
  p <- ggplot(res, aes(x = diff_MES_minus_ADRN, y = -log10(P.Value), color = sig)) +
    geom_point(alpha = 0.6, size = 1.3) +
    geom_text(aes(label = lbl), na.rm = TRUE, vjust = -0.6, size = 3, color = "black") +
  scale_color_manual(values = c("FDR < 0.05" = "red", "n.s." = "grey")) +
    labs(
      title = paste0("Aim 1: Differential TF dependency, ", label),
      x = "Dependency difference (MES - ADRN); positive = more essential in ADRN",
      y = expression(-log[10](p-value)),
      color = NULL
    ) +
    theme_bw(base_size = 12)
  ggsave(file.path("outputs/figures", fname), p, width = 8, height = 6, dpi = 150)
}

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
make_volcano(welch_res,  "Welch t-test",      "Figure_aim1_welch_volcano.png")
make_volcano(wilcox_res, "Wilcoxon rank-sum", "Figure_aim1_wilcoxon_volcano.png")

cat("\nSaved: outputs/aim1_TF_dependency_welch_results.csv\n")
cat("Saved: outputs/aim1_TF_dependency_wilcoxon_results.csv\n")
cat("Saved: outputs/figures/Figure_aim1_welch_volcano.png\n")
cat("Saved: outputs/figures/Figure_aim1_wilcoxon_volcano.png\n")
