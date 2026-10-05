## 02_aim1_differential_dependency.R
## Aim 1: identify ADRN- vs MES-specific transcription factor gene
## dependencies from DepMap CRISPR Chronos scores.
##
## Design: moderated t-test (limma) comparing Chronos dependency scores
## of ADRN-classified vs MES-classified neuroblastoma cell lines,
## restricted to the transcription-factor gene panel, with BH FDR
## correction. This directly implements Phase 1 / Aim 1 of the proposal
## (Section 6.3.I).

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

##  1. Load data 
## NOTE: depmap_GD.txt and depmap_TF_Activity.txt are indexed by the short
## cell-line "condition" name (e.g. "CHP212"), NOT the DepMap ModelID
## (e.g. "ACH-000120"). depmap_meta.txt carries both, so we join on
## "condition" rather than ModelID.
gd   <- read_tab("data_raw/depmap_GD.txt")          # 37 lines x 18,444 genes (dependency)
meta <- readRDS("outputs/depmap_meta.rds")           # canonical classification (from script 01)
tf_panel <- read_tab("data_raw/depmap_TF_Activity.txt") # 34 lines x 763 TFs -> use colnames as TF list

tf_genes <- setdiff(names(tf_panel), "condition")
cat("TF panel size:", length(tf_genes), "\n")

##  2. Match cell lines to classification 
meta_sub <- meta[, c("condition", "AM_class", "AM_class_stringent")]
gd_annot <- merge(meta_sub, gd, by = "condition")

cat("Cell lines in GD matrix:", nrow(gd), "\n")
cat("Cell lines with classification available:", nrow(gd_annot), "\n")
cat("Class breakdown (stringent):\n")
print(table(gd_annot$AM_class_stringent, useNA = "ifany"))

## Use the stringent classification for the primary comparison (drops
## ambiguous HYBRID lines); ADRN vs MES only.
keep <- gd_annot$AM_class_stringent %in% c("ADRN", "MES")
gd_use <- gd_annot[keep, ]
cat("\nCell lines used for ADRN vs MES comparison:", nrow(gd_use), "\n")
print(table(gd_use$AM_class_stringent))

##  3. Build TF-restricted dependency matrix 
tf_genes_present <- intersect(tf_genes, names(gd_use))
cat("\nTF genes present in dependency matrix:", length(tf_genes_present),
    "of", length(tf_genes), "\n")

expr_mat <- t(as.matrix(gd_use[, tf_genes_present]))   # genes x samples
colnames(expr_mat) <- gd_use$condition

## Drop genes with any missing dependency scores (limma needs complete data)
complete_genes <- rowSums(is.na(expr_mat)) == 0
expr_mat <- expr_mat[complete_genes, ]
cat("TF genes with complete data across all used lines:", nrow(expr_mat), "\n")

##  4. limma differential dependency analysis 
group <- factor(gd_use$AM_class_stringent, levels = c("ADRN", "MES"))
design <- model.matrix(~group)

fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)

## coefficient 2 = MES - ADRN. Chronos: negative = more essential.
## iff > 0 => MES less dependent (ADRN more dependent); diff < 0 => MES more dependent
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
res$Gene <- rownames(res)
res <- res[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
names(res)[names(res) == "logFC"] <- "diff_MES_minus_ADRN"

write.csv(res, "outputs/aim1_TF_dependency_limma_results.csv", row.names = FALSE)

cat("\n== Top 15 differentially dependent TFs (by p-value) ==\n")
print(head(res, 15))

##  5. Sanity check against literature-known CRC TFs 
known_ADRN <- c("PHOX2B", "GATA3", "HAND2", "ISL1", "TBX2", "ASCL1")
known_MES  <- c("WWTR1", "FOSL2", "TEAD4", "PRRX1", "RUNX1", "RUNX2")

cat("\n== Known ADRN-CRC TFs (expect positive diff = more essential in ADRN) ==\n")
print(res[res$Gene %in% known_ADRN, ])

cat("\n== Candidate MES-CRC TFs (expect negative diff = more essential in MES) ==\n")
print(res[res$Gene %in% known_MES, ])

##  6. Volcano plot 
res$sig <- ifelse(res$adj.P.Val < 0.05, "FDR < 0.05", "n.s.")
res$label <- ifelse(res$Gene %in% c(known_ADRN, known_MES), res$Gene, NA)

res$group <- ifelse(res$Gene %in% known_ADRN, "Known ADRN TF",
              ifelse(res$Gene %in% known_MES, "Candidate MES TF", "Other TF"))

library(ggrepel)   # install.packages("ggrepel") 필요
p <- ggplot(res, aes(x = diff_MES_minus_ADRN, y = -log10(P.Value), color = group)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_text_repel(aes(label = label), na.rm = TRUE, size = 3, color = "black",
                  max.overlaps = Inf) +
  scale_color_manual(values = c("Known ADRN TF" = "red",
                                "Candidate MES TF" = "blue",
                                "Other TF" = "grey")) +
  labs(
    title = "Aim 1: Differential TF dependency, MES vs ADRN ",
    x = "Dependency difference (MES - ADRN); positive = more essential in ADRN",
    y = expression(-log[10](p-value)),
    color = NULL
  ) +
  theme_bw (base_size = 12)

ggsave("outputs/aim1_volcano_TF_dependency.png", p, width = 8, height = 6, dpi = 150)

cat("\nSaved: outputs/aim1_TF_dependency_limma_results.csv\n")
cat("Saved: outputs/aim1_volcano_TF_dependency.png\n")
