## 02f_aim1_continuous_score.R
## SECONDARY / EXPLORATORY Aim 1 analysis.
##
## The primary Aim 1 analysis (02_aim1_differential_dependency.R) splits
## cell lines into ADRN vs MES (20 vs 5) and is underpowered. The literature
## describes ADRN-MES identity as a continuum with intermediate states
## (Bukkuri et al., 2024), and the supplied metadata already contains a
## continuous score (AM.score = MES-like minus ADRN-like). This script uses
## that score as a continuous predictor across ALL 28 classified lines with
## dependency data, including the 3 HYBRID lines, which avoids discarding
## information by dichotomising.
##
## DISCLOSURE: this analysis was run AFTER the binary analysis returned no
## FDR-significant hit, so it must be reported as secondary/exploratory and
## reported whatever it finds. The rationale (continuum, no dichotomisation)
## does not depend on the result.
##
## Interpretation of the slope: dependency change per +1 SD of AM.score
## (more MES-like). Chronos scores are negative = essential, so a POSITIVE
## slope means the gene is MORE essential in ADRN-like lines.

suppressMessages({ library(limma); library(ggplot2) })
dir.create("outputs", showWarnings = FALSE)

meta <- read.delim("data_raw/depmap_meta.txt", check.names = FALSE, stringsAsFactors = FALSE)
gd   <- read.delim("data_raw/depmap_GD.txt", row.names = 1, check.names = FALSE)
tfs  <- colnames(read.delim("data_raw/depmap_TF_Activity.txt", row.names = 1, check.names = FALSE))

mm <- meta[!is.na(meta$AM.score) & meta$condition %in% rownames(gd), ]
cat("Lines with AM.score and dependency data:", nrow(mm), "\n")
print(table(mm$AM_class_stringent))

panel <- intersect(tfs, colnames(gd))
X <- t(as.matrix(gd[mm$condition, panel]))
X <- X[rowSums(is.na(X)) == 0, ]
cat("TFs tested (same complete-case panel as primary analysis):", nrow(X), "\n")

## ---- limma with AM.score (per SD) as a continuous predictor ----------------
design <- model.matrix(~ scale(mm$AM.score))
fit <- eBayes(lmFit(X, design))
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
res$Gene <- rownames(res)
names(res)[names(res) == "logFC"] <- "slope_per_SD_AMscore"
res <- res[, c("Gene", "slope_per_SD_AMscore", "P.Value", "adj.P.Val")]
write.csv(res, "outputs/aim1_TF_dependency_continuous_results.csv", row.names = FALSE)

cat("\nFDR < 0.05:", sum(res$adj.P.Val < 0.05),
    "| raw p < 0.05:", sum(res$P.Value < 0.05),
    "(expected by chance:", round(0.05 * nrow(X)), ")\n")
print(head(res, 10), digits = 3, row.names = FALSE)

## ---- Robustness checks -------------------------------------------------------
hits <- res$Gene[res$adj.P.Val < 0.05]

# (a) without the HYBRID lines (25 ADRN/MES lines only)
keep <- mm$AM_class_stringent %in% c("ADRN", "MES")
fitA <- eBayes(lmFit(X[, keep], model.matrix(~ scale(mm$AM.score[keep]))))
rA <- topTable(fitA, coef = 2, number = Inf)
cat("\n(a) Excluding HYBRID lines: FDR<0.05 =", sum(rA$adj.P.Val < 0.05), "\n")

# (b) leave-one-line-out: does any single cell line drive a hit?
cat("\n(b) Leave-one-line-out: refits (of", nrow(mm), ") with raw p < 0.05\n")
for (g in hits) {
  ok <- sapply(seq_len(nrow(mm)), function(i) {
    d <- data.frame(y = X[g, -i], s = scale(mm$AM.score[-i])[, 1])
    summary(lm(y ~ s, d))$coefficients["s", "Pr(>|t|)"] < 0.05
  })
  cat(sprintf("  %-8s %2d / %d\n", g, sum(ok), nrow(mm)))
}

# (c) rank-based, outlier-resistant check
cat("\n(c) Spearman correlation with AM.score:\n")
for (g in hits) {
  ct <- suppressWarnings(cor.test(X[g, ], mm$AM.score, method = "spearman"))
  cat(sprintf("  %-8s rho = %+.2f, p = %.4f\n", g, ct$estimate, ct$p.value))
}
## ---- Figure: dependency against the continuous score for the hits ----------
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
plot_genes <- res$Gene[res$adj.P.Val < 0.05]
long <- do.call(rbind, lapply(plot_genes, function(g)
  data.frame(Gene = g, AM_score = mm$AM.score, Chronos = X[g, ],
             State = factor(mm$AM_class_stringent, levels = c("ADRN", "HYBRID", "MES")))))
long$Gene <- factor(long$Gene, levels = plot_genes)
p <- ggplot(long, aes(x = AM_score, y = Chronos)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  geom_point(aes(color = State), size = 2) +
  scale_color_manual(values = c(ADRN = "red", HYBRID = "grey", MES = "blue")) +
  facet_wrap(~Gene, scales = "free_y", ncol = 3) +
  labs(x = "AM score (more MES-like  >)", y = "Chronos dependency score", color = "Stringent class") +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave("outputs/figures/Figure_aim1_continuous_scatter.png", p, width = 8.5, height = 5.5, dpi = 300)

cat("\nSaved: outputs/aim1_TF_dependency_continuous_results.csv\n")
cat("Saved: outputs/figures/Figure_aim1_continuous_scatter.png\n")
