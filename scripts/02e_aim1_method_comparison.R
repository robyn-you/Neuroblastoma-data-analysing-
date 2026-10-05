## 02e_aim1_method_comparison.R
## Analyses the output of 02_aim1_differential_dependency.R (limma) against
## 02d_aim1_welch_wilcoxon.R (Welch t-test, Wilcoxon rank-sum). Produces:
## (1) counts of significant genes per method, at raw and FDR-adjusted level;
## (2) the theoretical floor on achievable significance at this sample size;
## (3) a merged per-gene table for the literature candidate TFs;
## (4) an agreement plot (limma vs each alternative);
## (5) a list of genes where the three methods disagree on significance.

suppressMessages({
  library(ggplot2)
})

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load each method's results -----------------------------------------
limma  <- read.csv("outputs/aim1_TF_dependency_limma_results.csv")
welch  <- read.csv("outputs/aim1_TF_dependency_welch_results.csv")
wilcox <- read.csv("outputs/aim1_TF_dependency_wilcoxon_results.csv")

## ---- 2. Overall counts: how many genes are significant per method? --------
cat("== Genes with raw p < 0.05 ==\n")
cat("  limma   :", sum(limma$P.Value  < 0.05), "/", nrow(limma),  "\n")
cat("  Welch t :", sum(welch$P.Value  < 0.05), "/", nrow(welch),  "\n")
cat("  Wilcoxon:", sum(wilcox$P.Value < 0.05), "/", nrow(wilcox), "\n")

cat("\n== Genes with FDR (adj.P.Val) < 0.05 ==\n")
cat("  limma   :", sum(limma$adj.P.Val  < 0.05), "\n")
cat("  Welch t :", sum(welch$adj.P.Val  < 0.05), "\n")
cat("  Wilcoxon:", sum(wilcox$adj.P.Val < 0.05), "\n")

## ---- 3. Theoretical floor at this sample size (20 ADRN vs 5 MES) -----------
## Even with perfect separation between groups, this is the smallest p-value
## the rank-sum test can produce -- a hard ceiling set by sample size alone,
## not by which test is used.
floor_w <- wilcox.test(1:20, 25 + 1:5)$p.value
cat("\nSmallest possible Wilcoxon p-value at n=20 vs n=5 (perfect separation):",
    signif(floor_w, 3), "\n")
cat("(", nrow(limma), "genes tested, so this floor sits close to the FDR",
    "significance boundary even in the best case)\n")

## ---- 4. Merge all three methods into one table -----------------------------
merged <- Reduce(function(x, y) merge(x, y, by = "Gene", suffixes = c("", "")),
                  list(
                    setNames(limma[,  c("Gene", "diff_MES_minus_ADRN", "P.Value", "adj.P.Val")],
                             c("Gene", "diff_limma",  "P_limma",  "FDR_limma")),
                    setNames(welch[,  c("Gene", "diff_MES_minus_ADRN", "P.Value", "adj.P.Val")],
                             c("Gene", "diff_welch",  "P_welch",  "FDR_welch")),
                    setNames(wilcox[, c("Gene", "diff_MES_minus_ADRN", "P.Value", "adj.P.Val")],
                             c("Gene", "diff_wilcox", "P_wilcox", "FDR_wilcox"))
                  ))
merged$n_methods_nominal_sig <- rowSums(cbind(merged$P_limma   < 0.05,
                                               merged$P_welch   < 0.05,
                                               merged$P_wilcox  < 0.05))
merged <- merged[order(merged$P_limma), ]
write.csv(merged, "outputs/aim1_TF_dependency_method_comparison.csv", row.names = FALSE)

known_ADRN <- c("PHOX2B", "GATA3", "HAND2", "ISL1", "TBX2", "ASCL1")
known_MES  <- c("WWTR1", "FOSL2", "TEAD4", "PRRX1", "RUNX1", "RUNX2")
cat("\n== Literature CRC TFs across all three methods (raw p-values) ==\n")
print(merged[merged$Gene %in% c(known_ADRN, known_MES),
             c("Gene", "P_limma", "P_welch", "P_wilcox", "n_methods_nominal_sig")],
      row.names = FALSE, digits = 3)

## ---- 5. Genes where the methods disagree on nominal significance ----------
disagree <- merged[merged$n_methods_nominal_sig > 0 & merged$n_methods_nominal_sig < 3, ]
cat("\n== Genes where methods disagree on nominal (p<0.05) significance:",
    nrow(disagree), "of", nrow(merged), "==\n")
cat("(This is expected at n=5 MES -- a single borderline observation can\n")
cat(" push one test's p-value just under 0.05 while another stays just above.)\n")

## ---- 6. Agreement plot: limma vs each alternative --------------------------
long <- rbind(
  data.frame(Gene = merged$Gene, limma = -log10(merged$P_limma),
             alt = -log10(merged$P_welch),  Method = "Welch t-test"),
  data.frame(Gene = merged$Gene, limma = -log10(merged$P_limma),
             alt = -log10(merged$P_wilcox), Method = "Wilcoxon rank-sum")
)
p <- ggplot(long, aes(x = limma, y = alt)) +
  geom_point(alpha = 0.35, size = 1.2, color = "#2C7FB8") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~Method) +
  labs(title = "Aim 1: agreement between limma and alternative tests",
       subtitle = "Dashed line = identical p-value; points cluster near the line",
       x = expression(-log[10](p)~"(limma)"), y = expression(-log[10](p)~"(alternative test)")) +
  theme_bw(base_size = 12)
ggsave("outputs/figures/Figure_aim1_method_agreement.png", p, width = 9, height = 4.5, dpi = 300)

cat("\nSaved: outputs/aim1_TF_dependency_method_comparison.csv\n")
cat("Saved: outputs/figures/Figure_aim1_method_agreement.png\n")
