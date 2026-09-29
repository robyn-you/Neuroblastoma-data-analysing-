## 04_crossdataset_integration.R
## Week 10 deliverable: consolidate Aim 1 (CRISPR dependency) and Aim 2
## (KOCAK survival) evidence into a single prioritised transcription
## factor list, per Section 8 of the proposal ("Rank candidate TFs based
## on convergent evidence: dependency significance, effect size,
## prognostic significance").

suppressMessages({
  library(ggplot2)
  library(reshape2)
})

dir.create("outputs", showWarnings = FALSE)

## ---- 1. Load results from Aim 1 and Aim 2 -----------------------------------
dep <- read.csv("outputs/aim1_TF_dependency_limma_results.csv")
efs <- read.csv("outputs/aim2_cox_EFS_results.csv")
os  <- read.csv("outputs/aim2_cox_OS_results.csv")

## Restrict to the candidate universe actually tested in Aim 2
## (literature ADRN/MES-CRC TFs + top Aim 1 hits)
candidate_TFs <- efs$Gene
dep <- dep[dep$Gene %in% candidate_TFs, c("Gene", "diff_MES_minus_ADRN", "P.Value", "adj.P.Val")]
names(dep) <- c("Gene", "dep_diff_MES_minus_ADRN", "dep_P", "dep_FDR")

efs <- efs[, c("Gene", "HR", "P.Value", "adj.P.Val")]
names(efs) <- c("Gene", "EFS_HR", "EFS_P", "EFS_FDR")

os <- os[, c("Gene", "HR", "P.Value", "adj.P.Val")]
names(os) <- c("Gene", "OS_HR", "OS_P", "OS_FDR")

## ---- 2. Merge into one integrated evidence table ----------------------------
integrated <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE),
                      list(dep, efs, os))

## ---- 3. Composite ranking score ----------------------------------------------
## Raw p-values (not FDR) are used for the composite score because the
## three FDR corrections were computed over different denominators
## (743 genome-wide TFs for Aim 1 vs. 22 candidates for Aim 2) and are
## not directly comparable; adjusted p-values per test are retained in
## the output table for transparency.
integrated$composite_score <- -log10(integrated$dep_P) +
                               -log10(integrated$EFS_P) +
                               -log10(integrated$OS_P)

integrated$n_nominal_sig <- rowSums(cbind(
  integrated$dep_P < 0.05,
  integrated$EFS_P < 0.05,
  integrated$OS_P  < 0.05
), na.rm = TRUE)

integrated$n_FDR_sig <- rowSums(cbind(
  integrated$dep_FDR < 0.05,
  integrated$EFS_FDR < 0.05,
  integrated$OS_FDR  < 0.05
), na.rm = TRUE)

integrated <- integrated[order(-integrated$composite_score), ]

write.csv(integrated, "outputs/integration_TF_ranking.csv", row.names = FALSE)

cat("== Integrated TF ranking (Aim 1 dependency + Aim 2 EFS/OS survival) ==\n")
print(integrated[, c("Gene", "composite_score", "n_nominal_sig", "n_FDR_sig",
                      "dep_P", "EFS_P", "OS_P")], row.names = FALSE, digits = 3)

## ---- 4. Convergence heatmap --------------------------------------------------
plot_df <- integrated[, c("Gene", "dep_P", "EFS_P", "OS_P")]
plot_df$Gene <- factor(plot_df$Gene, levels = rev(integrated$Gene))  # keep rank order
long <- melt(plot_df, id.vars = "Gene", variable.name = "Test", value.name = "P")
long$Test <- factor(long$Test, levels = c("dep_P", "EFS_P", "OS_P"),
                     labels = c("Dependency\n(ADRN vs MES)", "EFS", "OS"))
long$neg_log10_p <- -log10(long$P)
long$sig_label <- ifelse(long$P < 0.001, "***",
                   ifelse(long$P < 0.01, "**",
                   ifelse(long$P < 0.05, "*", "")))

p <- ggplot(long, aes(x = Test, y = Gene, fill = neg_log10_p)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 4, vjust = 0.7) +
  scale_fill_gradient(low = "white", high = "firebrick",
                       name = expression(-log[10](p))) +
  labs(title = "Integrated evidence: gene dependency vs. clinical survival",
       subtitle = "Genes ranked top-to-bottom by combined evidence strength (Aim 1 + Aim 2)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face = "bold"))

ggsave("outputs/integration_heatmap.png", p, width = 7, height = 8, dpi = 150)

cat("\nSaved: outputs/integration_TF_ranking.csv\n")
cat("Saved: outputs/integration_heatmap.png\n")
