## 03_aim2_survival_analysis.R
## Aim 2: evaluate clinical relevance of candidate transcription factors
## from Aim 1 using overall survival (OS) and event-free survival (EFS)
## in the KOCAK neuroblastoma cohort (Section 6.3.II of the proposal).
##
## For each candidate TF, activity is tested as a continuous predictor
## in a Cox proportional hazards model (per-SD scaled), for both EFS and
## OS, with BH FDR correction applied across all tests. Kaplan-Meier
## curves (median-split) are produced for the top hits for visualization.

suppressMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
})

dir.create("outputs", showWarnings = FALSE)

read_tab <- function(path, id_col = "condition") {
  df <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, quote = "\"")
  names(df)[1] <- id_col
  df
}

## ---- 1. Load and merge KOCAK data -------------------------------------------
meta     <- read_tab("data_raw/KOCAK_meta.txt")
survdat  <- read_tab("data_raw/KOCAK_survival.txt")
tf_act   <- read_tab("data_raw/KOCAK_TF_Activity.txt")

merged <- merge(meta, survdat, by = "condition")
merged <- merge(merged, tf_act, by = "condition")
cat("Patients with meta + survival + TF activity data:", nrow(merged), "\n")
cat("Class breakdown (stringent):\n")
print(table(merged$AM_class_stringent, useNA = "ifany"))

## ---- 2. Candidate TF list: literature CRC TFs + top Aim 1 hits --------------
literature_TFs <- c("PHOX2B", "GATA3", "HAND2", "ISL1", "TBX2", "ASCL1",
                     "WWTR1", "FOSL2", "TEAD4", "PRRX1", "RUNX1", "RUNX2")
aim1_top_hits  <- c("ATF5", "MAFA", "TFAP2D", "NFIC", "JUNB", "STAT1",
                     "HMGA2", "IRF7", "MAZ", "HHEX")
candidate_TFs <- unique(c(literature_TFs, aim1_top_hits))
candidate_TFs <- intersect(candidate_TFs, names(tf_act))
cat("\nCandidate TFs available in KOCAK TF activity data:", length(candidate_TFs), "\n")
print(candidate_TFs)

## ---- 3. Cox proportional hazards, per-SD-scaled activity, for EFS and OS ----
run_cox <- function(gene, time_col, event_col, data) {
  x <- scale(data[[gene]])[, 1]  # per-SD scale for comparable HRs across genes
  df <- data.frame(time = data[[time_col]], event = data[[event_col]], x = x)
  df <- df[complete.cases(df), ]
  fit <- tryCatch(coxph(Surv(time, event) ~ x, data = df), error = function(e) NULL)
  if (is.null(fit)) return(data.frame(Gene = gene, HR = NA, lower95 = NA,
                                       upper95 = NA, P.Value = NA, n = nrow(df)))
  s <- summary(fit)
  data.frame(
    Gene = gene,
    HR = s$coefficients[1, "exp(coef)"],
    lower95 = s$conf.int[1, "lower .95"],
    upper95 = s$conf.int[1, "upper .95"],
    P.Value = s$coefficients[1, "Pr(>|z|)"],
    n = nrow(df)
  )
}

efs_results <- do.call(rbind, lapply(candidate_TFs, run_cox,
                                      time_col = "EFS_d", event_col = "EFS_bin",
                                      data = merged))
efs_results$adj.P.Val <- p.adjust(efs_results$P.Value, method = "BH")
efs_results <- efs_results[order(efs_results$P.Value), ]

os_results <- do.call(rbind, lapply(candidate_TFs, run_cox,
                                     time_col = "OS_d", event_col = "OS_bin",
                                     data = merged))
os_results$adj.P.Val <- p.adjust(os_results$P.Value, method = "BH")
os_results <- os_results[order(os_results$P.Value), ]

cat("\n== Cox regression: Event-Free Survival (per-SD TF activity) ==\n")
print(efs_results, row.names = FALSE)

cat("\n== Cox regression: Overall Survival (per-SD TF activity) ==\n")
print(os_results, row.names = FALSE)

write.csv(efs_results, "outputs/aim2_cox_EFS_results.csv", row.names = FALSE)
write.csv(os_results, "outputs/aim2_cox_OS_results.csv", row.names = FALSE)

## ---- 4. Kaplan-Meier curves (median split) for top hits ---------------------
plot_km <- function(gene, time_col, event_col, data, label) {
  x <- data[[gene]]
  grp <- factor(ifelse(x > median(x, na.rm = TRUE), "High activity", "Low activity"),
                levels = c("Low activity", "High activity"))
  df <- data.frame(time = data[[time_col]], event = data[[event_col]], group = grp)
  df <- df[complete.cases(df), ]
  fit <- survfit(Surv(time, event) ~ group, data = df)
  p <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE,
                   title = paste0(gene, " — ", label),
                   xlab = "Days", legend.title = gene,
                   palette = c("#2C7FB8", "#D95F02"))
  p
}

top_efs_genes <- head(efs_results$Gene[!is.na(efs_results$P.Value)], 3)
for (g in top_efs_genes) {
  p <- plot_km(g, "EFS_d", "EFS_bin", merged, "Event-Free Survival")
  png(paste0("outputs/aim2_KM_EFS_", g, ".png"), width = 7, height = 7,
      units = "in", res = 150)
  print(p)
  dev.off()
}

cat("\nSaved: outputs/aim2_cox_EFS_results.csv\n")
cat("Saved: outputs/aim2_cox_OS_results.csv\n")
cat("Saved KM plots for top EFS hits:", paste(top_efs_genes, collapse = ", "), "\n")
