## 03b_aim2_adjusted_cox.R
## Robustness check for Aim 2: are TF-survival associations independent of
## established clinical risk factors? Each candidate TF is re-tested in a
## multivariable Cox model adjusted for MYCN amplification, INSS stage 4
## vs other, and age at diagnosis (<18 vs >=18 months).

suppressMessages({
  library(data.table)
  library(survival)
})

read_tab <- function(path, id_col = "condition") {
  df <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, quote = "\"")
  names(df)[1] <- id_col
  df
}

meta    <- read_tab("data_raw/KOCAK_meta.txt")
survdat <- read_tab("data_raw/KOCAK_survival.txt")
tf_act  <- read_tab("data_raw/KOCAK_TF_Activity.txt")
merged  <- merge(merge(meta, survdat, by = "condition"), tf_act, by = "condition")

merged$MYCN_bin      <- factor(merged$MYCN_bin)
merged$INSS_bin_4vall <- factor(merged$INSS_bin_4vall)
merged$Age_bin_18m   <- factor(merged$Age_bin_18m)

cat("Covariate coding check:\n")
print(table(merged$MYCN_status, merged$MYCN_bin, useNA = "ifany"))
print(table(merged$INSS, merged$INSS_bin_4vall, useNA = "ifany"))

candidates <- read.csv("outputs/aim2_cox_EFS_results.csv")$Gene

run_adj <- function(gene, time_col, event_col) {
  df <- data.frame(time = merged[[time_col]], event = merged[[event_col]],
                   x = scale(merged[[gene]])[, 1],
                   MYCN = merged$MYCN_bin, INSS = merged$INSS_bin_4vall,
                   AGE = merged$Age_bin_18m)
  df <- df[complete.cases(df), ]
  fit <- coxph(Surv(time, event) ~ x + MYCN + INSS + AGE, data = df)
  s <- summary(fit)
  data.frame(Gene = gene,
             adj_HR = s$coefficients["x", "exp(coef)"],
             lower95 = s$conf.int["x", "lower .95"],
             upper95 = s$conf.int["x", "upper .95"],
             P.Value = s$coefficients["x", "Pr(>|z|)"],
             n = nrow(df))
}

for (ep in list(c("EFS_d", "EFS_bin", "EFS"), c("OS_d", "OS_bin", "OS"))) {
  res <- do.call(rbind, lapply(candidates, run_adj, time_col = ep[1], event_col = ep[2]))
  res$adj.P.Val <- p.adjust(res$P.Value, method = "BH")
  res <- res[order(res$P.Value), ]
  cat("\n== Adjusted Cox (MYCN + INSS + age):", ep[3], "==\n")
  print(res, row.names = FALSE, digits = 3)
  write.csv(res, paste0("outputs/aim2_cox_", ep[3], "_adjusted.csv"), row.names = FALSE)
}
