## 03c_aim2_tcf7l1.R
## Extension to Aim 2 following the secondary (continuous-score) Aim 1 analysis.
## Rule declared before testing: every TF that reaches FDR < 0.05 in the secondary
## Aim 1 analysis must be tested in the patient cohort. Five of the six were
## already among the 22 Aim 2 candidates; TCF7L1 was not. This script tests it
## with the same models as 03 (unadjusted) and 03b (adjusted for MYCN, INSS
## stage 4, age >= 18 months). It is a single pre-declared test, kept outside
## the 22-candidate FDR family so the main Aim 2 results are unchanged.
suppressMessages({ library(data.table); library(survival) })
rt <- function(p){ d <- fread(p, sep="\t", header=TRUE, data.table=FALSE, quote="\""); names(d)[1] <- "condition"; d }
meta <- rt("data_raw/KOCAK_meta.txt"); sv <- rt("data_raw/KOCAK_survival.txt"); tf <- rt("data_raw/KOCAK_TF_Activity.txt")
d <- merge(merge(meta, sv, by="condition"), tf, by="condition")
stopifnot("TCF7L1" %in% names(d))
d$MYCN <- factor(d$MYCN_bin); d$INSS <- factor(d$INSS_bin_4vall); d$AGE <- factor(d$Age_bin_18m)
d$x <- as.numeric(scale(d$TCF7L1))
fit1 <- function(tm, ev, adj){
  f <- if (adj) as.formula(sprintf("Surv(%s,%s) ~ x + MYCN + INSS + AGE", tm, ev)) else as.formula(sprintf("Surv(%s,%s) ~ x", tm, ev))
  m <- coxph(f, data=d); s <- summary(m)
  c(HR=s$coefficients["x","exp(coef)"], lo=s$conf.int["x","lower .95"], hi=s$conf.int["x","upper .95"], p=s$coefficients["x","Pr(>|z|)"], n=m$n)
}
out <- rbind(EFS_unadj=fit1("EFS_d","EFS_bin",FALSE), OS_unadj=fit1("OS_d","OS_bin",FALSE),
             EFS_adj=fit1("EFS_d","EFS_bin",TRUE),   OS_adj=fit1("OS_d","OS_bin",TRUE))
print(round(out, 4))
write.csv(data.frame(model=rownames(out), out), "outputs/aim2_cox_TCF7L1.csv", row.names=FALSE)

## How different would the candidate list have been if it had been built from the continuous analysis?
bin  <- read.csv("outputs/aim1_TF_dependency_limma_results.csv")$Gene[1:10]
cont <- read.csv("outputs/aim1_TF_dependency_continuous_results.csv")$Gene[1:10]
cat("\nTop-10 overlap (binary vs continuous):", paste(intersect(bin,cont), collapse=", "), "\n")
cat("Binary-only :", paste(setdiff(bin,cont), collapse=", "), "\n")
cat("Continuous-only:", paste(setdiff(cont,bin), collapse=", "), "\n")
