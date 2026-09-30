## 01_data_import
## Purpose: confirm data_old (previously supplied "clean" files) and
## data_raw (newer supervisor export) agree with each other before we
## commit to one as the analysis source. We use data_raw as canonical
## going forward because it carries extra clinical covariates
## (INSS stage, Age_bin) and MYCN copy number that data_old lacks.

suppressMessages({
  library(data.table)
})

read_tab <- function(path) {
  fread(path, sep = "\t", header = TRUE, data.table = FALSE, quote = "\"")
}

##  1. DepMap cell line classification: old vs raw 
old_depmap_meta <- read_tab("data_old/depmap_meta_nb_clean.txt")
raw_depmap_meta <- read_tab("data_raw/depmap_meta.txt")

names(old_depmap_meta)[1] <- "ModelID"
names(raw_depmap_meta)[1] <- "ModelID"

cmp_depmap <- merge(
  old_depmap_meta[, c("ModelID", "AM_class", "AM_class_stringent")],
  raw_depmap_meta[, c("ModelID", "AM_class", "AM_class_stringent")],
  by = "ModelID", suffixes = c("_old", "_raw")
)
cat("== DepMap meta: rows matched =", nrow(cmp_depmap), "==\n")
cat("AM_class agreement (old vs raw):\n")
print(table(cmp_depmap$AM_class_old, cmp_depmap$AM_class_raw, useNA = "ifany"))

##  2. KOCAK patient classification: old vs raw 
old_kocak_meta <- read_tab("data_old/KOCAK_meta_nb_clean.txt")
raw_kocak_meta <- read_tab("data_raw/KOCAK_meta.txt")

cmp_kocak <- merge(
  old_kocak_meta[, c("condition", "AM_class", "AM_class_stringent")],
  raw_kocak_meta[, c("condition", "AM_class", "AM_class_stringent")],
  by = "condition", suffixes = c("_old", "_raw")
)
cat("\n== KOCAK meta: rows matched =", nrow(cmp_kocak), "==\n")
cat("AM_class agreement (old vs raw):\n")
print(table(cmp_kocak$AM_class_old, cmp_kocak$AM_class_raw, useNA = "ifany"))

##  3. Save canonical objects for downstream scripts 
dir.create("outputs", showWarnings = FALSE)
saveRDS(raw_depmap_meta, "outputs/depmap_meta.rds")
saveRDS(raw_kocak_meta, "outputs/kocak_meta.rds")

cat("\nSaved canonical meta tables to outputs/. Using data_raw as the source of truth.\n")
