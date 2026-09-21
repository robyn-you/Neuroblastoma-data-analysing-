#cellline data import
depmap_meta <- read.delim(
  "data_raw/depmap_meta.txt",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

#dimensions_preview_data_firsst_six_celllines
dim(depmap_meta)
head(depmap_meta)
#count_cell-state_classification
table(depmap_meta$AM_class_stringent, useNA = "ifany")


#dependency data import, full size
depmap_gd <- read.delim(
  "data_raw/depmap_GD.txt",
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(depmap_gd)

#example preview,
#ADRN(CHP134) # nolint
#HYBRID (SKINMM),
#MES (CHLA90) cell lines,
#PHOX2B gene
depmap_gd [
  c("CHP134", "SKNMM", "CHLA90"),
  "PHOX2B",
  drop = FALSE
]

#select metadata for clssified models with dependency
matched_meta <- depmap_meta[
  !is.na(depmap_meta$AM_class_stringent) &
    depmap_meta$condition %in% rownames(depmap_gd),
]

#select dependency
depmap_gd_classified <- depmap_gd[
  matched_meta$condition,
  ,
  drop = FALSE
]
#dimension and group counts
dim(depmap_gd_classified)
table(matched_meta$AM_class_stringent, useNA = "ifany")
identical(rownames(depmap_gd_classified),
          matched_meta$condition)

# check missing dependency measurments in the matched dataset
sum(is.na(depmap_gd_classified)) #identifies missing entries
mean(is.na(depmap_gd_classified)) * 100 #% check
missing_per_gene <- colSums(is.na(depmap_gd_classified))

#missing gene
gene_missing_summary <- data.frame(
  category = c("Complete", "Partly missing", "Entirely missing"),
  number_of_genes = c(
    sum(missing_per_gene == 0),
    sum(missing_per_gene > 0 &
          missing_per_gene < nrow(depmap_gd_classified))
  )
)

print(gene_missing_summary)

# missing measurement count
cell_missing_summary <- data.frame(
  cell_line = rownames(depmap_gd_classified),
  cell_state = matched_meta$AM_class_stringent,
  missing_genes = rowSums(is.na(depmap_gd_classified))
)
print(cell_missing_summary)
