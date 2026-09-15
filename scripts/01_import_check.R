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