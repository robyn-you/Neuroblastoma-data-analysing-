## 02c_dependency_heatmap.R
## Supplementary figure: raw Chronos dependency scores for the top candidate
## TFs across all 25 ADRN/MES cell lines (not just summary statistics), so the
## reader can see why n = 5 MES lines limits power.
## Uses base R heatmap(), as in the student's own version, with five fixes:
##  1. the output folder is created before saving;
##  2. the colour scale is centred on 0 (white = no effect) and saturates at
##     +/-0.8, so a few extreme values (MYC in CHLA15, -2.4) do not wash out
##     the rest. Clamping affects the plotted copy only, not the data;
##  3. a colour key and a cell-state legend are drawn;
##  4. the title has its own margin and is no longer cut off;
##  5. rows run ATF5 (top) to MYC (bottom), matching the report's ranking.

suppressMessages({ library(data.table) })

read_tab <- function(path, id_col = "condition") {
  df <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, quote = "\"")
  names(df)[1] <- id_col
  df
}

meta <- fread("data_raw/depmap_meta.txt", sep = "\t", header = TRUE,
              data.table = FALSE, quote = "\"")
gd   <- read_tab("data_raw/depmap_GD.txt")
rownames(gd) <- gd$condition; gd$condition <- NULL

pm <- meta[meta$AM_class_stringent %in% c("ADRN", "MES") &
             meta$condition %in% rownames(gd), ]
pm <- pm[order(pm$AM_class_stringent), ]                  # ADRN lines, then MES

genes <- c("ATF5", "HAND2", "RUNX1", "HHEX", "MAZ", "GATA3",
           "WWTR1", "RUNX2", "MYCN", "MYC")
genes <- intersect(genes, colnames(gd))

mat <- t(as.matrix(gd[pm$condition, genes, drop = FALSE]))
colnames(mat) <- pm$condition
cat("Lines:", ncol(mat), "| genes:", nrow(mat),
    "| lines with a missing value:", sum(colSums(is.na(mat)) > 0), "\n")
cat("Chronos range in this figure:", round(range(mat, na.rm = TRUE), 2), "\n")

## ---- colour scale ------------------------------------------------------------
lim  <- 0.8
cols <- colorRampPalette(c("darkblue", "white", "darkred"))(100)  # blue = essential
mat_plot <- pmax(pmin(mat, lim), -lim)
mat_plot <- mat_plot[nrow(mat_plot):1, , drop = FALSE]    # base heatmap draws row 1 at the bottom
state_col <- ifelse(pm$AM_class_stringent == "ADRN", "tomato", "blue")

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
png("outputs/figures/Figure_dependency_heatmap.png", width = 9, height = 5.8,
    units = "in", res = 300)
par(oma = c(0, 0, 2.5, 0))                                # room for the title
heatmap(mat_plot, Rowv = NA, Colv = NA, ColSideColors = state_col,
        col = cols, breaks = seq(-lim, lim, length.out = 101), scale = "none",
        cexRow = 1, cexCol = 0.75, margins = c(6, 6))
mtext("Dependency of top candidate TFs (25 ADRN/MES lines)", side = 3,
      outer = TRUE, line = 0.5, font = 2, cex = 1.1)

## colour key
par(fig = c(0.87, 0.90, 0.30, 0.62), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
image(1, seq(-lim, lim, length.out = 100), matrix(seq(-lim, lim, length.out = 100), nrow = 1),
      col = cols, axes = FALSE, xlab = "", ylab = "")
axis(4, at = c(-lim, -0.4, 0, 0.4, lim), las = 1, cex.axis = 0.75,
     labels = c("-0.8", "-0.4", "0", "0.4", "0.8"))
mtext("Chronos score", side = 3, line = 1.4, cex = 0.8, font = 2, adj = 0)
mtext("blue = essential", side = 3, line = 0.2, cex = 0.7, adj = 0)
mtext("beyond +/-0.8 shown\nas the end colour", side = 1, line = 1.2, cex = 0.65, adj = 0)

## cell-state legend
par(fig = c(0.84, 0.99, 0.70, 0.88), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("ADRN", "MES"), fill = c("tomato", "blue"),
       title = "Cell state", bty = "n", cex = 0.85)
dev.off()
cat("Saved: outputs/figures/Figure_dependency_heatmap.png\n")
