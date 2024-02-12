##================================
##  video 1
##================================

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

###### A Single Heatmap

set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)

mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))

Heatmap(mat)

#################################################################
### Legend title
Heatmap(mat)
# default color schema is “blue-white-red” 
# dendrograms, the row/column names and the heatmap legend
# legend title is assigned with matrix + an internal index number.

# legend title
Heatmap(mat, name = "mat")

### Row and Column titles
# row_title_gp and column_title_gp
Heatmap(mat, name = "mat", 
        column_title = "I am a column title", 
        row_title = "I am a row title")

# column title at the bottom
Heatmap(mat, name = "mat", 
        column_title = "column title at the bottom",
        column_title_side = "bottom")

# Rotations row/column titles, set by row_title_rot and column_title_rot
Heatmap(mat, name = "mat", row_title = "row title", row_title_rot = 0)
Heatmap(mat, name = "mat", column_title = "column title", column_title_rot = 90)

# use gpar () to specify graphics parameters
Heatmap(mat, name = "mat", 
        row_title = "big row title", 
        row_title_gp = gpar(fontsize = 20, fontface = "bold"),
        column_title = "big column title", 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"))

Heatmap(mat, name = "mat", 
        column_title = "I am a column title", 
        column_title_gp = gpar(col = "white", fill = "red", border = "blue"))

# adjust the space
ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
Heatmap(mat, name = "mat", column_title = "I am a column title", 
        column_title_gp = gpar(col = "white", fill = "red", border = "blue"))

ht_opt(RESET = TRUE)

# Title can be set as mathematical formulas
Heatmap(mat, name = "mat", 
        column_title = expression(hat(beta) == (X^t * X)^{-1} * X^t * y)) 

# using the gridtext package to draw more complicated text title
library(gridtext)
Heatmap(mat, name = "mat",
        column_title = gt_render(
          paste0("Some <span style='color:blue'>blue text **in bold.**</span><br>",
                 "And *italics text.*<br>And some ",
                 "<span style='font-size:18pt; color:black'>large</span> text."), 
          r = unit(2, "pt"), 
          padding = unit(c(2, 2, 2, 2), "pt")))



