################################################################################
# Write a list of data.frames in .csv tables.
# Input:
# tablelist: a list of data.frames.
# mydir: the directory where the tables will be saved.
# Output:
# .csv tables saved in mydir.
################################################################################

#' Write a list of data.frames in .csv tables
#'
#' @param tablelist a list of data.frames.
#' @param mydir the directory where the tables will be saved.
#'
#' @return
#' .csv tables saved in mydir
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


writeLRtable <- function(tablelist, mydir) {

  for (i in 1:length(tablelist)) {
    fwrite(tablelist[[i]], sprintf("%s/%s.csv", mydir, names(tablelist)[i]))
  }
}
