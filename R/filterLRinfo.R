################################################################################
# Filter LR info (expression, interaction score, and ligand score)
# Input:
# LRinfolist: output of "getLRinfosummary".
# cut.bounds: a list of lower- and upper-bounds of selected columns of data.frames in LRinfolist.
# Output:
# Filtered LRinfolist.
################################################################################

#' Filter LR info (expression, interaction score, and ligand score)
#'
#' @param LRinfolist output of "getLRinfosummary".
#' @param cut.bounds a list of lower- and upper-bounds of selected columns of data.frames in LRinfolist.
#'
#' @return
#' Filtered LRinfolist
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


filterLRinfo <- function(LRinfolist, cut.bounds) {
  LRinfolist_fil = list()
  lrs_fil_union = vector()
  for (k in 1:length(LRinfolist)) {
    mylrs = rownames(LRinfolist[[k]])
    for (m in 1:length(cut.bounds)) {
      myentry = names(cut.bounds)[m]
      myvec = LRinfolist[[k]][,myentry]
      names(myvec) = rownames(LRinfolist[[k]])
      myvec = myvec[!is.na(myvec)]
      mylowerbound = cut.bounds[[m]][1]
      myupperbound = cut.bounds[[m]][2]
      if (mylowerbound <= myupperbound) {
        mylrscutm = names(myvec)[myvec >= mylowerbound & myvec <= myupperbound]
      } else {
        mylrscutm = names(myvec)[myvec >= mylowerbound | myvec <= myupperbound]
      }
      mylrs = intersect(mylrs, mylrscutm)
    }
    LRinfolist_fil[[k]] = LRinfolist[[k]][mylrs,,drop=F]
    lrs_fil_union = union(lrs_fil_union, mylrs)
  }
  names(LRinfolist_fil) = names(LRinfolist)
  LRinfo_fil = list(LRinfolist_fil, lrs_fil_union)
  names(LRinfo_fil) = c("LRinfolist_fil", "lrs_fil_union")

  return(LRinfo_fil)
}
