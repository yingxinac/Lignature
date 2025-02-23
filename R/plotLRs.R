################################################################################
# Plot the average expression of ligands in sender cells, receptors in specified receiver cells, or the ligand scores calculated by Lignature signatures for identified LR interactions.
# Input:
# LRinfolist: output of "getLRinfosummary".
# LRs: vector of LR interactions "L_Lgene_R_Rgene" to plot.
# cls_from: sender cell-types/clusters.
# cl_to: specifies receiver cell-type/cluster.
# condpair: sample condition pair to plot, in the form of "cond1_vs_cond2".
# Lscore.by: how to score each ligand with respect to the scores of the set of signatures corresponding to it, one of "Lscore_mean", "Lscore_median", "Lscore_max", "Lscore_min", "Lscore_max_abs", and "Lscore_signed_maxabs".
# order.by: order the LR interactions by which column of the matrices in LRinfolist.
# order.decreasing: TRUE or FALSE, whether sort the LR interactions in a descending order or not, according to the values in the column specified by "order.by".
# plotwhich: plot average expression of ligand ("Lavgexpr") or receptor ("Ravgexpr") genes, or the ligand scores specified by "Lscore.by".
# col, scale, cluster_cols, cluster_rows, show_rownames, treeheight_row, treeheight_col, annotation_legend, border_color, fontsize_row, fontsize_col, fontsize: same parameters as in the function "pheatmap" in the "pheatmap" package.
# Output:
# Heatmap of the average expression of ligands in sender cells in each condition, average expression of receptors in specified receiver cells in each condition,
# or the corresponding ligand scores calculated by Lignature signatures, for identified LR interactions.
################################################################################

#' Plot the average expression of ligands in sender cells, receptors in specified receiver cells, or the ligand scores calculated by Lignature signatures for identified LR interactions
#'
#' @param LRinfolist output of "getLRinfosummary".
#' @param LRs vector of LR interactions "L_Lgene_R_Rgene" to plot.
#' @param cls_from sender cell-types/clusters.
#' @param cl_to specifies receiver cell-type/cluster.
#' @param condpair sample condition pair to plot, in the form of "cond1_vs_cond2".
#' @param Lscore.by how to score each ligand with respect to the scores of the set of signatures corresponding to it, one of "Lscore_mean", "Lscore_median", "Lscore_max", "Lscore_min", "Lscore_max_abs", and "Lscore_signed_maxabs".
#' @param order.by order the LR interactions by which column of the matrices in LRinfolist.
#' @param order.decreasing TRUE or FALSE, whether sort the LR interactions in a descending order or not, according to the values in the column specified by "order.by".
#' @param plotwhich plot average expression of ligand ("Lavgexpr") or receptor ("Ravgexpr") genes, or the ligand scores specified by "Lscore.by".
#' @param col same as in the function "pheatmap" in the "pheatmap" package.
#' @param scale same as in the function "pheatmap" in the "pheatmap" package.
#' @param cluster_cols same as in the function "pheatmap" in the "pheatmap" package.
#' @param cluster_rows same as in the function "pheatmap" in the "pheatmap" package.
#' @param show_rownames same as in the function "pheatmap" in the "pheatmap" package.
#' @param treeheight_row same as in the function "pheatmap" in the "pheatmap" package.
#' @param treeheight_col same as in the function "pheatmap" in the "pheatmap" package.
#' @param annotation_legend same as in the function "pheatmap" in the "pheatmap" package.
#' @param border_color same as in the function "pheatmap" in the "pheatmap" package.
#' @param fontsize_row same as in the function "pheatmap" in the "pheatmap" package.
#' @param fontsize_col same as in the function "pheatmap" in the "pheatmap" package.
#' @param fontsize same as in the function "pheatmap" in the "pheatmap" package.
#'
#' @return
#' Heatmap of the average expression of ligands in sender cells in each condition, average expression of receptors in specified receiver cells in each condition,
#' or the corresponding ligand scores calculated by Lignature signatures, for identified LR interactions
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


plotLRs <- function(LRinfolist, LRs, cls_from, cl_to, condpair, Lscore.by = "Lscore_max", order.by, order.decreasing = T, plotwhich = "LScore",
                    col, scale = 'none', cluster_cols = F, cluster_rows = F, show_rownames = T, treeheight_row = 0, treeheight_col = 0,
                    annotation_legend = TRUE, border_color = FALSE, fontsize_row = 13, fontsize_col = 13, fontsize = 13) {

  condsinpair = str_split(condpair, "_vs_")[[1]]
  order.vec = LRinfolist[[1]][,order.by]
  names(order.vec) = rownames(LRinfolist[[1]])
  order.vec = order.vec[LRs]
  order.vec = sort(order.vec, decreasing = order.decreasing)
  LRs = names(order.vec)
  clanno = cbind(rep(cls_from, each=length(condsinpair)), rep(condsinpair, length(cls_from)))
  rownames(clanno) = sprintf("%s_%s", clanno[,1], clanno[,2])
  colnames(clanno) = c("cl", "cond")
  clanno = data.frame(clanno)
  clanno$cond = factor(clanno$cond, levels = condsinpair)
  clanno$cl = factor(clanno$cl, levels = cls_from)
  if (plotwhich == "Lavgexpr") {
    mycols = vector()
    for (i in 1:length(cls_from)) {
      mycols = c(mycols, sprintf("%s_%s", cls_from[i], condsinpair))
    }
    mymatrix = matrix(0, nrow = length(LRs), ncol = length(mycols))
    rownames(mymatrix) = LRs
    colnames(mymatrix) = mycols
    mymatrix = data.frame(mymatrix)
    for (i in 1:length(cls_from)) {
      clpair = sprintf("%s_to_%s", cls_from[i], cl_to)
      mymatrix[,sprintf("%s_%s", cls_from[i], condsinpair)] = LRinfolist[[clpair]][LRs,sprintf("Lavg_%s", condsinpair)]
    }
    mymatrix = as.matrix(mymatrix)
    mode(mymatrix) = "numeric"
    rownames(mymatrix) = sprintf("%s_%s", str_split(LRs, "_", simplify = T)[,2], str_split(LRs, "_", simplify = T)[,3])
    pheatmap(mymatrix, scale = scale, cluster_cols = cluster_cols, cluster_rows = cluster_rows, show_rownames = show_rownames, color = col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             annotation_col = clanno,
             annotation_legend = annotation_legend,
             border_color = border_color,
             fontsize_row = fontsize_row,
             fontsize_col = fontsize_col,
             fontsize = fontsize)
  } else if (plotwhich == "Ravgexpr") {
    mycols = sprintf("%s_%s", cl_to, condsinpair)
    mymatrix = matrix(0, nrow = length(LRs), ncol = length(mycols))
    rownames(mymatrix) = LRs
    colnames(mymatrix) = mycols
    mymatrix = data.frame(mymatrix)
    clpair = sprintf("%s_to_%s", cls_from[1], cl_to)
    mymatrix[,sprintf("%s_%s", cl_to, condsinpair)] = LRinfolist[[clpair]][LRs,sprintf("Ravg_%s", condsinpair)]
    mymatrix = as.matrix(mymatrix)
    mode(mymatrix) = "numeric"
    rownames(mymatrix) = sprintf("%s_%s", str_split(LRs, "_", simplify = T)[,2], str_split(LRs, "_", simplify = T)[,3])
    pheatmap(mymatrix, scale = scale, cluster_cols = cluster_cols, cluster_rows = cluster_rows, show_rownames = show_rownames, color = col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             annotation_col = clanno,
             annotation_legend = annotation_legend,
             border_color = border_color,
             fontsize_row = fontsize_row,
             fontsize_col = fontsize_col,
             fontsize = fontsize)
  } else if (plotwhich == "LScore") {
    mymatrix = matrix(0, nrow = length(LRs), ncol = 1)
    rownames(mymatrix) = LRs
    colnames(mymatrix) = "LScore"
    mymatrix = data.frame(mymatrix)
    clpair = sprintf("%s_to_%s", cls_from[1], cl_to)
    mymatrix[,"LScore"] = LRinfolist[[clpair]][LRs,Lscore.by]
    mymatrix = as.matrix(mymatrix)
    mode(mymatrix) = "numeric"
    rownames(mymatrix) = sprintf("%s_%s", str_split(LRs, "_", simplify = T)[,2], str_split(LRs, "_", simplify = T)[,3])
    pheatmap(mymatrix, scale = scale, cluster_cols = cluster_cols, cluster_rows = cluster_rows, show_rownames = show_rownames, color = col,
             treeheight_row = treeheight_row,
             treeheight_col = treeheight_col,
             annotation_legend = annotation_legend,
             border_color = border_color,
             fontsize_row = fontsize_row,
             fontsize_col = fontsize_col,
             fontsize = fontsize)
  }
}
