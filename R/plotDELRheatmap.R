################################################################################
# Plot expression of differentially expressed ligand or receptor genes in each cell-type/cluster, each condition, in heatmap.
# Input:
# exprinfo: output of function "getExpInfo".
# siggroups: list of signatureid-vectors for each ligand in Lignature.
# lr_network: ligand_receptor network matrix with columns "L", "Lgene", "R", and "Rgene".
# whichLs: specifies the set of ligand genes to examine. "in.lrnetwork" (all ligand genes in lr_network) or "in.Lignature.only" (ligand genes in lr_network that are included in Lignature).
# cls: cell-types/clusters to plot.
# condpair: sample condition pair to plot, in the form of "cond1_vs_cond2".
# LorR: plot expression of differentially expressed ligand (L) or receptor (R) genes.
# pct.cut: cutoff value on the maximum detection rate across the sample conditions under comparison for each gene, for the identification of differentially expressed genes.
# padj.cut: cutoff value on the adjusted p.values of genes in the differential expression analysis, for the identification of differentially expressed genes.
# lfc.cut: cutoff value on the |avg_log2FC| values of genes in the differential expression analysis, for the identification of differentially expressed genes.
# col, scale, cluster_cols, cluster_rows, show_rownames, treeheight_row, treeheight_col, annotation_legend, border_color, fontsize_row, fontsize_col, fontsize: same parameters as in the function "pheatmap" in the "pheatmap" package.
# Output:
# Heatmap of the average expression values of differentially expressed ligand or receptor genes in each cell-type/cluster in each condition.
################################################################################

#' Plot expression of differentially expressed ligand or receptor genes in each cell-type/cluster, each condition, in heatmap
#'
#' @param exprinfo output of function "getExpInfo".
#' @param siggroups list of signatureid-vectors for each ligand in Lignature.
#' @param lr_network ligand_receptor network matrix with columns "L", "Lgene", "R", and "Rgene".
#' @param whichLs specifies the set of ligand genes to examine. "in.lrnetwork" (all ligand genes in lr_network) or "in.Lignature.only" (ligand genes in lr_network that are included in Lignature).
#' @param cls cell-types/clusters to plot.
#' @param condpair sample condition pair to plot, in the form of "cond1_vs_cond2".
#' @param LorR plot expression of differentially expressed ligand (L) or receptor (R) genes.
#' @param pct.cut cutoff value on the maximum detection rate across the sample conditions under comparison for each gene, for the identification of differentially expressed genes.
#' @param padj.cut cutoff value on the adjusted p.values of genes in the differential expression analysis, for the identification of differentially expressed genes.
#' @param lfc.cut cutoff value on the |avg_log2FC| values of genes in the differential expression analysis, for the identification of differentially expressed genes.
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
#' heatmap of the average expression values of differentially expressed ligand or receptor genes in each cell-type/cluster in each condition
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


plotDELRheatmap <- function(exprinfo, siggroups, lr_network, whichLs = "in.Lignature.only", cls, condpair, LorR = "L",
                            pct.cut = 0.1, padj.cut = 0.05, lfc.cut = 0.1,
                            col, scale = 'row', cluster_cols = F, cluster_rows = T, show_rownames = T, treeheight_row = 0, treeheight_col = 0,
                            annotation_legend = TRUE, border_color = FALSE, fontsize_row = 13, fontsize_col = 13, fontsize = 13) {

  if (whichLs == "in.Lignature.only") {
    myLgenes = unique(lr_network[lr_network[,"L"] %in% names(siggroups),"Lgene"])
    myRgenes = unique(lr_network[lr_network[,"L"] %in% names(siggroups),"Rgene"])
  } else if (whichLs == "in.lrnetwork") {
    myLgenes = unique(lr_network[,"Lgene"])
    myRgenes = unique(lr_network[,"Rgene"])
  }
  if (LorR == "L") {
    mygenes = myLgenes
  } else if (LorR == "R") {
    mygenes = myRgenes
  }
  candGs = vector()
  for (i in 1:length(mygenes)) {
    candGs = union(candGs, str_split(mygenes[i], "_")[[1]])
  }
  DEGs = vector()
  for (i in 1:length(cls)) {
    myDEmatrix = exprinfo$DElist[[cls[i]]][[condpair]]
    pctmax = apply(myDEmatrix[,c("pct.1", "pct.2"),drop=F], 1, max)
    DEGs = union(DEGs, rownames(myDEmatrix)[pctmax > pct.cut & myDEmatrix[,"p_val_adj"] < padj.cut & abs(myDEmatrix[,"avg_log2FC"]) > lfc.cut])
  }
  DEGs = intersect(DEGs, candGs)
  print(sprintf("num.DE%s = %s", LorR, length(DEGs)))
  mycols = vector()
  condsinpair = str_split(condpair, "_vs_")[[1]]
  for (i in 1:length(cls)) {
    mycols = c(mycols, sprintf("%s_%s", cls[i], condsinpair))
  }
  mymatrix = matrix(0, nrow = length(DEGs), ncol = length(mycols))
  rownames(mymatrix) = DEGs
  colnames(mymatrix) = mycols
  for (i in 1:length(cls)) {
    mymatrix[,sprintf("%s_%s", cls[i], condsinpair)] = exprinfo$avgexprlist[[cls[i]]][DEGs, condsinpair]
  }
  mode(mymatrix) = "numeric"
  mymatrix = data.frame(mymatrix)
  clanno = cbind(rep(cls, each=length(condsinpair)), rep(condsinpair, length(cls)))
  rownames(clanno) = colnames(mymatrix)
  colnames(clanno) = c("cl", "cond")
  clanno = data.frame(clanno)
  clanno$cond = factor(clanno$cond, levels = condsinpair)
  clanno$cl = factor(clanno$cl, levels = cls)
  pheatmap(mymatrix, scale = scale, cluster_cols = cluster_cols, cluster_rows = cluster_rows, show_rownames = show_rownames, color = col,
           treeheight_row = treeheight_row,
           treeheight_col = treeheight_col,
           annotation_col = clanno,
           annotation_legend = annotation_legend,
           border_color = border_color,
           fontsize_row = fontsize_row,
           fontsize_col = fontsize_col,
           fontsize = fontsize)
}

