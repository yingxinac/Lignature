
#' For ligands in identified LR interactions, create scatterplot of the ligand scores and the maximum avgexpr levels acorss all sender cells and conditions.
#'
#' @param LRinfolist output of "getLRinfosummary".
#' @param LRs vector of identified LR interactions "L_Lgene_R_Rgene".
#' @param cls_from sender cell-types/clusters.
#' @param cl_to specifies receiver cell-type/cluster.
#' @param condpair sample condition pair for the plot, in the form of "cond1_vs_cond2".
#' @param Lscore.by how to score each ligand with respect to the scores of the set of signatures corresponding to it, one of "Lscore_mean", "Lscore_median", "Lscore_max", "Lscore_min", "Lscore_max_abs", and "Lscore_signed_maxabs".
#' @param point.size point size.
#'
#' @return
#' Scatterplot of the ligand scores and the maximum avgexpr levels acorss all sender cells and conditions for ligands in identified LR interactions.
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


getScatter <- function(LRinfolist, LRs, cls_from, cl_to, condpair, Lscore.by = "Lscore_max", point.size = 3) {
  
  condsinpair = str_split(condpair, "_vs_")[[1]]
  
  clpair = sprintf("%s_to_%s", cls_from[1], cl_to)
  LScore = LRinfolist[[clpair]][LRs,Lscore.by]
  names(LScore) = LRs
  mode(LScore) = "numeric"
  
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
  maxLexpr = apply(mymatrix, 1, max)
  maxLexprfrom = colnames(mymatrix)[apply(mymatrix, 1, which.max)]
  names(maxLexpr) = rownames(mymatrix)
  names(maxLexprfrom) = rownames(mymatrix)
  
  myscattermx = cbind.data.frame(LScore, maxLexpr, maxLexprfrom)
  myscattermx = myscattermx[!is.na(myscattermx[,1]),,drop=F]
  myscattermx = myscattermx[!is.na(myscattermx[,1]),,drop=F]
  myscattermx$L = str_split(rownames(myscattermx), "_", simplify = T)[,2]
  myscattermx = unique(myscattermx)
  rownames(myscattermx) = myscattermx$L
  plotdf = myscattermx
  numcols = length(unique(plotdf$maxLexprfrom))
  mycols = distinctColorPalette(numcols)
  ggplot(plotdf, aes(x=LScore, y= maxLexpr, color= maxLexprfrom)) + 
    geom_point(size= point.size) +
    theme_ipsum() +
    geom_label_repel(data=plotdf,
                     aes(label = L),
                     label.padding = 0.1,
                     box.padding   = 0.1, 
                     point.padding = 0.3,
                     segment.color = 'grey50',show.legend = FALSE, size=6) +
    theme(axis.title.y = element_text(size = rel(2), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(2), angle = 00)) +
    scale_color_manual(values= mycols) +
    guides(col=guide_legend(title="Ligand maxExpr")) +
    labs(x="Lignature score", y="Ligand expression (max)", title="") +
    theme(axis.text = element_text(size = 15)) +
    theme(legend.text = element_text(size = 15)) +
    theme(axis.title = element_text(size = 10)) +
    theme(plot.title = element_text(size = 15)) +
    theme(legend.title = element_text(size = 15), 
          axis.line = element_line())
  
}
