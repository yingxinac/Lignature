
#' Get scatter plot of specified ligand scores and LR-interactoin scores (maximum from all sender cells) for LR pairs with these scores available.
#'
#' @param LRinfolist output of "getLRinfosummary".
#' @param cls_from sender cell-types/clusters.
#' @param cl_to specifies receiver cell-type/cluster.
#' @param Lscore.by how to score each ligand with respect to the scores of the set of signatures corresponding to it, one of "Lscore_mean", "Lscore_median", "Lscore_max", "Lscore_min", "Lscore_max_abs", and "Lscore_signed_maxabs".
#' @param SLRI.by specifies which LR-interaction score to plot, e.g. "SLRI_condmax" and "SLRIspec_condmax".
#' @param point.size point size.
#' @param Lscore.cut ligand score threshold.
#' @param SLRI.cut LR-interaction score threshold.
#'
#' @return
#' Scatterplot of specified ligand scores and LR-interactoin scores (maximum from all sender cells) for LR pairs with these scores available.
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


getScatterLRpairs <- function(LRinfolist, cls_from, cl_to, Lscore.by = "Lscore_max", SLRI.by = "SLRI_condmax",
                              point.size = 3, Lscore.cut = 0.1, SLRI.cut = 0.3) {

  LScore = LRinfolist[[1]][,Lscore.by]
  names(LScore) = rownames(LRinfolist[[1]])
  LScore = LScore[!is.na(LScore)]
  LScore = LScore[LScore > -2]
  mode(LScore) = "numeric"

  LRs = rownames(LRinfolist[[1]])
  SLRImatrix = matrix(0, nrow = length(LRs), ncol = length(cls_from))
  for (i in 1:length(cls_from)) {
    clpair = sprintf("%s_to_%s", cls_from[i], cl_to)
    SLRImatrix[,i] = LRinfolist[[clpair]][LRs,SLRI.by]
  }
  rownames(SLRImatrix) = LRs
  mode(SLRImatrix) = "numeric"
  SLRIvec = apply(SLRImatrix, 1, max)
  names(SLRIvec) = rownames(SLRImatrix)
  SLRIvec = SLRIvec[!is.na(SLRIvec)]
  SLRIvec = SLRIvec[SLRIvec > -2]

  LRs = intersect(names(LScore), names(SLRIvec))

  myplotmatrix = cbind.data.frame(LScore[LRs], SLRIvec[LRs])
  rownames(myplotmatrix) = LRs
  colnames(myplotmatrix) = c("Lscore", "SLRI")

  myplotmatrix$class = 0
  myplotmatrix[abs(myplotmatrix$Lscore) > Lscore.cut & myplotmatrix$SLRI > SLRI.cut,"class"] = "bothhigh"
  myplotmatrix[abs(myplotmatrix$Lscore) > Lscore.cut & myplotmatrix$SLRI < SLRI.cut,"class"] = "lowSLRI"
  myplotmatrix[abs(myplotmatrix$Lscore) < Lscore.cut & myplotmatrix$SLRI > SLRI.cut,"class"] = "lowLScore"
  myplotmatrix[abs(myplotmatrix$Lscore) < Lscore.cut & myplotmatrix$SLRI < SLRI.cut,"class"] = "bothlow"

  mycols = c("red", "yellowgreen", "skyblue", "lightgray")
  names(mycols) = c("lowSLRI", "bothhigh", "lowLScore", "bothlow")

  ggplot(myplotmatrix, aes(x= Lscore, y= SLRI, color= class)) +
    geom_point(size= point.size) +
    theme_ipsum() +
    geom_abline(intercept = SLRI.cut, slope = 0, color="gold",
                linetype="dashed", size=1) +
    geom_vline(xintercept = Lscore.cut, color="gold", linetype="dashed", size=1) +
    geom_vline(xintercept = (-1)*Lscore.cut, color="gold", linetype="dashed", size=1) +
    theme(axis.title.y = element_text(size = rel(2), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(2), angle = 00)) +
    scale_color_manual(values=mycols) +
    labs(x= "Lignature score", y= "SLRI (max)", title="") +
    theme(axis.text = element_text(size = 15)) +
    theme(legend.text = element_text(size = 15)) +
    theme(axis.title = element_text(size = 10)) +
    theme(plot.title = element_text(size = 15)) +
    theme(legend.title = element_text(size = 15),
          axis.line = element_line())

}




