################################################################################
# Summarize LR info (expression, interaction score, and ligand score)
# Input:
# exprinfo: output of "getExprInfo".
# LRIscorelist: output of "getLRIScore".
# Lscoredf: output of "getLScore".
# cl_to: receiver cell-type/cluster.
# condpair: sample conditions under comparison, in the form of "cond1_vs_cond2".
# lr_network: a ligand_receptor network matrix with columns "L", "Lgene", "R", and "Rgene".
# Output:
# A list of data.frames of LR-info for each sender cell-type/cluster to the specified receiver cell-type/cluster.
################################################################################

#' get signed values
#'
#' @param x absolute value vector.
#' @param y sign vector.
#'
#' @return
#' signed values
#' @import stringr
#' @export
signedmaxabs <- function(x, y) {
  z = vector()
  for (i in 1:length(y)) {
    yi = str_split(y[i], "and")[[1]]
    if (is.na(yi)) {
      z[i] = NA
    } else if (yi == "-1") {
      z[i] = x[i]*(-1)
    } else {
      z[i] = x[i]
    }
  }
  return(z)
}

#' Summarize LR info (expression, interaction score, and ligand score)
#'
#' @param exprinfo output of "getExprInfo".
#' @param LRIscorelist output of "getLRIScore".
#' @param Lscoredf output of "getLScore".
#' @param lr_network a ligand_receptor network matrix with columns "L", "Lgene", "R", and "Rgene".
#' @param cl_to receiver cell-type/cluster.
#' @param condpair sample conditions under comparison, in the form of "cond1_vs_cond2".
#'
#' @return
#' A list of data.frames of LR-info for each sender cell-type/cluster to the specified receiver cell-type/cluster
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export

getLRinfosummary <- function(exprinfo, LRIscorelist, Lscoredf, lr_network, cl_to, condpair) {
  conds = colnames(exprinfo$avgexprlist[[1]])
  condsinpair = str_split(condpair, "_vs_")[[1]]
  otherLs = setdiff(rownames(Lscoredf), lr_network[,"L"])
  otherlrs = matrix(NA, nrow = length(otherLs), ncol = ncol(lr_network))
  colnames(otherlrs) = colnames(lr_network)
  otherlrs[,"L"] = otherLs
  alllrs = rbind.data.frame(lr_network, otherlrs)
  rownames(alllrs) = sprintf("%s_%s_%s_%s", alllrs[,"L"], alllrs[,"Lgene"], alllrs[,"R"], alllrs[,"Rgene"])
  LRinfolist = list()
  clpairids = vector()
  for (k in 1:length(exprinfo$avgexprlist)) {
    print(sprintf("%s of %s, %s", k, length(exprinfo$avgexprlist), names(exprinfo$avgexprlist)[k]))
    cl_from = names(exprinfo$avgexprlist)[k]
    clpairk = sprintf("%s_to_%s", cl_from, cl_to)
    clpairids[k] = clpairk
    LRinfodf = matrix(NA, nrow = nrow(alllrs), ncol = ncol(alllrs) + 8 + length(conds)*6 + 6)
    rownames(LRinfodf) = rownames(alllrs)
    colnames(LRinfodf) = c(colnames(alllrs), sprintf("Lscore_%s", colnames(Lscoredf)[1:8]),
                           sprintf("SLRI_%s", conds), sprintf("SLRIspec_%s", conds),
                           sprintf("Lpct_%s", conds), sprintf("Rpct_%s", conds), sprintf("Lavg_%s", conds), sprintf("Ravg_%s", conds),
                           sprintf("Lpval_%s", condpair), sprintf("Lavglog2FC_%s", condpair), sprintf("Lpadj_%s", condpair),
                           sprintf("Rpval_%s", condpair), sprintf("Ravglog2FC_%s", condpair), sprintf("Rpadj_%s", condpair))
    LRinfodf = data.frame(LRinfodf)
    LRinfodf[,colnames(alllrs)] = alllrs
    for (i in 1:nrow(LRinfodf)) {
      Li = LRinfodf[i,"L"]
      if (Li %in% rownames(Lscoredf)) {
        LRinfodf[i,sprintf("Lscore_%s", colnames(Lscoredf)[1:8])] = Lscoredf[Li,1:8]
      }
    }
    LRIscoredf = LRIscorelist$score[[clpairk]]
    LRIspecdf = LRIscorelist$spec[[clpairk]]
    rownames(LRIscoredf) = sprintf("%s_%s_%s_%s", LRIscoredf[,"L"], LRIscoredf[,"Lgene"], LRIscoredf[,"R"], LRIscoredf[,"Rgene"])
    rownames(LRIspecdf) = sprintf("%s_%s_%s_%s", LRIspecdf[,"L"], LRIspecdf[,"Lgene"], LRIspecdf[,"R"], LRIspecdf[,"Rgene"])
    whichrows = intersect(rownames(LRIscoredf), rownames(alllrs))
    LRinfodf[whichrows, sprintf("SLRI_%s", conds)] = LRIscoredf[whichrows, sprintf("SLRI_%s", conds), drop = F]
    LRinfodf[whichrows, sprintf("SLRIspec_%s", conds)] = LRIspecdf[whichrows, sprintf("SLRIspec_%s", conds), drop = F]
    for (i in 1:nrow(LRinfodf)) {
      Lgenesi = str_split(LRinfodf[i,"Lgene"], "_")[[1]]
      Rgenesi = str_split(LRinfodf[i,"Rgene"], "_")[[1]]
      if (length(intersect(Lgenesi, rownames(exprinfo$avgexprlist[[1]]))) == length(Lgenesi)) {
        for (j in 1:length(conds)) {
          LRinfodf[i,sprintf("Lpct_%s", conds[j])] = min(exprinfo$pctexprlist[[cl_from]][Lgenesi, conds[j]])
          LRinfodf[i,sprintf("Lavg_%s", conds[j])] = min(exprinfo$avgexprlist[[cl_from]][Lgenesi, conds[j]])
        }
        Lgenesivec = apply(exprinfo$avgexprlist[[cl_from]][Lgenesi, condsinpair, drop=F], 1, mean)
        names(Lgenesivec) = Lgenesi
        Lmin = names(Lgenesivec)[which.min(Lgenesivec)[1]]
        LRinfodf[i,c(sprintf("Lpval_%s", condpair), sprintf("Lavglog2FC_%s", condpair), sprintf("Lpadj_%s", condpair))] = exprinfo$DElist[[cl_from]][[condpair]][Lmin,c("p_val", "avg_log2FC", "p_val_adj")]
      }
      if (length(intersect(Rgenesi, rownames(exprinfo$avgexprlist[[1]]))) == length(Rgenesi)) {
        for (j in 1:length(conds)) {
          LRinfodf[i,sprintf("Rpct_%s", conds[j])] = min(exprinfo$pctexprlist[[cl_to]][Rgenesi, conds[j]])
          LRinfodf[i,sprintf("Ravg_%s", conds[j])] = min(exprinfo$avgexprlist[[cl_to]][Rgenesi, conds[j]])
        }
        Rgenesivec = apply(exprinfo$avgexprlist[[cl_to]][Rgenesi, condsinpair, drop=F], 1, mean)
        names(Rgenesivec) = Rgenesi
        Rmin = names(Rgenesivec)[which.min(Rgenesivec)[1]]
        LRinfodf[i,c(sprintf("Rpval_%s", condpair), sprintf("Ravglog2FC_%s", condpair), sprintf("Rpadj_%s", condpair))] = exprinfo$DElist[[cl_to]][[condpair]][Rmin,c("p_val", "avg_log2FC", "p_val_adj")]
      }
    }
    LRinfodf$Lscore_signed_maxabs = signedmaxabs(x=LRinfodf$Lscore_max_abs, y=LRinfodf$Lscore_sign_maxabs)
    LRinfodf$SLRI_condmax = apply(LRinfodf[,sprintf("SLRI_%s", conds), drop=F], 1, max)
    LRinfodf$SLRIspec_condmax = apply(LRinfodf[,sprintf("SLRIspec_%s", conds), drop=F], 1, max)
    LRinfodf$Lpct_condmax = apply(LRinfodf[,sprintf("Lpct_%s", conds), drop=F], 1, max)
    LRinfodf$Rpct_condmax = apply(LRinfodf[,sprintf("Rpct_%s", conds), drop=F], 1, max)
    LRinfodf$Lavg_condmax = apply(LRinfodf[,sprintf("Lavg_%s", conds), drop=F], 1, max)
    LRinfodf$Ravg_condmax = apply(LRinfodf[,sprintf("Ravg_%s", conds), drop=F], 1, max)
    LRinfolist[[k]] = LRinfodf
  }
  names(LRinfolist) = clpairids

  return(LRinfolist)
}







