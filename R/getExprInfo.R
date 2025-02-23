################################################################################
# Get gene average expression, detection rate, and differential expression info
# Input:
# seuratobj: Seurat object of the input scRNA-seq data, with sample conditions and cell-types/clusters in the columns "cond" and "cl" of the meta.data.
# cls: cell-types/clusters.
# conds: sample conditions.
# condpairs: a matrix with two columns "cond1" and "cond2", each row is ordered, represents a comparison between two sample conditions of interest, cond1 vs cond2 in differential expression analysis.
# testmethod: the test.use option for differential expression analysis in the function "FindMarkers".
# Output:
# A list exprinfo:
# exprinfo$avgexprall: the average expression of all genes in all cells;
# exprinfo$avgexprlist: list of average expression matrix (genes in rows, sample conditions in columns) of each cell-type/cluser in cls;
# exprinfo$pctexprlist: list of gene detection rate matrix (genes in rows, sample conditions in columns) of each cell-type/cluser in cls;
# exprinfo$DElist: list of differential expression matrix (genes in rows) of each cell-type/cluser in cls, each sample condition pair in condpairs.
################################################################################

#' Get gene average expression, detection rate, and differential expression info
#'
#' @param seuratobj Seurat object of the input scRNA-seq data, with sample conditions and cell-types/clusters in the columns "cond" and "cl" of the meta.data.
#' @param cls cell-types/clusters.
#' @param conds sample conditions.
#' @param condpairs a matrix with two columns "cond1" and "cond2", each row is ordered, represents a comparison between two sample conditions of interest, cond1 vs cond2 in differential expression analysis.
#' @param testmethod the test.use option for differential expression analysis in the function "FindMarkers".
#'
#' @return
#' A list exprinfo:
#' exprinfo$avgexprall: the average expression of all genes in all cells
#' exprinfo$avgexprlist: list of average expression matrix (genes in rows, sample conditions in columns) of each cell-type/cluser in cls
#' exprinfo$pctexprlist: list of gene detection rate matrix (genes in rows, sample conditions in columns) of each cell-type/cluser in cls
#' exprinfo$DElist: list of differential expression matrix (genes in rows) of each cell-type/cluser in cls, each sample condition pair in condpairs
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


getExprInfo <- function(seuratobj, cls, conds, condpairs, testmethod = "wilcox") {
  avgall = log1p(mean(expm1(as.matrix(seuratobj[['RNA']]$data))))
  avgexprlist = list()
  pctexprlist = list()
  DElist = list()
  for (k in 1:length(cls)) {
    print(sprintf("k = %s (in 1:%s), cl = %s", k, length(cls), cls[k]))
    myobj = subset(seuratobj, subset = cl == cls[k])
    Idents(myobj) = 'cond'
    avg_expr_obj = AverageExpression(myobj, group.by = 'cond', return.seurat = TRUE)
    avg_expr = avg_expr_obj[['RNA']]$data[,conds,drop=FALSE]
    normdata = myobj[['RNA']]$data
    pct_expr = matrix(0, nrow = nrow(normdata), ncol = length(conds))
    rownames(pct_expr) = rownames(normdata)
    colnames(pct_expr) = conds
    for (i in 1:length(conds)) {
      subnormdata = normdata[,WhichCells(myobj, idents = conds[i]),drop=FALSE]
      subnormdata = as.matrix(subnormdata)
      if (ncol(subnormdata) > 0) {
        pct_expr[,i] = rowSums(subnormdata > 0)/ncol(subnormdata)
      }
    }
    DEinfolist = list()
    for (i in 1:nrow(condpairs)) {
      DEinfo = FindMarkers(myobj, ident.1 = condpairs[i,"cond1"], ident.2 = condpairs[i,"cond2"],
                           only.pos = FALSE, min.pct = 0, logfc.threshold = 0, test.use = testmethod)
      # use "wilcox_limma" for limma implementation of the Wilcoxon Rank Sum test
      DEinfolist[[i]] = DEinfo
    }
    names(DEinfolist) = sprintf("%s_vs_%s", condpairs[,"cond1"], condpairs[,"cond2"])
    avgexprlist[[k]] = avg_expr
    pctexprlist[[k]] = pct_expr
    DElist[[k]] = DEinfolist
  }
  names(avgexprlist) = cls
  names(pctexprlist) = cls
  names(DElist) = cls
  exprinfo = list(avgall, avgexprlist, pctexprlist, DElist)
  names(exprinfo) = c("avgexprall", "avgexprlist", "pctexprlist", "DElist")

  return(exprinfo)
}

