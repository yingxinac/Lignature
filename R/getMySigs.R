################################################################################
# Get signatures of input data.
# Input:
# DElist: output of function "getExpInfo", exprinfo$DElist, a list of differential expression matrix (genes in rows) of each cell-type/cluser in cls, each sample condition pair.
# cls: cell-types/clusters.
# condpairs: a matrix with two columns "cond1" and "cond2", each row is ordered, represents a comparison between two sample conditions of interest, cond1 vs cond2 in differential expression analysis.
# keggnodes: list of the gene-nodes (represented by gene symbols) of KEGG pathways.
# pctcut: gene detection rate cutoff-value.
# belowpctcut: specify how to deal with avg_log2FC values of genes with detection rate below pctcut, one of "as.is", "lowpct.as.zero", or "lowpct.remove".
# testdegpadjcut: gene p_val_adj (from differential expression analysis) cutoff-value, avg_log2FC of genes with p_val_adj above this value is zeroed out.
# Output:
# A list of signatures for each cell-type/cluser in cls, each sample condition pair in condpairs.
################################################################################

#' Get signatures of input data
#'
#' @param DElist output of function "getExpInfo", exprinfo$DElist, a list of differential expression matrix (genes in rows) of each cell-type/cluser in cls, each sample condition pair.
#' @param cls cell-types/clusters.
#' @param condpairs a matrix with two columns "cond1" and "cond2", each row is ordered, represents a comparison between two sample conditions of interest, cond1 vs cond2 in differential expression analysis.
#' @param keggnodes list of the gene-nodes (represented by gene symbols) of KEGG pathways.
#' @param pctcut gene detection rate cutoff-value.
#' @param belowpctcut specify how to deal with avg_log2FC values of genes with detection rate below pctcut, one of "as.is", "lowpct.as.zero", or "lowpct.remove".
#' @param testdegpadjcut gene p_val_adj (from differential expression analysis) cutoff-value, avg_log2FC of genes with p_val_adj above this value is zeroed out.
#'
#' @return
#' a list of signatures for each cell-type/cluser in cls, each sample condition pair in condpairs
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


getMySigs <- function(DElist, cls, condpairs, keggnodes,
                      pctcut = 0.1, belowpctcut = "as.is", testdegpadjcut = 1) {

  testsiglist = list()
  for (k in 1:length(cls)) {
    print(sprintf("k = %s (in 1:%s), cl = %s", k, length(cls), cls[k]))
    cl = cls[k]
    testsiglistk = list()
    for (i in 1:nrow(condpairs)) {
      condpair = sprintf("%s_vs_%s", condpairs[i,"cond1"], condpairs[i,"cond2"])
      DEGinfo = DElist[[cl]][[condpair]]
      DEGinfo = DEGinfo[!is.na(rownames(DEGinfo)),,drop=F]
      DEGinfo = DEGinfo[rownames(DEGinfo) != "",,drop=F]
      DEGinfo = DEGinfo[rownames(DEGinfo) != " ",,drop=F]
      DEGinfo = DEGinfo[rowSums(is.na(DEGinfo)) == 0,,drop=F]
      DEGinfo = DEGinfo[rowSums(DEGinfo == "") ==0,,drop=F]
      DEGinfo = DEGinfo[rowSums(DEGinfo == " ") ==0,,drop=F]
      ### genelfc
      testlfc = DEGinfo$avg_log2FC
      names(testlfc) = rownames(DEGinfo)
      testpadj = DEGinfo$p_val_adj
      names(testpadj) = rownames(DEGinfo)
      testmaxpct = apply(DEGinfo[,c("pct.1", "pct.2")], 1, max)
      names(testmaxpct) = rownames(DEGinfo)
      if (belowpctcut == "lowpct.remove") {
        testlfc = testlfc[names(testmaxpct)[testmaxpct >= pctcut]]
        testpadj = testpadj[names(testmaxpct)[testmaxpct >= pctcut]]
      } else if (belowpctcut == "lowpct.as.zero") {
        testlfc[names(testmaxpct)[testmaxpct < pctcut]] = 0
      }
      testlfc[names(testpadj)[testpadj > testdegpadjcut]] = 0
      ### fgsea
      if (sum(testlfc != 0) >= 3) {
        myvec = testlfc[order(testlfc, decreasing = T)]
        myvec[myvec > 1000] = 1001
        myvec[myvec < -1000] = -1001
        fgseaRes = fgsea(pathways = keggnodes,
                         stats = myvec,
                         minSize=0,
                         maxSize=1000,
                         eps = 0.0)
        fgseaRes = data.frame(fgseaRes)
        rownames(fgseaRes) = fgseaRes$pathway
      } else {
        fgseaRes = matrix(0, nrow = 2, ncol = 2)
        colnames(fgseaRes) = c("NES", "padj")
        rownames(fgseaRes) = c("NA1", "NA2")
        fgseaRes = data.frame(fgseaRes)
      }
      fgseaNES = fgseaRes$NES
      names(fgseaNES) = rownames(fgseaRes)
      fgseapadj = fgseaRes$padj
      names(fgseapadj) = rownames(fgseaRes)
      ### put together
      testsiglistki = list(testlfc, testpadj, fgseaNES, fgseapadj)
      names(testsiglistki) = c("lfc", "padj", "NES", "fgseapadj")
      testsiglistk[[i]] = testsiglistki
    }
    names(testsiglistk) = sprintf("%s_vs_%s", condpairs[,"cond1"], condpairs[,"cond2"])
    testsiglist[[k]] = testsiglistk
  }
  names(testsiglist) = cls

  return(testsiglist)
}





