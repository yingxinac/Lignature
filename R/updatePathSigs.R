################################################################################
# Update KEGG pathway signatures from fgsea.
# Input:
# siglist: the current siglist of Lignature.
# keggnodes: list of the gene-nodes (represented by gene symbols) of KEGG pathways.
# progress_interval: report progress very progress_interval runs in the for-loop.
# Output:
# siglist with updated pathway signatures.
################################################################################

#' Update KEGG pathway signatures from fgsea
#'
#' @param siglist the current siglist of Lignature.
#' @param keggnodes list of the gene-nodes (represented by gene symbols) of KEGG pathways.
#' @param progress_interval report progress very progress_interval runs in the for-loop.
#'
#' @return
#' siglist with updated pathway signatures
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


updatePathSigs <- function(siglist, keggnodes, progress_interval = 10) {

  for (i in 1:length(siglist)) {
    if (i %% progress_interval == 0) {
      print(sprintf("%s of %s", i, length(siglist)))
    }
    myvec = siglist[[i]][["lfc"]]
    myvec = myvec[!is.na(myvec)]
    myvec = myvec[myvec != ""]
    myvec = myvec[myvec != " "]
    myvec = myvec[!is.na(names(myvec))]
    myvec = myvec[names(myvec) != ""]
    myvec = myvec[names(myvec) != " "]
    myvec = myvec[order(myvec, decreasing = T)]
    myvec[myvec > 1000] = 1001
    myvec[myvec < -1000] = -1001
    fgseaRes = fgsea(pathways = keggnodes,
                     stats = myvec,
                     minSize=0,
                     maxSize=1000,
                     eps = 0.0)
    fgseaRes = data.frame(fgseaRes)
    rownames(fgseaRes) = fgseaRes$pathway
    myNES = fgseaRes$NES
    names(myNES) = rownames(fgseaRes)
    mypadj = fgseaRes$padj
    names(mypadj) = rownames(fgseaRes)
    siglist[[i]][["NES"]] = myNES
    siglist[[i]][["fgseapadj"]] = mypadj
  }

  return(siglist)
}
