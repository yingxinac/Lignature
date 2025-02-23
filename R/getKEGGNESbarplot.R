################################################################################
# get KEGG pathway NES barplot.
# Input:
# mysig: signature list of input data.
# padj.cut: specifies cutoff value on the adjusted p.values of the pathways from GSEA.
# num.decimal: number of decimal places to keep for the NES values in the plot.
# break.bound, break.step: breaks on the x-axis is set to seq(break.bound*(-1), break.bound, by = break.step).
# Output:
# barplot of the NES values of enriched KEGG pathways (with adjusted p.value < padj.cut) resulted from GSEA.
################################################################################

#' get KEGG pathway NES barplot
#'
#' @param mysig signature list of input data.
#' @param padj.cut specifies cutoff value on the adjusted p.values of the pathways from GSEA.
#' @param num.decimal number of decimal places to keep for the NES values in the plot.
#' @param break.bound breaks on the x-axis is set to seq(break.bound*(-1), break.bound, by = break.step).
#' @param break.step breaks on the x-axis is set to seq(break.bound*(-1), break.bound, by = break.step).
#'
#' @return
#' barplot of the NES values of enriched KEGG pathways (with adjusted p.value < padj.cut) resulted from GSEA
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


getKEGGNESbarplot <- function(mysig, padj.cut = 0.25, num.decimal = 2, break.bound = 2, break.step = 1) {
  mykeggsig = cbind.data.frame(mysig$NES, mysig$fgseapadj, names(mysig$NES))
  rownames(mykeggsig) = names(mysig$NES)
  colnames(mykeggsig) = c("NES", "padj", "pathway")
  subkeggsig = mykeggsig[mykeggsig$padj < padj.cut,,drop=F]
  subkeggsig = subkeggsig[order(subkeggsig$padj, decreasing = F),,drop=F]
  ggplot(subkeggsig, aes(x = reorder(pathway, NES), y = NES)) +
    geom_bar(stat = "identity",
             show.legend = FALSE,
             aes(fill = NES),
             color = "white") +
    geom_text(aes(label = round(NES, num.decimal),
                  hjust = ifelse(NES < 0, 1.5, -1),
                  vjust = 0.5),
              size = 3) +
    xlab("pathway") +
    ylab("NES") +
    scale_fill_gradient2(low = "pink",
                         mid = "aliceblue",
                         high = "lightblue") +
    coord_flip() +
    scale_y_continuous(breaks= seq(break.bound*(-1), break.bound, by = break.step),
                       limits = c(min(subkeggsig$NES) - 0.5,
                                  max(subkeggsig$NES) + 0.5)) +
    theme(panel.background = element_blank(),
          axis.line = element_line())
}
