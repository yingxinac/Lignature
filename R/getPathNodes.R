################################################################################
# Get gene-nodes in KEGG pathways.
# Input:
# species: one of the supported species.
# database: the name of the pathway database.
# numnode_min: minimum number of gene-nodes in the KEGG pathways to be collected (KEGG pathways with fewer number of gene nodes are removed).
# numnode_max: maximum number of gene-nodes in the KEGG pathways to be collected (KEGG pathways with larger number of gene nodes are removed).
# Output:
# A list of the gene-nodes (represented by gene symbols) of KEGG pathways with the number of included gene-nodes bounded by numnode_min and numnode_max.
################################################################################

#' Get gene-nodes in KEGG pathways
#'
#' @param species one of the supported species.
#' @param database the name of the pathway database.
#' @param numnode_min minimum number of gene-nodes in the KEGG pathways to be collected (KEGG pathways with fewer number of gene nodes are removed).
#' @param numnode_max maximum number of gene-nodes in the KEGG pathways to be collected (KEGG pathways with larger number of gene nodes are removed).
#'
#' @return
#' a list of the gene-nodes (represented by gene symbols) of KEGG pathways with the number of included gene-nodes bounded by numnode_min and numnode_max
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


getPathNodes <- function(species = "hsapiens", database = "kegg", numnode_min = 10, numnode_max = 500) {

  kegg = pathways(species, database)
  kegg = convertIdentifiers(kegg, "SYMBOL")
  kegg_node = list()
  numnode = vector()
  for (i in 1:length(kegg)) {
    kegg_node[[i]] = substr(nodes(kegg[[i]]), start = 8, stop = 1000)
    numnode[i] = length(kegg_node[[i]])
  }
  names(kegg_node) = names(kegg)
  names(numnode) = names(kegg)
  mykegg = kegg_node[numnode >= numnode_min & numnode <= numnode_max]

  return(mykegg)
}
