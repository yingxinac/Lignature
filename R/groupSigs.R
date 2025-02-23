################################################################################
# Group the signatures by treatment ligands.
# Input:
# siglist: siglist of Lignature
# sigmeta: sigmeta of Lignature
# Output:
# A list siggroups, with siggroups$sigLs = the vector of ligands in Lignature, siggroups$siggroups = the list of signatureid-vectors of each ligand
################################################################################

#' Group the signatures by treatment ligands
#'
#' @param siglist siglist of Lignature
#' @param sigmeta sigmeta of Lignature
#'
#' @return
#' a list siggroups:
#' siggroups$sigLs: the vector of ligands in Lignature
#' siggroups$siggroups: the list of signatureid-vectors of each ligand
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


groupSigs <- function(siglist, sigmeta) {

  sigmetaLs = unique(sigmeta$Ligand)
  sigLs = vector()
  for (i in 1:length(sigmetaLs)) {
    sigLs = union(sigLs, str_split(sigmetaLs[i], "and")[[1]])
  }

  siggroups = vector(mode = "list", length = length(sigLs))
  names(siggroups) = sigLs
  for (i in 1:length(siglist)) {
    sigLi = str_split(siglist[[i]]$sigL, "and")[[1]]
    for (j in 1:length(sigLi)) {
      siggroups[[sigLi[j]]] = union(siggroups[[sigLi[j]]], names(siglist)[i])
    }
  }

  siggroups = list(sigLs, siggroups)
  names(siggroups) = c("sigLs", "siggroups")
  return(siggroups)
}




