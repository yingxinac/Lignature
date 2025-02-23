################################################################################
# Score each ligand in Lignature for input data.
# Input:
# siggroups: list of signatureid-vectors for each ligand in Lignature.
# lr_network: a ligand_receptor network matrix with ligand names in the column "L".
# sigscorelist: output from the function "getSigScore".
# pvalcut: p-value cutoff, only signatures with score-pvalues below this cutoff are used for scoring the ligands in Lignature.
# Output:
# A data.frame of ligand scores and the signature sources.
# For each ligand, each column scores the ligand by the set of signatures of it in Lignature (with the score-pvalues <= specified pvalcut):
# "mean": mean value of the signature scores;
# "median": median of the signature scores;
# "max": maximum of the signature scores;
# "min": minimum of the signature scores;
# "max_abs": maximum absolute value of the signature scores;
# additionally, the following columns are also included:
# "pos_maxabs": maximum of the positive signature scores;
# "neg_maxabs": minimum of the negative signature scores;
# "sign_maxabs": sign(s) of the signature scores with the maximum absolute value;
# "sig_***": signature-ids corresponding to the maximum, min, maximum absolute value, maximum positive value, and the minimum negative value of the signature scores.
# Remark: for ligands in Lignature with no signature score-pvalue <= pvalcut, the the first seven columns are set to -2, and the other six columns are set to "none";
# for ligands with no positve/negative signature scores with score-pvalue <= pvalcut, "pos_maxabs"/"neg_maxabs" is set -2, and "sig_pos_maxabs"/"sig_neg_maxabs" is set to "none";
# for ligands that are not included in Lignature, all columns are set to NA.
################################################################################

#' Score the ligands in Lignature
#'
#' @param siggroups list of signatureid-vectors for each ligand in Lignature.
#' @param lr_network a ligand_receptor network matrix with ligand names in the column "L".
#' @param sigscorelist output from the function "getSigScore".
#' @param pvalcut p-value cutoff, only signatures with score-pvalues below this cutoff are used for scoring the ligands in Lignature.
#'
#' @return
#' a data.frame of ligand scores and the signature sources
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


getLScore <- function(siggroups, lr_network, sigscorelist, pvalcut = 1) {
  allLs = union(names(siggroups), lr_network[,"L"])
  scorematrix = matrix(NA, nrow = length(allLs), ncol = 7)
  rownames(scorematrix) = allLs
  colnames(scorematrix) = c("mean", "median", "max", "min", "max_abs", "pos_maxabs", "neg_maxabs")
  scorematrix_sig = matrix(NA, nrow = length(allLs), ncol = 6)
  rownames(scorematrix_sig) = allLs
  colnames(scorematrix_sig) = c("sign_maxabs", "sig_max", "sig_min", "sig_max_abs", "sig_pos_maxabs", "sig_neg_maxabs")
  for (i in 1:length(siggroups)) {
    sigL = names(siggroups)[i]
    sigs = siggroups[[sigL]]
    sigscorevec = sigscorelist$score[sigs]
    sigpvalvec = sigscorelist$pvalue[sigs]
    sigs = intersect(names(sigscorevec)[!is.na(sigscorevec)], names(sigpvalvec)[!is.na(sigpvalvec)])
    sigscorevec = sigscorevec[sigs]
    sigpvalvec = sigpvalvec[sigs]
    sigs = names(sigpvalvec)[sigpvalvec <= pvalcut]
    if (length(sigs) > 0) {
      sigscorevec = sigscorevec[sigs]
      vecpos = sigscorevec[sigscorevec >= 0]
      vecneg = sigscorevec[sigscorevec <= 0]
      scorematrix[sigL,"mean"] = mean(sigscorevec)
      scorematrix[sigL,"median"] = median(sigscorevec)
      scorematrix[sigL,"max"] = max(sigscorevec)
      scorematrix[sigL,"min"] = min(sigscorevec)
      scorematrix[sigL,"max_abs"] = max(abs(sigscorevec))
      if (length(vecpos) > 0) {
        scorematrix[sigL,"pos_maxabs"] = max(vecpos)
      } else {
        scorematrix[sigL,"pos_maxabs"] = -2
      }
      if (length(vecneg) > 0) {
        scorematrix[sigL,"neg_maxabs"] = min(vecneg)
      } else {
        scorematrix[sigL,"neg_maxabs"] = -2
      }
      scorematrix_sig[sigL,"sign_maxabs"] = paste(unique(sign(sigscorevec[which.max(abs(sigscorevec))])), collapse = "and")
      scorematrix_sig[sigL,"sig_max"] = paste(names(sigscorevec)[which.max(sigscorevec)], collapse = "and")
      scorematrix_sig[sigL,"sig_min"] = paste(names(sigscorevec)[which.min(sigscorevec)], collapse = "and")
      scorematrix_sig[sigL,"sig_max_abs"] = paste(names(sigscorevec)[which.max(abs(sigscorevec))], collapse = "and")
      if (length(vecpos) > 0) {
        scorematrix_sig[sigL,"sig_pos_maxabs"] = paste(names(vecpos)[which.max(vecpos)], collapse = "and")
      } else {
        scorematrix_sig[sigL,"sig_pos_maxabs"] = "none"
      }
      if (length(vecneg) > 0) {
        scorematrix_sig[sigL,"sig_neg_maxabs"] = paste(names(vecneg)[which.min(vecneg)], collapse = "and")
      } else {
        scorematrix_sig[sigL,"sig_neg_maxabs"] = "none"
      }
    } else {
      scorematrix[sigL,"mean"] = -2
      scorematrix[sigL,"median"] = -2
      scorematrix[sigL,"max"] = -2
      scorematrix[sigL,"min"] = -2
      scorematrix[sigL,"max_abs"] = -2
      scorematrix[sigL,"pos_maxabs"] = -2
      scorematrix[sigL,"neg_maxabs"] = -2
      scorematrix_sig[sigL,"sign_maxabs"] = "none"
      scorematrix_sig[sigL,"sig_max"] = "none"
      scorematrix_sig[sigL,"sig_min"] = "none"
      scorematrix_sig[sigL,"sig_max_abs"] = "none"
      scorematrix_sig[sigL,"sig_pos_maxabs"] = "none"
      scorematrix_sig[sigL,"sig_neg_maxabs"] = "none"
    }
  }
  scorematrixdf = cbind.data.frame(scorematrix, scorematrix_sig)
  colnames(scorematrixdf) = c(colnames(scorematrix), colnames(scorematrix_sig))
  rownames(scorematrixdf) = allLs

  return(scorematrixdf)
}
