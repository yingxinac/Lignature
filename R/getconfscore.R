################################################################################
# Get confidence scores for specified ligand scores.
# Input:
# myscores: a named vector of ligand scores.
# siglist: siglist of Lignature.
# siggroups: list of signatureid-vectors for each ligand in Lignature.
# numsigs_sameL: maximum number of signatures used for each ligand when constructing the distribution of signature similarity scores between different ligands.
# whichsig: "lfc" or "NES".
# sigcut: cutoff value applied to the signatures.
# belowsigcut: "as.is", "remove", or "as.zero", specify how to deal with signature values below sigcut.
# padjcut: cutoff values applied to the adjusted p.values of each entry of the signatures.
# minsiglength: minimum length of signatures under comparison.
# method: similarity measurement metric, one of "cor.pearson", "cor.spearman", "cosine.similarity", "euclidean", "manhattan", "euclidean_mean", "manhattan_mean".
# Output:
# A vector of confidence scores of the input ligand scores.
################################################################################

#' Format x as percentage
#'
#' @param x a number to be formatted as percentage
#'
#' @return
#' x formatted as percentage
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table circlize tidyverse randomcoloR hrbrthemes ggrepel caret glmnet
#' @export

percent <- function(x) {
  paste0(formatC(x * 100, format = 'f', digits = 2), "%")
}


#' Construct distribution of signature similarity scores between different ligands
#'
#' @param siglist siglist of Lignature.
#' @param siggroups list of signatureid-vectors for each ligand in Lignature.
#' @param numsigs_sameL maximum number of signatures used for each ligand.
#' @param whichsig "lfc" or "NES".
#' @param sigcut cutoff value applied to the signatures.
#' @param belowsigcut "as.is", "remove", or "as.zero", specify how to deal with signature values below sigcut.
#' @param padjcut cutoff values applied to the adjusted p.values of each entry of the signatures.
#' @param minsiglength minimum length of signatures under comparison.
#' @param method similarity measurement metric, one of "cor.pearson", "cor.spearman", "cosine.similarity", "euclidean", "manhattan", "euclidean_mean", "manhattan_mean".
#'
#' @return
#' A data.frame of similarity scores between signatures of the ligands in Lignature
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table circlize tidyverse randomcoloR hrbrthemes ggrepel caret glmnet
#' @export

diffLpwscores <- function(siglist, siggroups, numsigs_sameL = 2, whichsig = "lfc",
                          sigcut = 0.1, belowsigcut = "as.is", padjcut = 1,
                          minsiglength = 10, method = "cor.pearson") {

  for (i in 1:length(siglist)) {
    siglist[[i]]$sigLsep = str_split(siglist[[i]]$sigL, "and")[[1]]
  }
  selectsigs = vector()
  for (i in 1:length(siggroups)) {
    if (length(siggroups[[i]]) > numsigs_sameL) {
      selectsigs = union(selectsigs, sample(siggroups[[i]], numsigs_sameL, replace = F))
    } else {
      selectsigs = union(selectsigs, siggroups[[i]])
    }
  }
  siglist = siglist[selectsigs]
  print("Constructing distribution of signature similarity scores between different ligands")
  sigpairs = t(combn(names(siglist),2))
  colnames(sigpairs) = c("siga", "sigb")
  rownames(sigpairs) = sprintf("%s_vs_%s", sigpairs[,1], sigpairs[,2])
  scorematrix = matrix(0, nrow = nrow(sigpairs), ncol = 6)
  rownames(scorematrix) = rownames(sigpairs)
  colnames(scorematrix) = c("sig_a", "sig_b", "Ligand_a", "Ligand_b", "score", "sameL")
  scorematrix = data.frame(scorematrix)
  for (i in 1:nrow(scorematrix)) {
    if (i %% 1000 == 0) {
      print(sprintf("progress = %s", percent(i/nrow(scorematrix))))
    }
    sigida = sigpairs[i,"siga"]
    sigidb = sigpairs[i,"sigb"]
    scorematrix[i,"sig_a"] = sigida
    scorematrix[i,"sig_b"] = sigidb
    scorematrix[i,"Ligand_a"] = siglist[[sigida]][["sigL"]]
    scorematrix[i,"Ligand_b"] = siglist[[sigidb]][["sigL"]]
    scorematrix[i,"sameL"] = (length(intersect(siglist[[sigida]][["sigLsep"]], siglist[[sigidb]][["sigLsep"]])) > 0) *1
    sigveca = siglist[[sigida]][[whichsig]]
    sigvecb = siglist[[sigidb]][[whichsig]]
    if (whichsig == "lfc") {
      sigpadjveca = siglist[[sigida]][["padj"]]
      sigpadjvecb = siglist[[sigidb]][["padj"]]
    } else if (whichsig == "NES") {
      sigpadjveca = siglist[[sigida]][["fgseapadj"]]
      sigpadjvecb = siglist[[sigidb]][["fgseapadj"]]
    }
    mode(sigveca) = "numeric"
    mode(sigpadjveca) = "numeric"
    mode(sigvecb) = "numeric"
    mode(sigpadjvecb) = "numeric"
    sigveca = sigveca[!is.na(sigveca)]
    sigpadjveca = sigpadjveca[!is.na(sigpadjveca)]
    sigvecb = sigvecb[!is.na(sigvecb)]
    sigpadjvecb = sigpadjvecb[!is.na(sigpadjvecb)]
    keepsigentries_a = names(sigveca)[abs(sigveca) >= sigcut]
    keepsigentries_b = names(sigvecb)[abs(sigvecb) >= sigcut]
    rmsigentries_a = names(sigveca)[abs(sigveca) < sigcut]
    rmsigentries_b = names(sigvecb)[abs(sigvecb) < sigcut]
    if (belowsigcut == "remove") {
      sigveca = sigveca[keepsigentries_a]
      sigpadjveca = sigpadjveca[keepsigentries_a]
      sigvecb = sigvecb[keepsigentries_b]
      sigpadjvecb = sigpadjvecb[keepsigentries_b]
    } else if (belowsigcut == "as.zero") {
      sigveca[rmsigentries_a] = 0
      sigvecb[rmsigentries_b] = 0
    }
    sigveca[names(sigpadjveca)[sigpadjveca > padjcut]] = 0
    sigvecb[names(sigpadjvecb)[sigpadjvecb > padjcut]] = 0
    commonentries = intersect(names(sigveca), names(sigvecb))
    if (length(commonentries) >= minsiglength) {
      sigveca = sigveca[commonentries]
      sigvecb = sigvecb[commonentries]
      if (method %in% c("cor.pearson", "cor.spearman")) {
        scorei = mymetric(a = sigveca, b = sigvecb, mymethod = method)$estimate
      } else if (method %in% c("cosine.similarity", "euclidean", "manhattan", "euclidean_mean", "manhattan_mean")) {
        scorei = mymetric(a = sigveca, b = sigvecb, mymethod = method)
      }
    } else {
      scorei = NA
    }
    scorematrix[i, "score"] = scorei
  }

  return(scorematrix)
}


#' Calculate confidence scores for a vector of ligand scores
#'
#' @param myscores a named vector of ligand scores.
#' @param siglist siglist of Lignature.
#' @param siggroups list of signatureid-vectors for each ligand in Lignature.
#' @param numsigs_sameL maximum number of signatures used for each ligand when constructing the distribution of signature similarity scores between different ligands.
#' @param whichsig "lfc" or "NES".
#' @param sigcut cutoff value applied to the signatures.
#' @param belowsigcut "as.is", "remove", or "as.zero", specify how to deal with signature values below sigcut.
#' @param padjcut cutoff values applied to the adjusted p.values of each entry of the signatures.
#' @param minsiglength minimum length of signatures under comparison.
#' @param method similarity measurement metric, one of "cor.pearson", "cor.spearman", "cosine.similarity", "euclidean", "manhattan", "euclidean_mean", "manhattan_mean".
#'
#' @return
#' A vector of confidence scores of the input ligand scores
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table circlize tidyverse randomcoloR hrbrthemes ggrepel caret glmnet
#' @export

getconfscore <- function(myscores, siglist, siggroups, numsigs_sameL = 2, whichsig = "lfc",
                         sigcut = 0.1, belowsigcut = "as.is", padjcut = 1,
                         minsiglength = 10, method = "cor.pearson") {
  scorematrix = diffLpwscores(siglist = siglist, siggroups = siggroups, numsigs_sameL = numsigs_sameL, whichsig = whichsig,
                              sigcut = sigcut, belowsigcut = belowsigcut, padjcut = padjcut,
                              minsiglength = minsiglength, method = method)
  scorematrix = scorematrix[scorematrix$sameL == 0,,drop=F]
  diffscores = scorematrix[!is.na(scorematrix$score),"score"]
  myconfs = vector()
  for (i in 1:length(myscores)) {
    myconfs[i] = sum(myscores[i] >= diffscores)/length(diffscores)
  }
  names(myconfs) = names(myscores)

  return(myconfs)
}





