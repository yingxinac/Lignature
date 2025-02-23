################################################################################
# Score each signature in Lignature given signatures of input data.
# Input:
# siglist: siglist of Lignature.
# testsig: input signatures.
# whichsig: "lfc" or "NES".
# sigcut: cutoff value applied to the signatures.
# belowsigcut: "as.is", "remove", or "as.zero", specify how to deal with signature values below sigcut.
# padjcut: cutoff values applied to the adjusted p.values of each entry of the signatures.
# minsiglength: minimum length of signatures under comparison.
# method: similarity measurement metric, one of "cor.pearson", "cor.spearman", "cosine.similarity", "euclidean", "manhattan", "euclidean_mean", "manhattan_mean".
# nump: number of permutations of the input test signature in the calcuation of score-pvalues.
# corpby.permu: if TRUE, score-pvalues represent results from permutations when method is set to "cor.pearson" or "cor.spearman".
# Output:
# A list of two vectors: a vector of signature scores and a vector of score-pvalues.
################################################################################


#' Define similarity metrics
#'
#' @param a vector a
#' @param b vector b
#' @param mymethod similarity metric
#'
#' @return
#' similarity score
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export
mymetric <- function(a, b, mymethod) {
  if (mymethod == "cosine.similarity") {
    d = cosine(a, b)[1,1]
  } else if (mymethod == "euclidean") {
    d = sqrt(sum((a - b)^2))
  } else if (mymethod == "euclidean_mean") {
    d = sqrt(sum((a - b)^2))/length(a)
  } else if (mymethod == "manhattan") {
    d = sum(abs(a - b))
  } else if (mymethod == "manhattan_mean") {
    d = sum(abs(a - b))/length(a)
  } else if (mymethod == "cor.pearson") {
    d = cor.test(a, b, method = "pearson")
  } else if (mymethod == "cor.spearman") {
    d = cor.test(a, b, method = "spearman")
  }
  return(d)
}

#' Score each signature in Lignature given signatures of input data
#'
#' @param siglist siglist of Lignature.
#' @param testsig input signatures.
#' @param whichsig "lfc" or "NES".
#' @param sigcut cutoff value applied to the signatures.
#' @param belowsigcut "as.is", "remove", or "as.zero", specify how to deal with signature values below sigcut.
#' @param padjcut cutoff values applied to the adjusted p.values of each entry of the signatures.
#' @param minsiglength minimum length of signatures to compare.
#' @param method similarity measurement metric, one of "cor.pearson", "cor.spearman", "cosine.similarity", "euclidean", "manhattan", "euclidean_mean", "manhattan_mean".
#' @param nump number of permutations of the input test signature in the calcuation of score-pvalues.
#' @param corpby.permu if TRUE, score-pvalues represent results from permutations when method is set to "cor.pearson" or "cor.spearman".
#'
#' @return
#' a list of two vectors: a vector of signature scores and a vector of score-pvalues
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export

getSigScore <- function(siglist, testsig, whichsig = "lfc",
                        sigcut = 0.1, belowsigcut = "as.is", padjcut = 1,
                        minsiglength = 10, method = "cor.pearson", nump = 1000, corpby.permu = FALSE) {

  sigscores = vector()
  sigscoreps = vector()
  for (i in 1:length(siglist)) {
    if (i %% 50 == 0) {
      print(sprintf("i = %s (in 1:%s)", i, length(siglist)))
    }
    testsigvec = testsig[[whichsig]]
    if (whichsig == "lfc") {
      testpadjvec = testsig[["padj"]]
    } else if (whichsig == "NES") {
      testpadjvec = testsig[["fgseapadj"]]
    }
    mode(testsigvec) = "numeric"
    mode(testpadjvec) = "numeric"
    testsigvec = testsigvec[!is.na(testsigvec)]
    testpadjvec = testpadjvec[!is.na(testpadjvec)]

    refsigvec = siglist[[i]][[whichsig]]
    if (whichsig == "lfc") {
      refpadjvec = siglist[[i]][["padj"]]
    } else if (whichsig == "NES") {
      refpadjvec = siglist[[i]][["fgseapadj"]]
    }
    mode(refsigvec) = "numeric"
    mode(refpadjvec) = "numeric"
    refsigvec = refsigvec[!is.na(refsigvec)]
    refpadjvec = refpadjvec[!is.na(refpadjvec)]

    keepsigentries_test = names(testsigvec)[abs(testsigvec) >= sigcut]
    keepsigentries_ref = names(refsigvec)[abs(refsigvec) >= sigcut]
    rmsigentries_test = names(testsigvec)[abs(testsigvec) < sigcut]
    rmsigentries_ref = names(refsigvec)[abs(refsigvec) < sigcut]
    if (belowsigcut == "remove") {
      testsigvec = testsigvec[keepsigentries_test]
      testpadjvec = testpadjvec[keepsigentries_test]
      refsigvec = refsigvec[keepsigentries_ref]
      refpadjvec = refpadjvec[keepsigentries_ref]
    } else if (belowsigcut == "as.zero") {
      testsigvec[rmsigentries_test] = 0
      refsigvec[rmsigentries_ref] = 0
    }
    testsigvec[names(testpadjvec)[testpadjvec > padjcut]] = 0
    refsigvec[names(refpadjvec)[refpadjvec > padjcut]] = 0

    commonentries = intersect(names(testsigvec), names(refsigvec))
    if (length(commonentries) >= minsiglength) {
      testsigvec = testsigvec[commonentries]
      refsigvec = refsigvec[commonentries]
      if (method %in% c("cor.pearson", "cor.spearman") & corpby.permu == FALSE) {
        scorei = mymetric(a = testsigvec, b = refsigvec, mymethod = method)$estimate
        scorepi = mymetric(a = testsigvec, b = refsigvec, mymethod = method)$p.value
      } else if (method %in% c("cor.pearson", "cor.spearman") & corpby.permu == TRUE) {
        scorei = mymetric(a = testsigvec, b = refsigvec, mymethod = method)$estimate
        scorei_permu = vector()
        for (n in 1:nump) {
          testsigvec_permu = sample(testsigvec)
          names(testsigvec_permu) = names(testsigvec)
          scorei_permu[n] = mymetric(a = testsigvec_permu, b = refsigvec, mymethod = method)$estimate
        }
        scorepi = sum(abs(scorei_permu) > abs(scorei))/nump
      } else if (method == "cosine.similarity") {
        scorei = mymetric(a = testsigvec, b = refsigvec, mymethod = method)
        scorei_permu = vector()
        for (n in 1:nump) {
          testsigvec_permu = sample(testsigvec)
          names(testsigvec_permu) = names(testsigvec)
          scorei_permu[n] = mymetric(a = testsigvec_permu, b = refsigvec, mymethod = method)
        }
        scorepi = sum(abs(scorei_permu) > abs(scorei))/nump
      } else if (method %in% c("euclidean", "manhattan", "euclidean_mean", "manhattan_mean")) {
        scorei = mymetric(a = testsigvec, b = refsigvec, mymethod = method)
        scorei_permu = vector()
        for (n in 1:nump) {
          testsigvec_permu = sample(testsigvec)
          names(testsigvec_permu) = names(testsigvec)
          scorei_permu[n] = mymetric(a = testsigvec_permu, b = refsigvec, mymethod = method)
        }
        scorepi = sum(scorei_permu < scorei)/nump
      }
    } else {
      scorei = NA
      scorepi = NA
    }
    sigscores[i] = scorei
    sigscoreps[i] = scorepi
  }
  names(sigscores) = names(siglist)
  names(sigscoreps) = names(siglist)
  sigscorelist = list(sigscores, sigscoreps)
  names(sigscorelist) = c("score", "pvalue")

  return(sigscorelist)
}



