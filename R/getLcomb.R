################################################################################
# Perform combinatorial analysis and score the ligands.
# Input:
# siglist: siglist of Lignature.
# testsig: input signatures.
# siggroups: list of signatureid-vectors for each ligand in Lignature.
# sigLs: candidate ligands included in Lignature.
# whichsig: "lfc" or "NES".
# combmethod: "lasso", "ridge", or "elasticnet".
# Output:
# A list with two entries: a vector of ligand scores and a list of signature scores corresponding to each ligand.
################################################################################

#' collapse a vector
#'
#' @param x a vector
#'
#' @return
#' x collapsed with "_"
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table circlize tidyverse randomcoloR hrbrthemes ggrepel caret glmnet
#' @export

pastecollapse <- function(x) {
  return(paste(x, collapse = "_"))
}


#' Perform combinatorial analysis and score the ligands
#'
#' @param siglist siglist of Lignature.
#' @param testsig input signatures.
#' @param siggroups list of signatureid-vectors for each ligand in Lignature.
#' @param sigLs candidate ligands included in Lignature.
#' @param whichsig "lfc" or "NES".
#' @param combmethod "lasso", "ridge", or "elasticnet".
#'
#' @return
#' A vector of ligand scores, and a list of signature scores corresponding to each ligand
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table circlize tidyverse randomcoloR hrbrthemes ggrepel caret glmnet
#' @export


getLcomb <- function(siglist, testsig, siggroups, sigLs,
                     whichsig = "lfc", combmethod = "lasso") {
  for (i in 1:length(siglist)) {
    siglist[[i]]$sigLsep = str_split(siglist[[i]]$sigL, "and")[[1]]
  }
  allentries = names(testsig[[whichsig]])
  for (i in 1:length(siglist)) {
    allentries = union(allentries, names(siglist[[i]][[whichsig]]))
  }
  usematrix = matrix(0, nrow = length(allentries), ncol = (length(siglist)+1))
  rownames(usematrix) = allentries
  colnames(usematrix) = c("test", names(siglist))
  usematrix[names(testsig[[whichsig]]),"test"] = testsig[[whichsig]]
  for (i in 1:length(siglist)) {
    usematrix[names(siglist[[i]][[whichsig]]),names(siglist)[i]] = siglist[[i]][[whichsig]]
  }
  usematrix[is.na(usematrix)] = 0
  suppressWarnings({
    usematrix = data.frame(usematrix)
    y = usematrix %>% dplyr::select(test) %>% scale(center = TRUE, scale = FALSE) %>% as.matrix()
    X = usematrix %>% dplyr::select(-test) %>% as.matrix()
    lambdas_to_try = 10^seq(-3, 5, length.out = 100)
    if (combmethod == "ridge") {
      mycv = cv.glmnet(X, y, alpha = 0, lambda = lambdas_to_try, standardize = TRUE, nfolds = 5)
      lambda_cv = mycv$lambda.min
      model_cv = glmnet(X, y, alpha = 0, lambda = lambda_cv, standardize = TRUE)
      sigbeta = model_cv$beta
    } else if (combmethod == "lasso") {
      mycv = cv.glmnet(X, y, alpha = 1, lambda = lambdas_to_try, standardize = TRUE, nfolds = 5)
      lambda_cv = mycv$lambda.min
      model_cv = glmnet(X, y, alpha = 1, lambda = lambda_cv, standardize = TRUE)
      sigbeta = model_cv$beta
    } else if (combmethod == "elasticnet") {
      train_control = trainControl(method = "repeatedcv",
                                   number = 5,
                                   repeats = 5,
                                   search = "random",
                                   verboseIter = TRUE)
      model_cv = train(test ~ .,
                       data = cbind(y, X),
                       method = "glmnet",
                       preProcess = c("center", "scale"),
                       tuneLength = 25,
                       trControl = train_control)
      sigbeta = coef(model_cv$finalModel, s = model_cv$bestTune$lambda)
    }
    y_hat_cv = predict(model_cv, X)
    ssr_cv = t(y - y_hat_cv) %*% (y - y_hat_cv)
    rsq_cv = cor(y, y_hat_cv)^2
    print(sprintf("ssr_cv = %s, rsq_cv = %s", ssr_cv, rsq_cv))
    sigbeta = sigbeta[,1]
  })
  Lsigbeta = list()
  Lbeta = vector()
  for (i in 1:length(sigLs)) {
    Libeta = sigbeta[intersect(names(sigbeta), siggroups[[sigLs[i]]])]
    Libeta = Libeta[!is.na(Libeta)]
    Lsigbeta[[i]] = Libeta
    Lbeta[i] = max(Libeta)
  }
  names(Lsigbeta) = sigLs
  names(Lbeta) = sigLs

  mylist = list(Lbeta, Lsigbeta)
  names(mylist) = c("Lscores", "sigscores")

  return(mylist)
}






