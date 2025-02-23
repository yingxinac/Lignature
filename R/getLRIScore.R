################################################################################
# Score LR interactions in specified lr_network in each sample condition.
# Input:
# lr_network: a ligand_receptor network matrix with columns "L", "Lgene", "R", and "Rgene".
# exprinfo: output of function "getExpInfo".
# LRI.method: "cpdb" (mean value of the ligand and the receptor average expression--mean(Lexpr, Rexpr)),
#             "natmi" (product of the ligand and the receptor average expression--Lexpr * Rexpr),
#             or "scsigr" (sqrt(Lexpr*Rexpr)/(sqrt(Lexpr*Rexpr)+a), a is the average expression of all genes in all cells in the input data).
# Output:
# List of LR interaction score and score-specificity data.frames for each cell type/cluster-pair in each sample condition.
################################################################################

#' Score LR interactions in specified lr_network in each sample condition
#'
#' @param lr_network a ligand_receptor network matrix with columns "L", "Lgene", "R", and "Rgene".
#' @param exprinfo output of function "getExpInfo".
#' @param LRI.method "cpdb", "natmi", or "scsigr".
#'
#' @return
#' List of LR interaction score and score-specificity data.frames for each cell type/cluster-pair in each sample condition
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


getLRIScore <- function(lr_network, exprinfo, LRI.method = "scsigr") {
  conds = colnames(exprinfo$avgexprlist[[1]])
  ### get LRI score
  LRIscorelist = list()
  clpairidx = 0
  clpairids = vector()
  for (m in 1:length(exprinfo$avgexprlist)) {
    for(n in 1:length(exprinfo$avgexprlist)) {
      clpairidx = clpairidx + 1
      cl_from = names(exprinfo$avgexprlist)[m]
      cl_to = names(exprinfo$avgexprlist)[n]
      clpairids[clpairidx] = sprintf("%s_to_%s", cl_from, cl_to)
      print(sprintf("%s of %s, %s_to_%s", clpairidx, length(exprinfo$avgexprlist)^2, cl_from, cl_to))
      scorematrix = matrix(NA, nrow = nrow(lr_network), ncol = (ncol(lr_network) + length(conds)))
      colnames(scorematrix) = c(colnames(lr_network), sprintf("SLRI_%s", conds))
      scorematrix = data.frame(scorematrix)
      for (i in 1:nrow(lr_network)) {
        scorematrix[i,colnames(lr_network)] = lr_network[i,]
        Lgenes = str_split(lr_network[i,"Lgene"], "_")[[1]]
        Rgenes = str_split(lr_network[i,"Rgene"], "_")[[1]]
        if (length(intersect(Lgenes, rownames(exprinfo$avgexprlist[[1]]))) == length(Lgenes) &
            length(intersect(Rgenes, rownames(exprinfo$avgexprlist[[1]]))) == length(Rgenes)) {
          for (j in 1:length(conds)) {
            Lexpr = min(exprinfo$avgexprlist[[cl_from]][Lgenes, conds[j]])
            Rexpr = min(exprinfo$avgexprlist[[cl_to]][Rgenes, conds[j]])
            if (LRI.method == "cpdb") {
              scorematrix[i,sprintf("SLRI_%s", conds[j])] = mean(c(Lexpr, Rexpr))
            } else if (LRI.method == "natmi") {
              scorematrix[i,sprintf("SLRI_%s", conds[j])] = Lexpr * Rexpr
            } else if (LRI.method == "scsigr") {
              scorematrix[i,sprintf("SLRI_%s", conds[j])] = sqrt(Lexpr * Rexpr)/(sqrt(Lexpr * Rexpr) + exprinfo$avgexprall)
            }
          }
        }
      }
      LRIscorelist[[clpairidx]] = scorematrix
    }
  }
  names(LRIscorelist) = clpairids
  print("Done scoring.")

  ### get score-specificity
  print("Begin score-specificity.")

  ### Get specificity
  allscorelist = list()
  for (k in 1:length(conds)) {
    whichcol = sprintf("SLRI_%s", conds[k])
    allscoredf = lr_network
    clpairidx = 0
    clpairids = vector()
    for (m in 1:length(exprinfo$avgexprlist)) {
      for(n in 1:length(exprinfo$avgexprlist)) {
        clpairidx = clpairidx + 1
        cl_from = names(exprinfo$avgexprlist)[m]
        cl_to = names(exprinfo$avgexprlist)[n]
        clpairmn = sprintf("%s_to_%s", cl_from, cl_to)
        clpairids[clpairidx] = clpairmn
        allscoredf = cbind.data.frame(allscoredf, LRIscorelist[[clpairmn]][,whichcol,drop=F])
      }
    }
    colnames(allscoredf) = c(colnames(lr_network), clpairids)
    allscorelist[[k]] = allscoredf
  }
  names(allscorelist) = conds

  LRIspeclist = list()
  for (m in 1:length(LRIscorelist)) {
    clpairm = names(LRIscorelist)[m]
    print(sprintf("%s of %s, %s", m, length(LRIscorelist), clpairm))
    scorematrix = matrix(NA, nrow = nrow(lr_network), ncol = (ncol(lr_network) + length(conds)))
    colnames(scorematrix) = c(colnames(lr_network), sprintf("SLRIspec_%s", conds))
    scorematrix = data.frame(scorematrix)
    for (i in 1:nrow(lr_network)) {
      scorematrix[i,colnames(lr_network)] = lr_network[i,]
      for (j in 1:length(conds)) {
        scorevec = allscorelist[[conds[j]]][i,names(LRIscorelist)]
        names(scorevec) = names(LRIscorelist)
        scorevec = scorevec[!is.na(scorevec)]
        if (length(scorevec) > 0 & !is.na(allscorelist[[conds[j]]][i,clpairm])) {
          scorematrix[i,sprintf("SLRIspec_%s", conds[j])] = sum(allscorelist[[conds[j]]][i,clpairm] > scorevec)/length(scorevec)
        }
      }
    }
    LRIspeclist[[m]] = scorematrix
  }
  names(LRIspeclist) = names(LRIscorelist)
  print("Done.")

  LRIscorelist = list(LRIscorelist, LRIspeclist)
  names(LRIscorelist) = c("score", "spec")

  return(LRIscorelist)
}

