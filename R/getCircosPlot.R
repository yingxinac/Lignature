################################################################################
# Create circos plot of identified LR interactions.
# Inputs:
# plotLRinfolist: the "LRinfolist_fil" entry of the output from the function "filterLRinfo".
# plotcond: specifies from which condition the interaction scores will be plot,"condmax" (the maximum interaction score across conditions) or one of the conditions under comparison.
# cls_from: sender cell-types/clusters.
# cl_to: specifies receiver cell-type/cluster.
# mycols: a vector of colors for related cell-types/clusters.
# width_same_ligand_cluster: gap width between ligand nodes in the same cluster.
# width_different_ligand_cluster: gap width between different ligand clusters.
# width_ligand_receptor: gap width between the ligand nodes and the receptor nodes.
# width_same_receptor_cluster: gap width between receptor nodes in the same cluster.
# width_different_receptor_cluster: gap width between different receptor clusters.
# cplotthresh: specifies threshold of the SLRI values in the plot. Only plot LR interactions with SLRI > cplotthresh.
# cex: font size.
# transupperbase, transidx: parameters of edge transparency adjustment.
# Output:
# Circos plot of identified LR interactions.
################################################################################

#' Create circos plot of identified LR interactions
#'
#' @param plotLRinfolist the "LRinfolist_fil" entry of the output from the function "filterLRinfo".
#' @param plotcond specifies from which condition the interaction scores will be plot,"condmax" (the maximum interaction score across conditions) or one of the conditions under comparison.
#' @param cls_from sender cell-types/clusters.
#' @param cl_to specifies receiver cell-type/cluster.
#' @param mycols a vector of colors for related cell-types/clusters.
#' @param width_same_ligand_cluster gap width between ligand nodes in the same cluster.
#' @param width_different_ligand_cluster gap width between different ligand clusters.
#' @param width_ligand_receptor gap width between the ligand nodes and the receptor nodes.
#' @param width_same_receptor_cluster gap width between receptor nodes in the same cluster.
#' @param width_different_receptor_cluster gap width between different receptor clusters.
#' @param cplotthresh specifies threshold of the SLRI values in the plot. Only plot LR interactions with SLRI > cplotthresh.
#' @param cex font size.
#' @param transupperbase parameter of edge transparency adjustment.
#' @param transidx parameter of edge transparency adjustment.
#'
#' @return
#' Circos plot of identified LR interactions
#' @import dplyr Seurat stringr RColorBrewer pheatmap ggplot2 fgsea graphite lsa devtools data.table geneconverter circlize tidyverse randomcoloR hrbrthemes ggrepel
#' @export


getCircosPlot <- function(plotLRinfolist, plotcond, cls_from, cl_to, mycols,
                        width_same_ligand_cluster, width_different_ligand_cluster, width_ligand_receptor,
                        width_same_receptor_cluster, width_different_receptor_cluster,
                        cplotthresh, cex, transupperbase, transidx) {

  mymatrix = vector()
  for (i in 1:length(cls_from)) {
    clpair = sprintf("%s_to_%s", cls_from[i], cl_to)
    matrixi = plotLRinfolist[[clpair]][,c("Lgene", "R", sprintf("SLRI_%s", plotcond)), drop=F]
    matrixi$clfrom = cls_from[i]
    matrixi$clto = cl_to
    matrixi$clfromL = sprintf("%s_%s", cls_from[i], matrixi$Lgene)
    matrixi$cltoR = sprintf("%s_%s", cl_to, matrixi$R)
    matrixi$LR = sprintf("%s_to_%s", matrixi$clfromL, matrixi$cltoR)
    rownames(matrixi) = matrixi$LR
    mymatrix = rbind.data.frame(mymatrix, matrixi)
  }
  myLclusterdf = unique(mymatrix[,c("clfromL", "clfrom")])
  myLcluster = myLclusterdf$clfrom
  names(myLcluster) = myLclusterdf$clfromL
  myRclusterdf = unique(mymatrix[,c("cltoR", "clto")])
  myRcluster = myRclusterdf$clto
  names(myRcluster) = myRclusterdf$cltoR
  myLRclusterdf = cbind.data.frame(rownames(mymatrix), mymatrix$clfrom)
  myLRcluster = myLRclusterdf$`mymatrix$clfrom`
  names(myLRcluster) = myLRclusterdf$`rownames(mymatrix)`
  clsunion = union(cls_from, cl_to)
  cltonum = cbind(clsunion, as.character(c(0:(length(clsunion)-1))))
  rownames(cltonum) = cltonum[,1]
  colnames(cltonum) = c("cl", "num")
  for (i in 1:nrow(cltonum)) {
    myLcluster[myLcluster == cltonum[i,"cl"]] = cltonum[i,"num"]
    myRcluster[myRcluster == cltonum[i,"cl"]] = cltonum[i,"num"]
    myLRcluster[myLRcluster == cltonum[i,"cl"]] = cltonum[i,"num"]
  }
  myplotmatrix = mymatrix[,c("clfromL", "cltoR", sprintf("SLRI_%s", plotcond), "LR"), drop=F]
  colnames(myplotmatrix)[colnames(myplotmatrix)==sprintf("SLRI_%s", plotcond)] = "SLRI"
  names(mycols) = as.character(c(0:(length(clsunion)-1)))

  ######################## Circos data
  LR_df = unique(myplotmatrix)
  colnames(LR_df)[1:2] = c('L_name', 'R_name')
  LR_tb = tibble(LR_df)

  ### The user pre-defined gene and LR pair clusters
  LRP_clusters = cbind(myplotmatrix[,'clfromL'], myplotmatrix[,'cltoR'], myLcluster[myplotmatrix[,'clfromL']], myRcluster[myplotmatrix[,'cltoR']], myLRcluster[myplotmatrix[,'LR']], myplotmatrix[,'LR'])
  rownames(LRP_clusters) = NULL
  colnames(LRP_clusters) = c('L', 'R', 'L_cluster', 'R_cluster', 'LRcluster', 'LR')

  # L cluster
  L_clusterset = unique(LRP_clusters[,"L_cluster"])
  L_cluster_list = list()
  for (i in 1:length(L_clusterset)) {
    L_cluster_list[[i]] = unique(LRP_clusters[LRP_clusters[,"L_cluster"]==L_clusterset[i],"L"])
  }
  names(L_cluster_list) = L_clusterset
  L_cluster_vec = vector()
  for (i in 1:length(L_clusterset)) {
    L_cluster_vec = c(L_cluster_vec, rep(L_clusterset[i], length(L_cluster_list[[i]])))
  }
  L_name_vec = vector()
  for (i in 1:length(L_clusterset)) {
    L_name_vec = c(L_name_vec, L_cluster_list[[i]])
  }
  Lcluster_indication_tb = tibble(L_cluster = L_cluster_vec, L_name = L_name_vec)

  # R cluster
  R_clusterset = unique(LRP_clusters[,"R_cluster"])
  R_cluster_list = list()
  for (i in 1:length(R_clusterset)) {
    R_cluster_list[[i]] = unique(LRP_clusters[LRP_clusters[,"R_cluster"]==R_clusterset[i],"R"])
  }
  names(R_cluster_list) = R_clusterset
  R_cluster_vec = vector()
  for (i in 1:length(R_clusterset)) {
    R_cluster_vec = c(R_cluster_vec, rep(R_clusterset[i], length(R_cluster_list[[i]])))
  }
  R_name_vec = vector()
  for (i in 1:length(R_clusterset)) {
    R_name_vec = c(R_name_vec, R_cluster_list[[i]])
  }
  Rcluster_indication_tb = tibble(R_cluster = R_cluster_vec, R_name = R_name_vec)

  # LR cluster
  LR_clusterset = unique(LRP_clusters[,"LRcluster"])
  LR_cluster_list = list()
  for (i in 1:length(LR_clusterset)) {
    LR_cluster_list[[i]] = unique(LRP_clusters[LRP_clusters[,"LRcluster"]==LR_clusterset[i], 'LR'])
  }
  names(LR_cluster_list) = LR_clusterset
  LR_cluster_vec = vector()
  for (i in 1:length(LR_clusterset)) {
    LR_cluster_vec = c(LR_cluster_vec, rep(LR_clusterset[i], length(LR_cluster_list[[i]])))
  }
  LR_name_vec = vector()
  for (i in 1:length(LR_clusterset)) {
    LR_name_vec = c(LR_name_vec, LR_cluster_list[[i]])
  }
  LRcluster_indication_tb = tibble(LR_cluster = LR_cluster_vec, LR = LR_name_vec)

  LR_addcluster = unique(inner_join(LR_tb, Lcluster_indication_tb))
  LR_addcluster = unique(inner_join(LR_addcluster, Rcluster_indication_tb))
  LR_addcluster = unique(inner_join(LR_addcluster, LRcluster_indication_tb))

  ### Colors
  grid_col_tbl_L = tibble(L_cluster = mycols %>% names(), color_L_cluster = mycols)
  grid_col_tbl_R = tibble(R_cluster = mycols %>% names(), color_R_cluster = mycols)
  grid_col_tbl_LR = tibble(LR_cluster = mycols %>% names(), color_LR_cluster = mycols)

  LR_addcluster = LR_addcluster %>% mutate(L_name = paste(L_name," "))

  LR_addcluster = LR_addcluster %>%
    inner_join(grid_col_tbl_L) %>%
    inner_join(grid_col_tbl_R) %>%
    inner_join(grid_col_tbl_LR)
  LR_circle = LR_addcluster %>% dplyr::select(L_name, R_name, SLRI)

  L_color = LR_addcluster %>% distinct(L_name,color_L_cluster)
  grid_L_color = L_color$color_L_cluster %>% set_names(L_color$L_name)

  R_color = LR_addcluster %>% distinct(R_name,color_R_cluster)
  grid_R_color = R_color$color_R_cluster %>% set_names(R_color$R_name)

  LR_color = LR_addcluster %>% distinct(LR,color_LR_cluster)
  grid_LR_color = LR_color$color_LR_cluster %>% set_names(LR_color$LR)

  grid_col =c(grid_L_color,grid_R_color)

  ### Edge transparency
  transparency = LR_addcluster %>% mutate(transparency = (transupperbase- SLRI)*transidx) %>% .$transparency

  ### Order the ligands and the receptors
  L_order = L_name_vec %>% c(paste(.," ")) %>% intersect(LR_addcluster$L_name)
  R_order = R_name_vec %>% c(paste(.,"")) %>% intersect(LR_addcluster$R_name)
  order = c(L_order,R_order)

  ### Gaps
  gaps_L = vector()
  for (i in 1:length(L_clusterset)) {
    gaps_L = c(gaps_L,
               rep(width_same_ligand_cluster,
                   times = (LR_addcluster %>% filter(L_cluster == L_clusterset[i]) %>% distinct(L_name) %>% nrow() -1)),
               width_different_ligand_cluster)
  }
  gaps_L = gaps_L[1:(length(gaps_L)-1)]

  gaps_R = vector()
  for (i in 1:length(R_clusterset)) {
    gaps_R = c(gaps_R,
               rep(width_same_receptor_cluster,
                   times = (LR_addcluster %>% filter(R_cluster == R_clusterset[i]) %>% distinct(R_name) %>% nrow() -1)),
               width_different_receptor_cluster)
  }
  gaps_R = gaps_R[1:(length(gaps_R)-1)]

  gaps = c(gaps_L, width_ligand_receptor, gaps_R, width_ligand_receptor)

  ### Plot
  LR_circle = LR_addcluster %>% dplyr::select(L_name,R_name, SLRI, LR_cluster, color_LR_cluster)
  circos.clear()
  circos.par(gap.degree = gaps)
  chordDiagram(LR_circle, directional = 1, order=order, link.sort = TRUE, link.decreasing = FALSE, col = LR_circle$color_LR_cluster,
               grid.col = grid_col, transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow", link.visible = LR_circle$SLRI > cplotthresh, annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.2))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = cex)
  }, bg.border = NA)
}

