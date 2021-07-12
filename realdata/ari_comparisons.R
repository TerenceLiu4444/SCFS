library(DuoClustering2018)
library(ExperimentHub)
require(mclust)
library(sparcl)
library(Seurat)
library(scater)
library(FEAST)
library(devtools)
install_github("terenceliu4444/SCFS")
library(SCFS)

current_folderpath <- dirname(sys.frame(1)$ofile)
## Convert data into Seurat object
get_seur_new <- function(sce) {
    rownames(sce) <- rowData(sce)$symbol
    seur <- CreateSeuratObject(counts = assays(sce)$counts, min.cells = 3, min.features = 0)
    seur$truth <- colData(sce)$Truth
    Idents(seur) <- seur$truth
    return(seur)
}
get_seur_old <- function(data) {
    colnames(data$in_X) <- 1:ncol(data$in_X)
    seur <- CreateSeuratObject(counts = exp(data$in_X) - 1, min.cells = 3, min.features = 0)
    seur$truth <- data$true_labs[[1]]
    Idents(seur) <- data$true_labs[[1]]
    seur <- NormalizeData(seur)
    seur <- FindVariableFeatures(object = seur, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seur)
    seur <- ScaleData(seur, features = all.genes)
    seur <- RunPCA(object = seur, features = VariableFeatures(object = seur))
    seur <- JackStraw(seur, num.replicate = 100, verbose = FALSE)
    seur <- FindNeighbors(seur, dims = 1:10)
    seur <- FindClusters(seur)
    return(seur)
}

## Subset T cells
get_tcell <- function(seur, num) {
   if(num == 'sil') {
       seur_sub <- subset(seur, idents = c('CD4+ T Helper2', 'CD4+/CD25 T Reg', 'CD4+/CD45RA+/CD25- Naive T', 'CD4+/CD45RO+ Memory', 'CD8+ Cytotoxic T', 'CD8+/CD45RA+ Naive Cytotoxic'))
   } else {
       seur_sub <- subset(seur, idents = c('cd.t.helper', 'memory.t', 'naive.cytotoxic', 'naive.t', 'regulatory.t'))
   }
   seur_sub <- NormalizeData(seur_sub, normalization.method = "LogNormalize", scale.factor = 10000)
   all.genes <- rownames(seur_sub)
   seur_sub <- ScaleData(seur_sub, features = all.genes)
   seur_sub <- FindVariableFeatures(object = seur_sub, selection.method = "vst", nfeatures = 2000)
   seur_sub <- RunPCA(object = seur_sub, features = VariableFeatures(object = seur_sub))
   seur_sub <- JackStraw(seur_sub, num.replicate = 100, verbose = FALSE)
   seur_sub <- FindNeighbors(seur_sub, dims = 1:10)
   seur_sub <- FindClusters(seur_sub)
   return(seur_sub)
}

## Functions for FEAST
get_cluster <- function(ixs, Y, k) {
	tops = c(500, 1000, 2000)
	cluster_res = NULL
	for (top in tops){
		tmp_ixs = ixs[1:top]
		tmp_markers = rownames(Y)[tmp_ixs]
		tmp_res = TSCAN_Clust(Y, k = k, input_markers = tmp_markers)
		cluster_res[[toString(top)]] = tmp_res
	}
	return(cluster_res)
}

get_mse <- function(Y, cluster_res) {
	## validation step
	Ynorm = Norm_Y(Y)
	mse_res = NULL
	for (top in names(cluster_res)){
		tmp_res = cluster_res[[top]]
		tmp_cluster = tmp_res$cluster
		tmp_mse = cal_MSE(Ynorm = Ynorm, cluster = tmp_cluster)
		mse_res = c(mse_res, tmp_mse)
	}
	names(mse_res) = names(cluster_res)
	return(mse_res)
}


for (num in c('Test_2_Kolod', 'Test_4_Usoskin', 'sil', 'zheng')) {
  result.ari = rep(0, 6)
  
  ## Read in data
  if(num %in% c('Test_2_Kolod', 'Test_4_Usoskin')) {
      data <- load(file = paste0(current_folderpath, '/datasets/', num, '.RData'))
      data <- get(data)
      data$in_X[data$in_X < 0] <- 0
      data_seur <- get_seur_old(data)
  } else if (num == 'sil'){
      data <- load(file = paste0(current_folderpath, '/datasets/Sce_Dataset2.RDat'))
  	  data <- get(data)
      data_seur <- get_seur_new(data)
      data_seur <- get_tcell(data_seur, 'sil')
  } else {
      eh <- ExperimentHub()
      data <- eh[["EH1532"]]
      colData(data)$Truth <- gsub('[[:digit:]]+', '', colnames(data))
      data_seur <- get_seur_new(data)
      data_seur <- get_tcell(data_seur, 'zheng')
  }
  
  ## Get data matrix, cluster, k
  if(num %in% c('Test_2_Kolod', 'Test_4_Usoskin')) {
      x <- scale(t(data$in_X))
      clust_Ind <- as.matrix(data$true_labs)
  } else{
      x <- t(data_seur[['RNA']]@scale.data)
      clust_Ind <- as.matrix(data_seur$truth)
  }
  
  k_vec <- c(3, 4, 5, 5)
  names(k_vec) <- c('Test_2_Kolod', 'Test_4_Usoskin', 'sil', 'zheng')
  k <- k_vec[num]
  
  ## SC-FS1 ##
  scafs1 <- SpectralClusterFeatureSelection(x, k, NULL, FALSE)
  result.ari[1] <- adjustedRandIndex(scafs1$cluster_ids,clust_Ind)
  print('SCAFS1 ARI')
  print(result.ari[1])
  
  ## SC-FS2 ##
  scafs2 <- SpectralClusterFeatureSelection(x, k, NULL, TRUE)
  result.ari[2] <- adjustedRandIndex(scafs2$cluster_ids,clust_Ind)
  print('SCAFS2 ARI')
  print(result.ari[2])
  
  ## spectral ##
  cluster <- SpectralClustering(x, k)
  result.ari[3] =adjustedRandIndex(cluster,clust_Ind)
  print('spectral ARI')
  print(result.ari[3])
  
  ## Seurat
  result.ari[4] = adjustedRandIndex(as.matrix(data_seur$seurat_clusters), clust_Ind)
  print('Seurat ARI')
  print(result.ari[4])
  
  ## FEAST
  y <- process_Y(data_seur[['RNA']]@counts, thre = 3)
  con_res <- Consensus(y, k)
  f_res <- cal_F2(y, con_res$cluster)
  ixs <- order(f_res$F_scores, decreasing = TRUE)
  cluster <- get_cluster(ixs, y, k)
  mse <- get_mse(y, cluster)
  id <- which.min(mse)
  eval_res <- eval_Cluster(cluster[[id]]$cluster, data_seur$truth)
  result.ari[5] = eval_res[[1]]
  print('FEAST ARI')
  print(result.ari[5])
  
  ## sparse kmeans
  km.perm <- KMeansSparseCluster.permute(x,K=k,wbounds=3*10^(0:3),nperms=5)
  spkm = KMeansSparseCluster(x,K=k,wbounds=km.perm$bestw)
  result.ari[6] =adjustedRandIndex(spkm[[1]]$Cs,clust_Ind)
  print('sparse kmeans ARI')
  print(result.ari[6])
  
  names(result.ari) <- c('SCAFS1', 'SCAFS2', 'spectral', 'Seurat', 'FEAST', 'sparse kmeans')
  write.table(result.ari, paste0('/tmp/result_ari_', num, '.txt'))
}