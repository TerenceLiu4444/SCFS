library(DuoClustering2018)
library(ExperimentHub)
require(mclust)
library(sparcl)
library(Seurat)
library(scater)
library(devtools)
install_github("terenceliu4444/SCFS")
library(SCFS)
library(cowplot)
library(ggplot2)

## Convert data into Seurat object
get_seur_old <- function(data) {
    
    colnames(data$in_X) <- 1:ncol(data$in_X)
    seur <- CreateSeuratObject(counts = data$in_X, min.cells = 3, min.features = 0)
    seur$truth <- data$true_labs[[1]]
    Idents(seur) <- data$true_labs[[1]]
    seur <- NormalizeData(seur)
    all.genes <- rownames(seur)
    seur <- ScaleData(seur, features = all.genes)
    seur_data <- GetAssayData(seur, slot = "data")
    seur_dist <- dist(t(as.matrix(seur_data)))
    seur[["umap"]] <- RunUMAP(seur_dist)

    return(seur)
}

get_seur_new <- function(sce) {
    rownames(sce) <- rowData(sce)$symbol
    seur <- CreateSeuratObject(counts = assays(sce)$counts, min.cells = 3, min.features = 0)
    seur$truth <- colData(sce)$Truth
    Idents(seur) <- seur$truth
    return(seur)
}

seur_rerun <- function(seur) {
    seur <- NormalizeData(seur)
    seur_data <- GetAssayData(seur, slot = "data")
    seur_dist <- dist(t(as.matrix(seur_data)))
    seur[["umap"]] <- RunUMAP(seur_dist)
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
   seur_data <- GetAssayData(seur_sub, slot = "data")
   seur_dist <- dist(t(as.matrix(seur_data)))
   seur_sub[["umap"]] <- RunUMAP(seur_dist)
   return(seur_sub)
}

## Read in data
current_folderpath <- dirname(sys.frame(1)$ofile)
test2 <- load(file = paste0(current_folderpath, '/datasets/Test_2_Kolod.RData'))
test2 <- get(test2)
test2$in_X[test2$in_X < 0] <- 0

test4 <- load(file = paste0(current_folderpath, '/datasets/Test_4_Usoskin.RData'))
test4 <- get(test4)
test4$in_X[test4$in_X < 0] <- 0

sil <- load(file = paste0(current_folderpath, '/datasets/Sce_Dataset2.RData'))
sil <- get(sil)
sil_seur <- get_seur_new(sil)
sil_seur <- get_tcell(sil_seur, 'sil')

eh <- ExperimentHub()
zheng <- eh[["EH1532"]]
colData(zheng)$Truth <- gsub('[[:digit:]]+', '', colnames(zheng))
zheng_seur <- get_seur_new(zheng)
zheng_seur <- get_tcell(zheng_seur, 'zheng')

## Get data matrix, k
x_lst <- list(
    test2 = scale(t(test2$in_X)),
    test4 = scale(t(test4$in_X)),
    silver2_tcell = t(sil_seur[['RNA']]@scale.data),
    zheng_tcell = t(zheng_seur[['RNA']]@scale.data)
)
clust_Ind_lst <- list(
    test2 = as.matrix(test2$true_labs),
    test4 = as.matrix(test4$true_labs),
    silver2_tcell = as.matrix(sil_seur$truth),
    zheng_tcell = as.matrix(zheng_seur$truth)
)
k_vec <- c(3, 4, 5, 5)
names(k_vec) <- c('Test_2_Kolod', 'Test_4_Usoskin', 'sil', 'zheng')

test2$in_X <- exp(test2$in_X) - 1
test4$in_X <- exp(test4$in_X) - 1
test2_seur <- get_seur_old(test2)
test4_seur <- get_seur_old(test4)

seur_lst <- list(
    test2 = test2_seur,
    test4 = test4_seur,
    sil = sil_seur,
    zheng = zheng_seur
)

## Plot UMAP with all features
titles <- c('Kolod', 'Usoskin', '10x', 'Zheng')
plotlist <- list()
for(i in 1:length(seur_lst)) {
    plotlist[[i]] <- DimPlot(seur_lst[[i]], reduction = "umap") +
                        ggtitle(paste(titles[i], 'all features'))
}
       
## SC-FS ##
res_lst <- mapply(SpectralClusterFeatureSelection, data = x_lst, num_clusters = k_vec, MoreArgs = list(init_cluster_ids = NULL, use_lloyd_iteration = TRUE), SIMPLIFY = FALSE)

## Get selected features
festures_lst <- lapply(res_lst, function(res) res$info_feat_ids)
                        
## Subset Seurat
seur_lst_sub <- mapply(function(seur, features) seur[features,], seur = seur_lst, features = festures_lst)                    
seur_lst_sub <- lapply(seur_lst_sub, seur_rerun)

## Plot UMAP with selected features
for(i in 1:length(seur_lst)) {
    plotlist[[i + 4]] <- DimPlot(seur_lst_sub[[i]], reduction = "umap") +
                        ggtitle(paste(titles[i], 'selcted features'))
} 
pll <- plot_grid(plotlist = plotlist, ncol=4, rel_widths = rep(c(1, 1, 1.5, 1.3), 2))
ggsave('/tmp/compare_umap.eps', pll, width = 25, height = 10)
                       