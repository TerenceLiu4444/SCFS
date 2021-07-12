library(Seurat)
library(scater)
require(mclust)
require(clue)
library(sparcl)
library(DuoClustering2018)
library(ExperimentHub)

num_cluster_plot = function(data,title){
    x.svd = svd(data)
    k.ch = 1:20
    rat = NULL
    for (k in k.ch){
        print(paste0("title=", title, "k=",k))
        x.cluster = kmeans(x.svd$u[,1:k], centers=k, iter.max = 200, nstart = 100)
        rat = c(rat, 1-sum(x.cluster$withinss)/x.cluster$totss)
    }

    plot(k.ch,rat,xlab = "k",ylab = "Variation explained",main = title)
}

## Convert data into Seurat object
get_seur_new <- function(sce) {
    rownames(sce) <- rowData(sce)$symbol
    seur <- CreateSeuratObject(counts = assays(sce)$counts, min.cells = 3, min.features = 0)
    seur$truth <- colData(sce)$Truth
    Idents(seur) <- seur$truth
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
   return(seur_sub)
}

## Read in data
current_folderpath <- dirname(sys.frame(1)$ofile)
silver2 <- load(file = paste0(current_folderpath, '/datasets/Sce_Dataset2.RData'))
silver2 <- get(silver2)
eh <- ExperimentHub()
zheng <- eh[["EH1532"]]
colData(zheng)$Truth <- gsub('[[:digit:]]+', '', colnames(zheng))
test2 <- load(file = paste0(current_folderpath, '/datasets/Test_2_Kolod.RData'))
test2 <- get(test2)
test4 <- load(file = paste0(current_folderpath, '/datasets/Test_4_Usoskin.RData'))
test4 <- get(test4)

test2$in_X[test2$in_X < 0] <- 0

seur_silver2 <- get_seur_new(silver2)
seur_zheng <- get_seur_new(zheng)
seur_silver2_tcell <- get_tcell(seur_silver2, 'sil')
seur_zheng_tcell <- get_tcell(seur_zheng, 'zheng')

## Get data matrix, k
x_lst <- list(
    test2 = scale(t(test2$in_X)),
    test4 = scale(t(test4$in_X)),
    silver2_tcell = t(seur_silver2_tcell[['RNA']]@scale.data),
    zheng_tcell = t(seur_zheng_tcell[['RNA']]@scale.data)
)
k_vec <- c(test2$n_clust, test4$n_clust, length(unique(seur_silver2_tcell$truth)), length(unique(seur_zheng_tcell$truth)))
names(k_vec) <- names(x_lst)
titles <- paste0(c('Kolod', 'Usoskin', '10x', 'Zheng'), ' k=', k_vec)

setEPS()
postscript('/tmp/choosek.eps')
par(mfrow=c(2,2))
mapply(num_cluster_plot, data = x_lst, title = titles)   
dev.off()
