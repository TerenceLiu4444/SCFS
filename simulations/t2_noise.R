library(sparcl)
library(ggplot2)
library(parallel)
library(clue)
library(devtools)
install_github("terenceliu4444/SCFS")
library(SCFS)
library(dplyr)
library(foreach)
library(doParallel)
library(doSNOW)
library(tcltk)

N_REP = 50
k = 4
s = 500
p = 8000
lp = log(p)
n.ch = ceiling(c(15, 20, 25, 30) * lp)
signal_strength = 6

ComputeErrorRate <- function(zhat, z) {
    n = length(zhat)
    tab = table(zhat, z)
    k = dim(tab)[1]
    match = solve_LSAP(tab, maximum = T)[1:k]
    a = rep(0, k)
    for (i in 1:k) {
        a[i] = tab[i, match[i]]
    }
    return (1 - sum(a) / n)
}

ConductSimulation <- function(n, p, s, k, signal_strength, rr) {
    print(sprintf("n=%d, r=%d", n, rr))
    synthetic_data = GenerateSyntheticData(n, p, s, k,
                                           signal_strength,
                                           noise_type = "t2")
    data = synthetic_data$data
    labels = synthetic_data$labels
    labels.spec = SpectralClustering(data, k)
    sp.res = ComputeErrorRate(labels, labels.spec)
    labels.spec_lloyd = LloydIteration(data, k, labels.spec, 3 * log(n) + 1)
    sp_lloyd.res = ComputeErrorRate(labels, labels.spec_lloyd)
    km.perm <-
        KMeansSparseCluster.permute(data,
                                    K = k,
                                    wbounds = 3 * 10 ^ (0:3),
                                    nperms = 5)
    wit = KMeansSparseCluster(data, K = k, wbounds = km.perm$bestw)
    wit.res = ComputeErrorRate(labels, wit[[1]]$Cs)
    scafs = SpectralClusterFeatureSelection(data, k, NULL, use_lloyd_iteration =
                                                FALSE)
    scafs.res = ComputeErrorRate(labels, scafs$cluster_ids)
    scafs_lloyd = SpectralClusterFeatureSelection(data, k, NULL, use_lloyd_iteration =
                                                      TRUE)
    scafs_lloyd.res = ComputeErrorRate(labels, scafs_lloyd$cluster_ids)
    df = data.frame(
        n = n,
        spectral = sp.res,
        spec_lloyd = sp_lloyd.res,
        sp_kmeans = wit.res,
        scfs1 = scafs.res,
        scfs2 = scafs_lloyd.res
    )
    return(df)
}

df_all = NULL

for (n.ind in 1:length(n.ch)) {
    n = n.ch[n.ind]
    cores = detectCores()
    cl <- makeSOCKcluster(cores[1] - 1)
    registerDoSNOW(cl)
    pb <- tkProgressBar(max=N_REP)
    progress <- function(n) setTkProgressBar(pb, n)
    opts <- list(progress=progress)
    df_n <-
        foreach(
            rr = 1:N_REP,
            .packages = c("sparcl", "SCFS", "clue"),
            .combine = rbind,
            .options.snow=opts
        ) %dopar% {
            df = ConductSimulation(n, p, s, k, signal_strength, rr)
            df
        }
    stopCluster(cl)
    df_all = rbind(df_all, df_n)
}

print(df_all %>% group_by(n) %>%  summarise(across(everything(), list(
    mean = mean
))), n=100)

print(df_all %>% group_by(n) %>%  summarise(across(everything(), list(
    sd = sd
))), n=100)
