# Reference: https://ourcodingclub.github.io/tutorials/writing-r-package/
require(mclust)
require(clue)
library(sparcl)

GenerateSyntheticData <-
  function(n, p, s, k, signal_strength, noise_type) {
    #' Generates synthetic data.
    #'
    #' @description This function generates synthetic clustering data.
    #'
    #' It generates synthetic by the following steps:
    #'
    #'   1. generates cluster centers with s informative and p total features according to the following steps:
    #'      1). generates orthonormal matrix of dimension k by s.
    #'      2). multiplies the matrix by signal_strength.
    #'      3). binds (p-s) zero columns to the above matrix.
    #'
    #'   2. randomly generates n by k one-hot cluster assignments.
    #'
    #'   3. generates n by p signal matrix and add scaled standard gaussian or t2 noise.
    #'
    #' @param n int. Number of observations.
    #' @param p int. Number of features.
    #' @param s int. Number of informative features.
    #' @param k int. Number of clusters.
    #' @param signal_strength float. Signal strength.
    #' @param noise_type character. Noise type. Must be either "gaussian" or "t2".
    #' @usage GenerateSyntheticData(n, p, s, k, signal_strength, noise_type)
    #' @return list. The result contains two attributes: $data is the data matrix, $labels contain the cluster ids.
    #'
    #' @references T. Liu, Y. Lu, B. Zhu, H. Zhao (2021). High-dimensional Clustering via Feature Selection with
    #' Applications to Single Cell RNA-seq Data.
    #' @examples
    #' GenerateSyntheticData(n=10, p=10, s=5, k=2, signal_strength=1, noise_type="gaussian")
    #' @export
    # Construct signal matrix.
    rand_mat = matrix(rnorm(k * s), k, s)
    rand_mat.svd = svd(rand_mat)
    cluster_means = t(rand_mat.svd$v)
    cluster_means = cluster_means * signal_strength
    cluster_means = cbind(cluster_means, matrix(0, k, p - s))
    
    # Generate cluster ids.
    labels = sample(1:k, n, replace = T)
    one_hot_cluster_ids = matrix(0, n, k)
    for (i in 1:n) {
      one_hot_cluster_ids[i, labels[i]] = 1
    }
    
    # Generate data matrix.
    signal_mat = (one_hot_cluster_ids %*% cluster_means)
    if (noise_type == "gaussian") {
      data = signal_mat + scale((matrix(rnorm(n * p), n, p)))
    } else if (noise_type == "t2") {
      data = signal_mat + scale((matrix(rt(n * p, 2), n, p)))
    } else {
      stop("noise_type must by either 'gaussian' or 't2'.")
    }
    result = list(data = scale(data), labels = labels)
    return(result)
  }


LloydIteration <-
  function(data,
           num_clusters,
           init_cluster_ids,
           num_iterations) {
    #' Conducts Lloyd iterations
    #'
    #' @description This function conducts Lloyd iterations from an initial clusters assignment.
    #'
    #' @param data matrix. Input data matrix for clustering.
    #' @param num_clusters int. Number of clusters.
    #' @param init_cluster_ids vector. The initial clusters assignment.
    #' @param num_iterations int. Number of Lloyd iterations.
    #' @usage LloydIteration(data, num_clusters, init_cluster_ids, num_iterations)
    #' @return A vector containing the cluster assignment after iterations.
    #' @export
    data = data.matrix(data)
    n = nrow(data)
    d = ncol(data)
    zhat = init_cluster_ids
    for (N in 1:num_iterations) {
      mhat = matrix(0, num_clusters, d)
      for (i in 1:num_clusters) {
        if (length(which(zhat == i)) > 1) {
          mhat[i,] = colMeans(data[which(zhat == i),])
        }
      }
      for (i in 1:n) {
        zhat[i] = which.min(rowSums((rbind(data[i,])[rep(1, num_clusters), ] - mhat) ^ 2))
      }
    }
    return(zhat)
  }

SpectralClustering <- function(data, num_clusters) {
  #' Conducts spectral clustering algorithm.
  #'
  #' @description This function implements spectral clustering algorithm.
  #'
  #' It conducts k-means on top k left singular vectors of the data matrix.
  #'
  #' @param data matrix. Input data matrix for clustering.
  #' @param num_clusters int. Number of clusters.
  #' @usage SpectralClustering(data, num_clusters)
  #' @return A vector containing the cluster assignment after iterations.
  #' @references T. Liu, Y. Lu, B. Zhu, H. Zhao (2021). High-dimensional Clustering via Feature Selection with
  #' Applications to Single Cell RNA-seq Data.
  #' @examples
  #' synthetic_data <- GenerateSyntheticData(n=10, p=10, s=5, k=2, signal_strength=1, noise_type="gaussian")
  #' label.est <- SpectralClustering(synthetic_data$data, 2)
  #' @export
  data.svd = svd(data)
  kmeans.result = kmeans(
    data.svd$u[, 1:num_clusters],
    centers = num_clusters,
    iter.max = 200,
    nstart = 100
  )
  return(kmeans.result$cluster)
}

SumOfSquares <- function(x) {
  return(sum((x - mean(x)) ^ 2))
}

SpectralClusterFeatureSelection <-
  function(data,
           num_clusters,
           init_cluster_ids,
           use_lloyd_iteration) {
    #' Conducts spectral clustering followed by feature selection and another clustering.
    #'
    #' @description This function implements the main algorithm of SCFS.
    #'
    #' It conducts clustering by the following steps:
    #'
    #'   1. conducts spectral clustering to obtain initial guess of cluster assignments.
    #'
    #'   2. uses the initial guess to compute R^2 and select informative features.
    #'
    #'   3. conducts spectral clustering again (SCFS1) followed by Lloyd iteration (optionally for SCFS2).
    #'
    #' @param data matrix. Input data matrix for clustering.
    #' @param num_clusters int. Number of clusters.
    #' @param init_cluster_ids vector. The initial clusters assignment. If NULL, the initial guess will
    #' be obtained by spectral clustering. Otherwise, step 1 is skipped and the feature selection will be
    #' applied based on the initial guess.
    #' @param use_lloyd_iteration bool. If TRUE, conduct Lloyd iteration at step 3. Otherwise, no Lloyd iteration.
    #' @usage SpectralClusterFeatureSelection(data, num_clusters, init_cluster_ids, use_lloyd_iteration)
    #' @return A list containing two attributes: $cluster_ids and $info_feat_ids, where $cluster_ids
    #' contains the estimated cluster assignments and $info_feat_ids contains the selected feature indices.
    #' @references T. Liu, Y. Lu, B. Zhu, H. Zhao (2021). High-dimensional Clustering via Feature Selection with
    #' Applications to Single Cell RNA-seq Data.
    #' @examples
    #' synthetic_data <- GenerateSyntheticData(n=10, p=10, s=5, k=2, signal_strength=1, noise_type="gaussian")
    #' scfs <- SpectralClusterFeatureSelection(data=synthetic_data$data, num_clusters=2, init_cluster_ids=NULL, use_lloyd_iteration=FALSE)
    #' @export
    n = nrow(data)
    p = ncol(data)
    k = num_clusters
    if (is.null(init_cluster_ids)) {
      # Step 1: use spectral clustering for initial guess.
      init_cluster_ids = SpectralClustering(data, num_clusters)
    }
    
    # Step 2: feature selection using R^2.
    cond_ssq = matrix(0, k, p)
    total_var = apply(data, 2, var)
    
    for (m in 1:k) {
      if (length(which(init_cluster_ids == m)) > 1) {
        cond_ssq[m,] = apply(data[which(init_cluster_ids == m),], 2, SumOfSquares)
      } else {
        cond_ssq[m,] = 0.
      }
    }
    info.score = apply(cond_ssq, 2, sum) /  total_var
    info.ind = which(info.score < 0.9 * n)
    
    # In case selected features are too few.
    if (length(info.ind) < k) {
      info.ind = order(info.score)[1:min(10 * k, p)]
    }
    data.filtered = data[, info.ind]
    
    # Step 3: additional spectral clustering with (optional) Lloyd iterations.
    scfs.cluster_ids = SpectralClustering(data.filtered, num_clusters)
    if (use_lloyd_iteration) {
      scfs.cluster_ids = LloydIteration(data.filtered, num_clusters, scfs.cluster_ids, 3 * log(n))
    }
    data.filtered.svd = svd(data.filtered)
    scfs.result = list(cluster_ids = scfs.cluster_ids, info_feat_ids = info.ind)
    return(scfs.result)
  }

ComputeErrorRate <- function(zhat, z) {
  #' Computes error rate given estimated and real cluster ids.
  #'
  #' @param zhat vector. Estimated cluster ids.
  #' @param z vector. True cluster ids.
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