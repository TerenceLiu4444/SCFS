library(sparcl)
library(ggplot2)
library(parallel)
library(devtools)
install_github("terenceliu4444/SCFS")
library(SCFS)
library(dplyr)

N_REP = 50

k = 4
s = 100
p = 500
lp = log(p)
signal.ch = c(5., 10.)
n.ch = ceiling(c(10, 50, 100) * lp)
err.ch = c(0.05, 0.1, 0.15, 0.2, 0.3)
res = array(0, c(length(n.ch), length(err.ch), length(signal.ch), N_REP))

df_all = NULL
for (ind.sig in 1:length(signal.ch)) {
  for (ind.n in 1:length(n.ch)) {
    for (ind.err in 1:length(err.ch)) {
      f1 = rep(0, N_REP)
      for (rr in 1:N_REP) {
        signal = signal.ch[ind.sig]
        n = n.ch[ind.n]
        err = err.ch[ind.err]
        synthetic_data = GenerateSyntheticData(n, p, s, k, signal, noise_type = "gaussian")
        data = scale(synthetic_data$data)
        labels = synthetic_data$labels
        PB = matrix(err / (k - 1), k, k)
        for (i in 1:k) {
          PB[i, i] = 1 - err
        }
        init = NULL
        for (i in 1:n) {
          init = c(init, which(as.vector(rmultinom(1, 1, PB[labels[i], ]) == 1)))
        }
        
        scfs = SpectralClusterFeatureSelection(data, k, init, FALSE)
        est_info_feats = scfs$info_feat_ids
        true_info_feats = seq(1, s)
        n_true_pos = length(intersect(est_info_feats, true_info_feats))
        prec = n_true_pos / length(est_info_feats)
        recall = n_true_pos / length(true_info_feats)
        f1[rr] = (2 * prec * recall) / (prec + recall)
      }
      df = data.frame(n=n, err=err, signal=signal, f1=f1)
      df_all = rbind(df_all, df)
    }
  }
}

result = df_all %>% group_by(n, err, signal) %>%  summarise(across(c("f1"), list(
  mean = mean,
  sd = sd
)))
print(Unstack(as.data.frame(result), 4, 2, numeric(0),c(1,3)))
print(Unstack(as.data.frame(result), 5, 2, numeric(0),c(1,3)))
