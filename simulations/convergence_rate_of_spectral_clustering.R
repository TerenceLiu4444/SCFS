library(sparcl)
library(ggplot2)
library(parallel)
library(devtools)
install_github("terenceliu4444/SCFS")
library(SCFS)
N_REP = 50

####################
# Checks rate of p #
####################
CheckErrorRateP <- function() {
    k = 4
    pch = seq(100, 1000, 10)
    n = 100
    signal_strength = 4
    res_sp = NULL
    for (p in pch) {
        res_tmp = NULL
        for (rr in 1:N_REP) {
            # Generates synthetic data.
            synthetic_data = GenerateSyntheticData(n, p, p, k, signal_strength, noise_type="gaussian")
            data = synthetic_data$data
            true_labels = synthetic_data$labels
    
            # Conducts spectral clustering.
            est_labels = SpectralClustering(data, k)
            spec_err = ComputeErrorRate(est_labels, true_labels)
            res_tmp = c(res_tmp, spec_err)
        }
        res_sp = c(res_sp, mean(res_tmp))
    }
    
    data = data.frame(pch = pch, res_sp = res_sp)
    ggplot(data, aes(x = pch, y = res_sp)) +
        geom_point(shape = 1) +
        xlab("p") + ylab("error rate") +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14, face = "bold"))
    ggsave("figures/err_rate_p.eps")
}

CheckErrorRateP()

####################
# Checks rate of n #
####################
CheckErrorRateN <- function() {
    k = 4
    nch = seq(100, 1000, 10)
    p = 100
    signal_strength = 4
    res_sp = NULL
    for (n in nch) {
        res_tmp = NULL
        for (rr in 1:N_REP) {
            # Generates synthetic data.
            synthetic_data = GenerateSyntheticData(n, p, p, k, signal_strength, noise_type="gaussian")
            data = synthetic_data$data
            true_labels = synthetic_data$labels
            
            # Conducts spectral clustering.
            est_labels = SpectralClustering(data, k)
            spec_err = ComputeErrorRate(est_labels, true_labels)
            res_tmp = c(res_tmp, spec_err)
        }
        res_sp = c(res_sp, mean(res_tmp))
    }
    
    data = data.frame(nch = nch, res_sp = res_sp)
    ggplot(data, aes(x = nch, y = res_sp)) +
        geom_point(shape = 1) +
        xlab("n") + ylab("error rate") +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14, face = "bold"))
    ggsave("figures/err_rate_n.eps")
}

CheckErrorRateN()

#########################
# Checks rate of signal #
#########################
CheckErrorRateSignal <- function() {
    k = 4
    n = 100
    p = 100
    sig_ch = seq(2, 5, 0.05)
    res_sp = NULL
    for (signal_strength in sig_ch) {
        res_tmp = NULL
        for (rr in 1:N_REP) {
            # Generates synthetic data.
            synthetic_data = GenerateSyntheticData(n, p, p, k, signal_strength, noise_type="gaussian")
            data = synthetic_data$data
            true_labels = synthetic_data$labels
            
            # Conducts spectral clustering.
            est_labels = SpectralClustering(data, k)
            spec_err = ComputeErrorRate(est_labels, true_labels)
            res_tmp = c(res_tmp, spec_err)
        }
        res_sp = c(res_sp, mean(res_tmp))
    }
    
    data = data.frame(sig_ch = sig_ch, res_sp = res_sp)
    ggplot(data, aes(x = sig_ch, y = res_sp)) +
        geom_point(shape = 1) +
        xlab("signal strength") + ylab("error rate") +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14, face = "bold"))
    ggsave("figures/err_rate_sig.eps")
}

CheckErrorRateSignal()
