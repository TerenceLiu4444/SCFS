test_that("Lloyd returns correct value gaussian noise", {
  k=3
  s=10
  p=50
  n=100
  signal_strength=10
  synthetic_data = GenerateSyntheticData(n, p, s, k, signal_strength, noise_type="gaussian")
  init_cluster_ids = synthetic_data$labels
  init_cluster_ids[1:5] = 3
  label.est = LloydIteration(init_cluster_ids, num_iteration=10, data=synthetic_data$data, num_cluster=k)
  expect_equal(synthetic_data$labels, label.est)
})

test_that("Lloyd returns correct value t2 noise", {
  k=3
  s=10
  p=50
  n=100
  signal_strength=10
  synthetic_data = GenerateSyntheticData(n, p, s, k, signal_strength, noise_type="t2")
  init_cluster_ids = synthetic_data$labels
  init_cluster_ids[1:5] = 3
  label.est = LloydIteration(init_cluster_ids, num_iteration=10, data=synthetic_data$data, num_cluster=k)
  expect_equal(synthetic_data$labels, label.est)
})

test_that("Spectral clustering", {
  k=3
  s=10
  p=50
  n=100
  signal_strength=100
  synthetic_data = GenerateSyntheticData(n, p, s, k, signal_strength, noise_type="gaussian")
  label.est = SpectralClustering(data=synthetic_data$data, num_cluster=k)
  expect_lt(ComputeErrorRate(label.est, synthetic_data$labels), 0.01)
})

test_that("SCFS", {
  k=3
  s=10
  p=50
  n=100
  signal_strength=20
  synthetic_data = GenerateSyntheticData(n, p, s, k, signal_strength, noise_type="gaussian")
  scfs = SpectralClusterFeatureSelection(data=synthetic_data$data, num_cluster=k, init_cluster_ids = NULL, use_lloyd_iteration=FALSE)
  expect_lt(ComputeErrorRate(scfs$cluster_ids, synthetic_data$labels), 0.01)
})

test_that("Error rate computation", {
  expect_equal(ComputeErrorRate(c(0,1), c(0,1)), 0)
  expect_equal(ComputeErrorRate(c(0,1), c(1,0)), 0)
  expect_equal(ComputeErrorRate(c(0,1,1,1), c(0,0,1,1)), 0.25)
})


