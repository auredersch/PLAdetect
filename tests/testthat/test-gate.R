test_that("dual gate finds only the cells that are high in both scores", {
  set.seed(1)
  plat <- c(rnorm(800, 0.05, 0.01), rnorm(100, 0.20, 0.02), rnorm(100, 0.20, 0.02))
  imm  <- c(rnorm(800, 0.10, 0.02), rnorm(100, 0.10, 0.02), rnorm(100, 0.30, 0.03))
  g <- gate_dual_gmm(plat, imm)
  expect_gt(mean(g$PLA[901:1000]), 0.9)    # real PLAs are found
  expect_lt(mean(g$PLA[1:900]), 0.05)      # everything else stays negative
})

test_that("gate_dual_gmm stops for too few cells", {
  expect_error(gate_dual_gmm(rnorm(10), rnorm(10)), "cells")
})