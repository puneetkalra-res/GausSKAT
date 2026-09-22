#!/usr/bin/env Rscript

# Verifies that the independent distance implementation reproduces the exact
# weighted-Gaussian formula used in the original research code:
#   weighted_Z <- sweep(Z, 2, weights, "*")
#   K <- exp(-as.matrix(dist(weighted_Z))^2 / ell)

library(GausSKAT)

Z <- rbind(
  c(0, 0, 1, 2),
  c(0, 1, 1, 1),
  c(1, 0, 2, 0),
  c(2, 1, 0, 0)
)
weights <- c(25, 10, 3, 1)
ell <- 100

reference_weighted_Z <- sweep(Z, 2L, weights, `*`)
reference_kernel <- exp(
  -(as.matrix(stats::dist(reference_weighted_Z))^2) / ell
)

D <- weighted_squared_distances(Z, weights)
independent_kernel <- gausskat_kernel(D, ell)

stopifnot(isTRUE(all.equal(
  independent_kernel,
  reference_kernel,
  tolerance = 1e-12
)))
message("Weighted-Gaussian kernel equivalence check passed.")
