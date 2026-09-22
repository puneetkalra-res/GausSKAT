#!/usr/bin/env Rscript

library(SKAT)
library(GausSKAT)

data("SKAT.example", package = "SKAT")
null_model <- SKAT_Null_Model(
  y.c ~ X,
  data = SKAT.example,
  out_type = "C"
)

fit <- GausSKAT(
  Z = SKAT.example$Z,
  null_model = null_model,
  epsilon = 0.05,
  number_grid_points = 5L
)

print(fit)
print(data.frame(
  grid_position = seq_along(fit$ell.grid),
  ell = fit$ell.grid,
  p_value = fit$component.p.values
))
