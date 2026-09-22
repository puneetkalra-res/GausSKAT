test_that("ACAT returns the component p-value when all components agree", {
  for (p in c(1e-12, 0.01, 0.25, 0.9)) {
    expect_equal(acat_combine(rep(p, 5)), p, tolerance = 1e-10)
  }
})

test_that("ACAT validates inputs", {
  expect_error(acat_combine(numeric()), "At least one")
  expect_error(acat_combine(c(-0.1, 0.2)), "between zero and one")
  expect_error(acat_combine(c(0.1, 0.2), weights = 1), "same length")
  expect_true(is.na(acat_combine(c(0.1, NA_real_))))
})
