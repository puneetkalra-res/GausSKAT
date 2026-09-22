#' Aggregate dependent p-values by ACAT
#'
#' @param p_values Numeric vector of component p-values.
#' @param weights Optional non-negative aggregation weights. Equal weights are
#'   used by default.
#' @return The aggregated Cauchy association test p-value.
#' @export
acat_combine <- function(p_values, weights = NULL) {
  p_values <- as.numeric(p_values)
  if (length(p_values) == 0L) {
    stop("At least one p-value is required.", call. = FALSE)
  }
  if (any(!is.finite(p_values))) {
    return(NA_real_)
  }
  if (any(p_values < 0 | p_values > 1)) {
    stop("All p-values must lie between zero and one.", call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- rep(1 / length(p_values), length(p_values))
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != length(p_values)) {
      stop("weights and p_values must have the same length.", call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
      stop("weights must be finite, non-negative, and have a positive sum.",
           call. = FALSE)
    }
    weights <- weights / sum(weights)
  }

  active <- weights > 0
  p_values <- p_values[active]
  weights <- weights[active]
  weights <- weights / sum(weights)
  if (any(p_values == 0)) {
    return(0)
  }

  p_work <- pmin(p_values, 1 - .Machine$double.eps)
  very_small <- p_work < 1e-15
  terms <- numeric(length(p_work))
  terms[very_small] <- 1 / (pi * p_work[very_small])
  terms[!very_small] <- tan((0.5 - p_work[!very_small]) * pi)
  statistic <- sum(weights * terms)

  if (statistic > 1e15) {
    p_value <- 1 / (pi * statistic)
  } else {
    p_value <- stats::pcauchy(statistic, lower.tail = FALSE)
  }
  min(max(p_value, 0), 1)
}
