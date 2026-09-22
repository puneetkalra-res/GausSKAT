\name{acat_combine}
\alias{acat_combine}
\title{Aggregate Dependent P-values by ACAT}
\usage{
acat_combine(p_values, weights = NULL)
}
\arguments{
  \item{p_values}{Numeric vector of component p-values.}
  \item{weights}{Optional non-negative aggregation weights.}
}
\value{The aggregated Cauchy association test p-value.}
\description{
Combines dependent p-values using the aggregated Cauchy association test.
Equal weights are used by default.
}
