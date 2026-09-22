\name{weighted_squared_distances}
\alias{weighted_squared_distances}
\title{Weighted Squared Genetic Distances}
\usage{
weighted_squared_distances(Z, weights)
}
\arguments{
  \item{Z}{Complete genotype dosage matrix with individuals in rows.}
  \item{weights}{Non-negative per-variant weights supplied to SKAT.}
}
\value{A symmetric matrix of weighted squared genetic distances.}
\description{
Computes pairwise squared Euclidean distances after multiplying genotype
column \code{k} by \code{weights[k]}. The effective distance coefficient is
therefore \code{weights[k]^2}.
}
