\name{make_adaptive_grid}
\alias{make_adaptive_grid}
\title{Construct the Data-Adaptive GausSKAT Grid}
\usage{
make_adaptive_grid(D, X, epsilon = 0.05, number_grid_points = 5L)
}
\arguments{
  \item{D}{Matrix of weighted squared genetic distances.}
  \item{X}{Null-model design matrix whose column space contains an intercept.}
  \item{epsilon}{Positive stabilization tolerance smaller than
  \code{sqrt(2)}.}
  \item{number_grid_points}{Number of logarithmically spaced grid points.}
}
\value{
A list containing \code{ell_min}, \code{ell_max_approx}, the grid, and
projected-kernel diagnostics.
}
\description{
Uses the median strictly positive distance for the lower endpoint and the
first-order, covariate-projected stabilization approximation for the practical
upper endpoint.
}
