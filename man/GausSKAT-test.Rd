\name{GausSKAT-test}
\alias{GausSKAT-test}
\alias{GausSKAT}
\alias{print.GausSKAT}
\title{GausSKAT Test for a Continuous Trait}
\usage{
GausSKAT(
  Z,
  null_model,
  X = NULL,
  weights = NULL,
  weights_beta = c(1, 25),
  epsilon = 0.05,
  number_grid_points = 5L,
  method = "davies",
  acat_weights = NULL,
  keep_component_fits = FALSE,
  warn_on_endpoint = TRUE,
  parallel = FALSE,
  n_cores = NULL
)

\method{print}{GausSKAT}(x, ...)
}
\arguments{
  \item{Z}{Complete genotype dosage matrix with individuals in rows.}
  \item{null_model}{Continuous-trait null model returned by
  \code{SKAT::SKAT_Null_Model()}.}
  \item{X}{Optional null-model design matrix. The matrix stored in the null
  model is used by default.}
  \item{weights}{Optional non-negative per-variant weights.}
  \item{weights_beta}{Two beta-weight parameters used when \code{weights} is
  \code{NULL}.}
  \item{epsilon}{Positive stabilization tolerance smaller than
  \code{sqrt(2)}.}
  \item{number_grid_points}{Number of logarithmically spaced kernels.}
  \item{method}{SKAT p-value method: \code{"davies"}, \code{"liu"}, or
  \code{"liu.mod"}.}
  \item{acat_weights}{Optional ACAT component weights.}
  \item{keep_component_fits}{Whether to retain the complete component fits.}
  \item{warn_on_endpoint}{Whether to warn when the observed change at the
  approximate upper endpoint exceeds \code{epsilon}.}
  \item{parallel}{Whether to evaluate the component Gaussian-kernel tests in
  parallel. The default is \code{FALSE}.}
  \item{n_cores}{Number of component-test workers when \code{parallel = TRUE}.
  When \code{NULL}, at most one fewer than the detected logical cores is used,
  capped by the number of grid points.}
  \item{x}{A fitted \code{GausSKAT} object.}
  \item{...}{Additional arguments reserved for print methods.}
}
\value{
An object of class \code{GausSKAT} containing the ACAT p-value, component
p-values, adaptive grid, endpoint estimates, and diagnostics.
}
\description{
Constructs the phenotype-independent adaptive weighted-Gaussian grid,
evaluates every component using \pkg{SKAT}, and combines the component
p-values with ACAT.
}
\examples{
\dontrun{
data("SKAT.example", package = "SKAT")
obj <- SKAT::SKAT_Null_Model(y.c ~ X, data = SKAT.example, out_type = "C")
GausSKAT(SKAT.example$Z, obj)
}
}
