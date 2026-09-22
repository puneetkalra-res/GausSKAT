\name{gausskat_kernel}
\alias{gausskat_kernel}
\title{Construct a Weighted Gaussian Kernel}
\usage{
gausskat_kernel(D, ell)
}
\arguments{
  \item{D}{Matrix of weighted squared genetic distances.}
  \item{ell}{Positive Gaussian-kernel hyperparameter.}
}
\value{A weighted Gaussian kernel matrix.}
\description{
Returns the matrix whose \code{(i, i')} entry is
\code{exp(-D[i, i'] / ell)}.
}
