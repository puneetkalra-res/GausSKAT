\name{gausskat_simulation_config}
\alias{gausskat_simulation_config}
\title{Define a GausSKAT Simulation Configuration}
\usage{
gausskat_simulation_config(
  sample_size = 1000L,
  region_length = 200L,
  causal_percent = 50,
  causal_maf_cutoff = 0.01,
  negative_percent = 0,
  effect_size_constant = 0.15,
  alpha = 1e-06,
  weights_beta = c(1, 25),
  epsilon = 0.05,
  number_grid_points = 5L,
  covariate_effects = c(0.5, 0.5),
  residual_sd = 1,
  haplotype_sampling = c("manuscript", "full_pool")
)
}
\arguments{
  \item{sample_size}{Number of diploid individuals.}
  \item{region_length}{Number of consecutive variants per region.}
  \item{causal_percent}{Percentage of eligible variants that are causal.}
  \item{causal_maf_cutoff}{Upper MAF cutoff for candidate causal variants.}
  \item{negative_percent}{Percentage of causal effects with negative signs.}
  \item{effect_size_constant}{Constant in
  \code{abs(log10(MAF)) * effect_size_constant}.}
  \item{alpha}{Significance threshold used to estimate power.}
  \item{weights_beta}{Two beta-weight parameters.}
  \item{epsilon}{Upper-endpoint stabilization tolerance.}
  \item{number_grid_points}{Number of GausSKAT grid points.}
  \item{covariate_effects}{Effects of the continuous and binary covariates.}
  \item{residual_sd}{Standard deviation of the normal error.}
  \item{haplotype_sampling}{Haplotype sampling scheme. The default reproduces
  the manuscript research scripts; \code{"full_pool"} samples from all
  available reference haplotypes.}
}
\value{A validated object of class \code{GausSKAT_simulation_config}.}
\description{Defines the parameters for a reproducible power simulation.}
