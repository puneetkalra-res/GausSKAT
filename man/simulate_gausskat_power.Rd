\name{simulate_gausskat_power}
\alias{simulate_gausskat_power}
\alias{print.GausSKAT_simulation}
\title{Simulate GausSKAT Power}
\usage{
simulate_gausskat_power(
  configuration = gausskat_simulation_config(),
  n_replications = 100L,
  seed = 7761L,
  haplotypes = NULL,
  comparison_kernels = c("linear.weighted", "IBS.weighted"),
  method = "davies",
  progress = interactive(),
  fail_fast = TRUE,
  parallel = FALSE,
  n_cores = NULL
)

\method{print}{GausSKAT_simulation}(x, ...)
}
\arguments{
  \item{configuration}{Configuration returned by
  \code{gausskat_simulation_config()}.}
  \item{n_replications}{Number of simulation replications.}
  \item{seed}{Non-negative integer random seed.}
  \item{haplotypes}{Optional haplotype matrix; the dataset supplied by
  \pkg{SKAT} is used by default.}
  \item{comparison_kernels}{SKAT kernels evaluated as comparators.}
  \item{method}{SKAT p-value method.}
  \item{progress}{Whether to print progress updates.}
  \item{fail_fast}{Whether the first failed replication stops the run.}
  \item{parallel}{Whether to evaluate the component Gaussian-kernel tests in
  parallel within each replication. Replications remain sequential so the
  seeded simulation sequence is unchanged.}
  \item{n_cores}{Number of component-test workers when \code{parallel = TRUE}.
  When \code{NULL}, at most one fewer than the detected logical cores is used,
  capped by the number of grid points.}
  \item{x}{A \code{GausSKAT_simulation} object.}
  \item{...}{Additional arguments reserved for print methods.}
}
\value{
An object of class \code{GausSKAT_simulation} containing replicate-level
results, summaries, the seed, configuration, and session information.
}
\description{
Runs a seed-controlled continuous-trait simulation and retains every grid
point and p-value needed for detailed supplementary summaries.
}
