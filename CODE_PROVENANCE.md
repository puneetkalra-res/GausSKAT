# Code provenance and separation from SKAT

The research prototype was developed by importing the SKAT source into an
RStudio project and adding a weighted Gaussian kernel with a user-supplied
hyperparameter. The three changes essential to that prototype were:

1. construction of `weighted_Z` by multiplying genotype column `k` by the
   supplied SKAT weight;
2. construction of `K = exp(-dist(weighted_Z)^2 / ell)`; and
3. propagation of `ell` through the local SKAT call stack.

This repository implements the same kernel formula externally and passes the
resulting matrix to the official `SKAT` package. The official package already
supports a user-supplied kernel matrix, so the second and third source-code
changes are not needed in the independent implementation.

The repository deliberately excludes all other changes present in the local
research copy, including alternative binary-trait links, GPU experiments,
package renaming, altered default kernels, obsolete endpoint searches, and
local compiled objects. It also excludes all participant-level Dallas Heart
Study files and intermediate simulation output.

The script `scripts/verify_kernel_equivalence.R` compares the independent
kernel constructor with the prototype formula on a fixed example. Automated
tests additionally cover endpoint construction, weight-scale equivariance,
the large-hyperparameter weighted-linear limit, the SKAT interface, ACAT, and
the simulation output structure.

The default simulation RNG sequence follows `revision3W_simulation.R`, the
final full-simulation script: a fresh-session Mersenne--Twister stream is seeded
once and used continuously across replications. The normal and sampling
generators are fixed explicitly as `Inversion` and `Rejection` so that prior
session settings cannot silently change the sequence.
