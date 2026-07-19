// Package randommatrix implements the core objects of random matrix theory
// (RMT) using only the Go standard library.
//
// Random matrix theory studies the statistical behaviour of the eigenvalues of
// matrices whose entries are drawn at random. Despite the randomness of the
// individual entries, the collective spectrum obeys sharp, universal laws. This
// package provides seeded generators for the classical Gaussian and Wishart
// ensembles together with the analytic machinery used to describe their spectra.
//
// The package is organised around a small number of themes:
//
//   - Dense real and complex matrix types (Matrix and CMatrix) with the linear
//     algebra needed to extract spectra: a Jacobi eigensolver for real
//     symmetric matrices and a Hermitian eigensolver built on the real
//     symmetric embedding.
//
//   - Seeded generators for the Gaussian ensembles GOE (beta = 1), GUE
//     (beta = 2) and GSE (beta = 4), for the Wishart / Laguerre ensembles, and
//     for the Ginibre ensembles. Every generator takes an explicit int64 seed
//     and uses math/rand, so results are reproducible and never depend on the
//     wall clock or on a cryptographic source.
//
//   - The limiting spectral laws: the Wigner semicircle law for Wigner matrices
//     and the Marchenko-Pastur law for sample covariance (Wishart) matrices,
//     with densities, cumulative distributions, supports and moments.
//
//   - Edge fluctuations described by the Tracy-Widom distributions for
//     beta = 1, 2, 4, provided through the accurate shifted-gamma approximation
//     of Chiani (2014) together with tabulated moments and the soft-edge
//     rescaling of the largest eigenvalue.
//
//   - Bulk fluctuations described by the nearest-neighbour level-spacing
//     distributions: the Poisson law for uncorrelated spectra and the Wigner
//     surmise for the Gaussian ensembles, plus the level-spacing ratio
//     statistics that avoid the need for spectral unfolding.
//
//   - Spectral moments and elementary free probability: the Stieltjes / Cauchy
//     transform, the R-transform, free cumulants via non-crossing partitions
//     and free additive convolution of the semicircle law.
//
// All randomness is seeded by the caller. Only the packages math, math/big,
// math/cmplx, sort, errors, fmt and strings are imported; there is no cgo and
// no third-party dependency.
package randommatrix
