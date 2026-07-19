// Package spectralpde implements spectral methods for the numerical
// solution of differential equations.
//
// The package provides Fourier and Chebyshev collocation machinery:
// spectral differentiation matrices, Chebyshev-Gauss-Lobatto and periodic
// Fourier grids, discrete cosine/sine and (fast) Fourier transforms,
// barycentric interpolation, Clenshaw-Curtis and Gauss quadrature, and
// end-to-end solvers for the Poisson, Helmholtz, heat and
// advection-diffusion equations in one and two space dimensions.
//
// All routines are written against the Go standard library only. Grids and
// operators follow the conventions of Trefethen, Spectral Methods in MATLAB:
// the Chebyshev-Gauss-Lobatto nodes are x_j = cos(pi*j/N) for j = 0..N
// (ordered from +1 down to -1), and the periodic Fourier nodes are
// x_j = 2*pi*j/N for j = 0..N-1 on [0, 2*pi).
//
// Spectral methods converge exponentially fast for smooth data. The helper
// SpectralConvergenceRate and the error metrics (L2Error, LinfError, ...)
// make it easy to observe that behaviour in tests and applications.
//
// Matrices are stored row-major as [][]float64 unless wrapped by the Matrix
// type, and vectors are plain []float64. Functions never mutate their inputs
// unless explicitly documented to do so.
package spectralpde
