// Package operatortheory implements finite-dimensional linear operator and
// spectral theory using only the Go standard library.
//
// The central type is [Matrix], a dense complex matrix that is interpreted as a
// bounded linear operator on the finite-dimensional Hilbert space C^n endowed
// with the standard inner product <x,y> = sum conj(x_i) y_i. On top of the
// usual matrix arithmetic (addition, multiplication, adjoint, Kronecker
// product, powers) the package provides the machinery of operator theory:
//
//   - Classification predicates: [Matrix.IsHermitian], [Matrix.IsSkewHermitian],
//     [Matrix.IsUnitary], [Matrix.IsNormal], [Matrix.IsProjection],
//     [Matrix.IsPositiveDefinite], [Matrix.IsContraction], [Matrix.IsIsometry]
//     and many more.
//   - Norms: the Frobenius, induced 1-, 2- (spectral) and infinity norms, the
//     nuclear (trace) norm and the numerical radius.
//   - Spectral theory: [Matrix.Eigenvalues], [Matrix.HermitianEigen],
//     [Matrix.Spectrum], [Matrix.SpectralRadius], the [Matrix.Resolvent],
//     the [Matrix.NumericalRange] and [Matrix.NumericalAbscissa].
//   - Decompositions: [Matrix.QR], [Matrix.Hessenberg], [Matrix.SVD],
//     [Matrix.PolarDecomposition] and the Hermitian spectral decomposition
//     [Matrix.SpectralDecomposition].
//   - Functional calculus on normal operators: [Matrix.ApplyFunction] together
//     with the ready-made [Matrix.Exp], [Matrix.Log], [Matrix.Sqrt],
//     [Matrix.Sign] and trigonometric functions.
//   - Pseudospectra estimates: [Matrix.ResolventNorm],
//     [Matrix.PseudospectralAbscissa] and [Matrix.PseudospectrumGrid].
//
// Numerical methods. Symmetric/Hermitian eigenproblems are solved with the
// cyclic Jacobi method applied to the real 2n-by-2n symmetric embedding of a
// Hermitian matrix, which is backward stable and returns a full orthonormal set
// of eigenvectors. General (non-normal) spectra are computed with the
// explicitly shifted QR algorithm on the upper-Hessenberg form. The singular
// value decomposition is obtained from the Hermitian eigendecomposition of the
// Gram matrix. Every routine depends only on math, math/cmplx and sort.
//
// Randomised constructors ([RandomMatrix], [RandomHermitian], [RandomUnitary])
// take a caller-supplied seed and use math/rand so that results are fully
// reproducible.
package operatortheory
