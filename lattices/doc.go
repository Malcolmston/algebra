// Package lattices implements algorithms for integer and real lattices.
//
// A lattice is the set of all integer linear combinations of a set of
// linearly independent basis vectors b_0, ..., b_{m-1} living in R^n. This
// package represents a basis as an ordered list of row vectors ([Basis]) and
// provides the classical toolbox of lattice algorithms:
//
//   - Gram matrices, Gram determinants and lattice covolume (determinant).
//   - Gram-Schmidt orthogonalization with the mu coefficient matrix.
//   - LLL (Lenstra-Lenstra-Lovasz) basis reduction and standalone size
//     reduction.
//   - Shortest-vector (SVP) and closest-vector (CVP) search by
//     Fincke-Pohst enumeration.
//   - Babai's rounding and nearest-plane approximate CVP solvers.
//   - Minkowski and Hermite bounds and successive-minima estimates.
//   - The dual (reciprocal) lattice.
//   - Hermite normal form of integer matrices.
//
// Vectors and dense matrices come in three flavours: [Vec] and [Matrix] use
// float64 for speed and are used by the numeric reduction and enumeration
// routines; [RatVec] and [RatMatrix] use math/big.Rat for exact rational
// arithmetic (exact Gram determinants, exact inverses, exact duals); and
// [IntMatrix] uses math/big.Int for exact integer work such as the Hermite
// normal form. Everything is implemented with the Go standard library only.
//
// Conventions: bases are lists of rows, so a basis of m vectors in R^n is an
// m-element [Basis] whose entries are n-dimensional [Vec] values. The Gram
// matrix of a basis is G with G[i][j] = <b_i, b_j>, and the lattice
// determinant (covolume) is sqrt(det G).
package lattices
