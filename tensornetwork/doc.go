// Package tensornetwork implements dense N-dimensional real tensors together
// with the numerical machinery of tensor decompositions and tensor networks,
// using only the Go standard library.
//
// The package is fully self contained: it defines its own dense, row-major
// tensor type ([Tensor]) and a small dense float64 linear-algebra backend
// ([Matrix]) with singular-value, QR and symmetric-eigenvalue solvers, and
// builds every higher-level algorithm on top of those primitives. Nothing
// outside the standard library is imported and no cgo is used. Any algorithm
// that needs randomness (for example the initialization of CP-ALS) takes an
// explicit int64 seed and draws from a deterministic math/rand source, so
// results are fully reproducible.
//
// # Tensors
//
// A [Tensor] is a rank-r array of float64 values stored contiguously in
// C (row-major) order. The package covers construction ([New], [Zeros],
// [Ones], [Full], [NewWithData], [Scalar], [FromVector], [FromMatrix],
// [ARange], [LinSpace], [RandTensor]), inspection ([Tensor.Shape],
// [Tensor.Rank], [Tensor.Size], [Tensor.At], [Tensor.Set]), shape
// manipulation ([Tensor.Reshape], [Tensor.Permute], [Tensor.Transpose],
// [Tensor.SwapAxes], [Tensor.MoveAxis], [Tensor.ExpandDims],
// [Tensor.Squeeze], [Concatenate], [Stack]) and elementwise arithmetic and
// reductions ([Tensor.Add], [Tensor.Mul], [Tensor.Scale], [Tensor.Sum],
// [Tensor.Norm], [Tensor.SumAxis]).
//
// # Products and contractions
//
// Multilinear products include the outer product ([Outer]), the Kronecker and
// Khatri-Rao products ([Kronecker], [KhatriRao], [KroneckerMatrix],
// [KhatriRaoMatrix]), the n-mode product with a matrix ([ModeProduct],
// [MultiModeProduct]), general tensor contraction ([TensorDot], [Contract]),
// and a small Einstein-summation evaluator ([Einsum]). Matricization is
// provided by [Tensor.Unfold] and [Fold].
//
// # Decompositions
//
// The package implements the three workhorse tensor factorizations:
//
//   - CP / PARAFAC via alternating least squares ([CPALS], [CPDecomposition]);
//   - Tucker via the higher-order SVD and higher-order orthogonal iteration
//     ([HOSVD], [HOOI], [TuckerDecomposition]);
//   - tensor-train / matrix-product-state factorization ([TTSVD],
//     [TTSVDRank], [TensorTrain]).
//
// Each decomposition can reconstruct a dense approximation and report the
// relative reconstruction error.
//
// # Matrix backend
//
// [Matrix] is a dense row-major float64 matrix with the operations needed by
// the decompositions: [Matrix.Mul], [Matrix.Transpose], [Matrix.SVD],
// [Matrix.QR], [Matrix.EigSym], [Matrix.Pinv], [Matrix.Rank] and
// [Matrix.Solve]. The SVD uses a one-sided Jacobi sweep and the symmetric
// eigensolver a cyclic Jacobi sweep; both are accurate for the small to
// moderate matrices that arise inside tensor decompositions.
//
// # Network contraction
//
// A [Network] describes a set of tensors sharing labelled indices. Its total
// contraction can be evaluated ([Network.Contract]) and the cheapest pairwise
// contraction order found by dynamic programming ([Network.OptimalOrder]) or
// approximated greedily ([Network.GreedyOrder]).
package tensornetwork
