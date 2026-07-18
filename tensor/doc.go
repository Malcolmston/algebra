// Package tensor implements dense, row-major N-dimensional real tensors and the
// core operations of multilinear algebra on them.
//
// The central type is [Tensor], a rank-r array of float64 values stored in a
// single contiguous slice using C (row-major) ordering. Rank 0 tensors are
// scalars, rank 1 are vectors and rank 2 are matrices; the same code paths
// handle arbitrary rank.
//
// The package covers:
//
//   - construction and inspection ([New], [Zeros], [Ones], [Full],
//     [NewWithData], [FromVector], [FromMatrix], [Tensor.At], [Tensor.Set],
//     [Tensor.Shape], [Tensor.Rank], [Tensor.Size]);
//   - shape manipulation ([Tensor.Reshape], [Tensor.Ravel], [Tensor.Transpose],
//     [Tensor.Permute], [Tensor.SwapAxes], [Tensor.MoveAxis],
//     [Tensor.ExpandDims], [Tensor.Squeeze], [Concatenate], [Stack]);
//   - elementwise arithmetic and reductions ([Tensor.Add], [Tensor.Mul],
//     [Tensor.Scale], [Tensor.Apply], [Tensor.Sum], [Tensor.Norm],
//     [Tensor.SumAxis]);
//   - multilinear products and contractions ([Outer], [Kronecker], [TensorDot],
//     [MatMul], [Dot], [Tensor.Contract], [Tensor.Trace], [Einsum]);
//   - index gymnastics on (pseudo-)Riemannian tensors ([KroneckerDelta],
//     [LeviCivita], [LeviCivitaTensor], [EuclideanMetric], [MinkowskiMetric],
//     [LowerIndex], [RaiseIndex]).
//
// All computation is deterministic and uses only the Go standard library.
// Methods that cannot express a failure in their signature (for example
// [Tensor.At]) panic on misuse; every other operation returns a sentinel error
// such as [ErrShape], [ErrAxis] or [ErrIndex].
package tensor
