# tensornetwork — implements dense N-dimensional real tensors together with the numerical machinery

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/tensornetwork.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/tensornetwork)

Package tensornetwork implements dense N-dimensional real tensors together
with the numerical machinery of tensor decompositions and tensor networks,
using only the Go standard library.

The package is fully self contained: it defines its own dense, row-major
tensor type (`Tensor`) and a small dense float64 linear-algebra backend
(`Matrix`) with singular-value, QR and symmetric-eigenvalue solvers, and
builds every higher-level algorithm on top of those primitives. Nothing
outside the standard library is imported and no cgo is used. Any algorithm
that needs randomness (for example the initialization of CP-ALS) takes an
explicit int64 seed and draws from a deterministic math/rand source, so
results are fully reproducible.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/tensornetwork"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func MultilinearRank(t *Tensor, tol float64) []int` | MultilinearRank returns the numerical multilinear rank of t: the rank of each mode-n unfolding, using the given relative singular-value tolerance. |

### Types

| Type | What it is |
| --- | --- |
| `CPDecomposition` | CPDecomposition is a rank-R CANDECOMP/PARAFAC factorization of a tensor. |
| `CPOptions` | CPOptions configures `CPALS`. |
| `ContractStep` | ContractStep is one pairwise contraction in an order: the two operand indices (referring to the running list of intermediate tensors) that are… |
| `Matrix` | Matrix is a dense, row-major matrix of float64 values. |
| `Network` | Network is a tensor network: a collection of tensors whose axes carry string labels. |
| `Tensor` | Tensor is a dense, row-major (C-ordered) N-dimensional array of float64 values. |
| `TensorTrain` | TensorTrain is a tensor-train (matrix-product-state) factorization. |
| `TuckerDecomposition` | TuckerDecomposition is a Tucker factorization: a core tensor together with one orthonormal factor matrix per mode. |

<details>
<summary><code>CPDecomposition</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CPALS(t *Tensor, rank int, opts CPOptions) (*CPDecomposition, error)` | CPALS computes a rank-R CP decomposition of t by alternating least squares. |
| `func (cp *CPDecomposition) Reconstruct() *Tensor` | Reconstruct returns the dense tensor represented by the CP decomposition. |
| `func (cp *CPDecomposition) RelError(t *Tensor) (float64, error)` | RelError returns the relative Frobenius reconstruction error of the CP decomposition against t. |

</details>

<details>
<summary><code>CPOptions</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DefaultCPOptions() CPOptions` | DefaultCPOptions returns reasonable defaults for `CPALS`: 200 iterations, a tolerance of 1e-10 and seed 1. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DiagMatrix(v []float64) *Matrix` | DiagMatrix returns a square diagonal matrix whose diagonal is v. |
| `func IdentityMatrix(n int) *Matrix` | IdentityMatrix returns the n x n identity matrix. |
| `func KhatriRao(mats ...*Matrix) (*Matrix, error)` | KhatriRao returns the Khatri-Rao product of a non-empty list of matrices, all with the same number of columns, combined left to right so that the… |
| `func KhatriRaoMatrix(a, b *Matrix) (*Matrix, error)` | KhatriRaoMatrix returns the column-wise Khatri-Rao product a⊙b. |
| `func KroneckerMatrix(a, b *Matrix) *Matrix` | KroneckerMatrix returns the Kronecker product a⊗b, an (a.rows*b.rows) x (a.cols*b.cols) matrix. |
| `func MatMul(a, b *Matrix) *Matrix` | MatMul returns a*b or panics on a dimension mismatch. |
| `func MatrixFromRows(rows [][]float64) (*Matrix, error)` | MatrixFromRows builds a matrix from a slice of equal-length rows. |
| `func NewMatrix(rows, cols int) *Matrix` | NewMatrix returns a new rows x cols matrix with all entries set to zero. |
| `func NewMatrixData(rows, cols int, data []float64) (*Matrix, error)` | NewMatrixData wraps the given row-major data as a rows x cols matrix. |
| `func RandMatrix(seed int64, rows, cols int) *Matrix` | RandMatrix returns a rows x cols matrix filled with independent samples from the standard normal distribution, drawn from a deterministic source… |
| `func (m *Matrix) Add(n *Matrix) (*Matrix, error)` | Add returns the entrywise sum m+n. |
| `func (m *Matrix) At(i, j int) float64` | At returns the entry at row i, column j. |
| `func (m *Matrix) Clone() *Matrix` | Clone returns an independent deep copy of m. |
| `func (m *Matrix) Col(j int) []float64` | Col returns a copy of column j. |
| `func (m *Matrix) Cols() int` | Cols returns the number of columns. |
| `func (m *Matrix) Cond() float64` | Cond returns the 2-norm condition number of m, the ratio of its largest to its smallest singular value. |
| `func (m *Matrix) Data() []float64` | Data returns the underlying row-major slice. |
| `func (m *Matrix) Diag() []float64` | Diag returns the main diagonal of m as a slice of length min(rows, cols). |
| `func (m *Matrix) EigSym() (values []float64, vectors *Matrix, err error)` | EigSym computes the eigenvalues and eigenvectors of a real symmetric matrix using the cyclic Jacobi method. |
| `func (m *Matrix) Equal(n *Matrix, tol float64) bool` | Equal reports whether m and n have the same shape and all corresponding entries differ by at most tol in absolute value. |
| `func (m *Matrix) FrobeniusNorm() float64` | FrobeniusNorm returns the Frobenius norm of m, sqrt of the sum of squares of all entries. |
| `func (m *Matrix) GramColumns() *Matrix` | GramColumns returns the Gram matrix mᵀ·m of the columns of m. |
| `func (m *Matrix) Hadamard(n *Matrix) (*Matrix, error)` | Hadamard returns the entrywise (Hadamard) product m∘n. |
| `func (m *Matrix) Mul(n *Matrix) (*Matrix, error)` | Mul returns the matrix product m*n. |
| `func (m *Matrix) MulVec(v []float64) ([]float64, error)` | MulVec returns the matrix-vector product m*v. |
| `func (m *Matrix) Pinv(tol float64) *Matrix` | Pinv returns the Moore-Penrose pseudoinverse of m computed from its SVD. |
| `func (m *Matrix) QR() (Q *Matrix, R *Matrix, err error)` | QR computes a reduced QR factorization of m using Householder reflections. |
| `func (m *Matrix) Rank(tol float64) int` | Rank returns the numerical rank of m, the number of singular values that exceed tol times the largest singular value. |
| `func (m *Matrix) Row(i int) []float64` | Row returns a copy of row i. |
| `func (m *Matrix) Rows() int` | Rows returns the number of rows. |
| `func (m *Matrix) SVD() (U *Matrix, s []float64, V *Matrix, err error)` | SVD computes a thin singular value decomposition of m using one-sided Jacobi rotations. |
| `func (m *Matrix) Scale(s float64) *Matrix` | Scale returns m with every entry multiplied by s. |
| `func (m *Matrix) Set(i, j int, v float64)` | Set stores v at row i, column j. |
| `func (m *Matrix) Solve(b []float64) ([]float64, error)` | Solve returns the least-squares solution x minimizing ‖m·x − b‖ using the pseudoinverse. |
| `func (m *Matrix) SpectralNorm() float64` | SpectralNorm returns the largest singular value of m. |
| `func (m *Matrix) Sub(n *Matrix) (*Matrix, error)` | Sub returns the entrywise difference m-n. |
| `func (m *Matrix) SubMatrix(r0, r1, c0, c1 int) *Matrix` | SubMatrix returns the r0:r1 x c0:c1 block of m (half-open ranges) as a fresh matrix. |
| `func (m *Matrix) Trace() float64` | Trace returns the sum of the diagonal entries of m. |
| `func (m *Matrix) Transpose() *Matrix` | Transpose returns the transpose of m. |

</details>

<details>
<summary><code>Network</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewNetwork() *Network` | NewNetwork returns an empty tensor network. |
| `func (nw *Network) Add(t *Tensor, labels ...string) error` | Add inserts a tensor into the network with the given axis labels. |
| `func (nw *Network) Contract() (*Tensor, []string, error)` | Contract evaluates the whole network in the greedy contraction order and returns the resulting tensor together with its (sorted) free labels. |
| `func (nw *Network) ContractOptimal() (*Tensor, []string, error)` | ContractOptimal contracts the network following the order from `Network.OptimalOrder` and returns the resulting tensor and its labels. |
| `func (nw *Network) ContractionCost(order []ContractStep) (float64, error)` | ContractionCost returns the total contraction cost of the given order under the same cost model as `Network.OptimalOrder`. |
| `func (nw *Network) FreeLabels() []string` | FreeLabels returns the sorted list of labels that appear on exactly one tensor: the free output indices of the network. |
| `func (nw *Network) GreedyOrder() ([]ContractStep, float64, error)` | GreedyOrder returns a contraction order chosen greedily: at each step it contracts the pair whose contraction cost (the product of the dimensions of… |
| `func (nw *Network) Len() int` | Len returns the number of tensors currently in the network. |
| `func (nw *Network) OptimalOrder() ([]ContractStep, float64, error)` | OptimalOrder returns a minimum-cost pairwise contraction order found by dynamic programming over subsets of tensors, together with the total cost. |

</details>

<details>
<summary><code>Tensor</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ARange(start, stop, step float64) (*Tensor, error)` | ARange returns a rank-1 tensor with values start, start+step, … up to but not including stop. |
| `func Concatenate(axis int, tensors ...*Tensor) (*Tensor, error)` | Concatenate joins tensors along the given axis. |
| `func Contract(a, b *Tensor, axisA, axisB int) (*Tensor, error)` | Contract contracts a single axis of a with a single axis of b, a convenience wrapper over `TensorDot`. |
| `func Einsum(spec string, operands ...*Tensor) (*Tensor, error)` | Einsum evaluates an Einstein-summation expression over the given operands. |
| `func Fold(m *Matrix, mode int, shape []int) (*Tensor, error)` | Fold reconstructs a tensor of the given shape from its mode-n matricization, inverting `Tensor.Unfold`. |
| `func FromMatrix(m *Matrix) *Tensor` | FromMatrix returns a rank-2 tensor holding a copy of the entries of m. |
| `func FromVector(v []float64) *Tensor` | FromVector returns a rank-1 tensor holding a copy of v. |
| `func Full(v float64, shape ...int) *Tensor` | Full returns a tensor of the given shape with every entry equal to v. |
| `func KhatriRaoTensors(tensors ...*Tensor) (*Tensor, error)` | KhatriRaoTensors returns the column-wise Khatri-Rao product of a list of rank-2 tensors interpreted as matrices, a convenience wrapper over… |
| `func Kronecker(a, b *Tensor) (*Tensor, error)` | Kronecker returns the Kronecker product of two tensors of equal rank. |
| `func LinSpace(start, stop float64, n int) (*Tensor, error)` | LinSpace returns a rank-1 tensor of n evenly spaced values from start to stop inclusive. |
| `func ModeProduct(t *Tensor, m *Matrix, mode int) (*Tensor, error)` | ModeProduct returns the n-mode product t ×ₙ m, contracting axis mode of t with the columns of m. |
| `func MultiModeProduct(t *Tensor, mats []*Matrix, modes []int) (*Tensor, error)` | MultiModeProduct applies a sequence of n-mode products, contracting each given mode of t with the corresponding matrix. |
| `func New(shape ...int) *Tensor` | New returns a new tensor of the given shape with all entries zero. |
| `func NewWithData(data []float64, shape ...int) (*Tensor, error)` | NewWithData returns a tensor of the given shape wrapping a copy of data, which must be in row-major order. |
| `func Ones(shape ...int) *Tensor` | Ones returns a tensor of the given shape with every entry equal to 1. |
| `func Outer(a, b *Tensor) *Tensor` | Outer returns the outer (tensor) product of a and b: a tensor whose shape is the concatenation of the shapes of a and b, with out[i…, j…] =… |
| `func RandTensor(seed int64, shape ...int) *Tensor` | RandTensor returns a tensor of the given shape filled with independent standard-normal samples from a deterministic source seeded by seed. |
| `func Scalar(v float64) *Tensor` | Scalar returns a rank-0 tensor holding the value v. |
| `func Stack(axis int, tensors ...*Tensor) (*Tensor, error)` | Stack joins tensors of identical shape along a new axis inserted at position axis. |
| `func TensorDot(a, b *Tensor, axesA, axesB []int) (*Tensor, error)` | TensorDot contracts a and b over the axis lists axesA and axesB, summing over the paired axes. |
| `func Zeros(shape ...int) *Tensor` | Zeros is an alias for `New`: it returns a zero-filled tensor of the given shape. |
| `func (t *Tensor) Abs() *Tensor` | Abs returns the entrywise absolute value of t. |
| `func (t *Tensor) Add(o *Tensor) (*Tensor, error)` | Add returns the entrywise sum t+o. |
| `func (t *Tensor) AddScalar(s float64) *Tensor` | AddScalar returns t with s added to every entry. |
| `func (t *Tensor) Apply(f func(float64) float64) *Tensor` | Apply returns a new tensor with f applied to every entry of t. |
| `func (t *Tensor) ArgMax() int` | ArgMax returns the flat index of the largest entry of t. |
| `func (t *Tensor) ArgMin() int` | ArgMin returns the flat index of the smallest entry of t. |
| `func (t *Tensor) At(idx ...int) float64` | At returns the element at the given multi-index. |
| `func (t *Tensor) AtErr(idx ...int) (float64, error)` | AtErr returns the element at the given multi-index or an error if the index is invalid. |
| `func (t *Tensor) AtFlat(i int) float64` | AtFlat returns the element at flat offset i in row-major order. |
| `func (t *Tensor) Clone() *Tensor` | Clone returns an independent deep copy of t. |
| `func (t *Tensor) Data() []float64` | Data returns the underlying row-major slice. |
| `func (t *Tensor) Dim(i int) int` | Dim returns the size of axis i. |
| `func (t *Tensor) Div(o *Tensor) (*Tensor, error)` | Div returns the entrywise quotient t/o. |
| `func (t *Tensor) Dot(o *Tensor) (float64, error)` | Dot returns the full inner product of t and o, the sum over all positions of the entrywise products. |
| `func (t *Tensor) Equal(o *Tensor, tol float64) bool` | Equal reports whether t and o have the same shape and all corresponding entries differ by at most tol in absolute value. |
| `func (t *Tensor) Exp() *Tensor` | Exp returns the entrywise exponential of t. |
| `func (t *Tensor) ExpandDims(axis int) (*Tensor, error)` | ExpandDims returns a tensor with a new axis of length 1 inserted at position axis. |
| `func (t *Tensor) Flatten() *Tensor` | Flatten is an alias for `Tensor.Ravel`. |
| `func (t *Tensor) Log() *Tensor` | Log returns the entrywise natural logarithm of t. |
| `func (t *Tensor) Max() float64` | Max returns the largest entry of t. |
| `func (t *Tensor) MaxAxis(axis int) (*Tensor, error)` | MaxAxis returns the maximum of t along the given axis, which is removed from the result. |
| `func (t *Tensor) Maximum(o *Tensor) (*Tensor, error)` | Maximum returns the entrywise maximum of t and o. |
| `func (t *Tensor) Mean() float64` | Mean returns the arithmetic mean of all entries of t. |
| `func (t *Tensor) MeanAxis(axis int) (*Tensor, error)` | MeanAxis returns the mean of t along the given axis, which is removed from the result. |
| `func (t *Tensor) Min() float64` | Min returns the smallest entry of t. |
| `func (t *Tensor) MinAxis(axis int) (*Tensor, error)` | MinAxis returns the minimum of t along the given axis, which is removed from the result. |
| `func (t *Tensor) Minimum(o *Tensor) (*Tensor, error)` | Minimum returns the entrywise minimum of t and o. |
| `func (t *Tensor) MoveAxis(src, dst int) (*Tensor, error)` | MoveAxis returns the tensor with the axis at position src moved to position dst, other axes keeping their relative order. |
| `func (t *Tensor) Mul(o *Tensor) (*Tensor, error)` | Mul returns the entrywise (Hadamard) product t∘o. |
| `func (t *Tensor) MultiIndex(flat int) []int` | MultiIndex converts a flat row-major offset into the corresponding multi-index. |
| `func (t *Tensor) Neg() *Tensor` | Neg returns the entrywise negation of t. |
| `func (t *Tensor) Norm() float64` | Norm returns the Frobenius (Euclidean) norm of t: the square root of the sum of squares of all entries. |
| `func (t *Tensor) NormP(p float64) float64` | NormP returns the entrywise p-norm of t for p >= 1. |
| `func (t *Tensor) Normalize() (*Tensor, float64)` | Normalize returns t scaled so that its Frobenius norm is 1, together with the original norm. |
| `func (t *Tensor) Permute(perm ...int) (*Tensor, error)` | Permute returns a tensor whose axes are reordered according to perm, a permutation of 0..rank-1. |
| `func (t *Tensor) ProdAxis(axis int) (*Tensor, error)` | ProdAxis returns t multiplied along the given axis, which is removed from the result. |
| `func (t *Tensor) Product() float64` | Product returns the product of all entries of t. |
| `func (t *Tensor) Rank() int` | Rank returns the number of dimensions (axes) of the tensor. |
| `func (t *Tensor) Ravel() *Tensor` | Ravel returns a rank-1 tensor with the same elements in row-major order. |
| `func (t *Tensor) RelError(o *Tensor) (float64, error)` | RelError returns the relative Frobenius error ‖t-o‖/‖t‖ between t and o. |
| `func (t *Tensor) Reshape(shape ...int) (*Tensor, error)` | Reshape returns a tensor sharing t's data reinterpreted with the new shape, which must have the same total number of elements. |
| `func (t *Tensor) Scale(s float64) *Tensor` | Scale returns t with every entry multiplied by s. |
| `func (t *Tensor) Set(v float64, idx ...int)` | Set stores v at the given multi-index. |
| `func (t *Tensor) SetFlat(i int, v float64)` | SetFlat stores v at flat offset i in row-major order. |
| `func (t *Tensor) Shape() []int` | Shape returns a copy of the tensor's shape. |
| `func (t *Tensor) Size() int` | Size returns the total number of elements in the tensor. |
| `func (t *Tensor) Sqrt() *Tensor` | Sqrt returns the entrywise square root of t. |
| `func (t *Tensor) Squeeze(axes ...int) (*Tensor, error)` | Squeeze returns a tensor with all length-1 axes removed. |
| `func (t *Tensor) Strides() []int` | Strides returns a copy of the tensor's row-major strides. |
| `func (t *Tensor) String() string` | String returns a short human-readable description of the tensor's shape and, for small tensors, its entries. |
| `func (t *Tensor) Sub(o *Tensor) (*Tensor, error)` | Sub returns the entrywise difference t-o. |
| `func (t *Tensor) Sum() float64` | Sum returns the sum of all entries of t. |
| `func (t *Tensor) SumAxis(axis int) (*Tensor, error)` | SumAxis returns t summed along the given axis, which is removed from the result. |
| `func (t *Tensor) SwapAxes(a, b int) (*Tensor, error)` | SwapAxes returns the tensor with axes a and b exchanged. |
| `func (t *Tensor) ToMatrix() (*Matrix, error)` | ToMatrix converts a rank-2 tensor into a `Matrix`. |
| `func (t *Tensor) TraceAxes(a, b int) (*Tensor, error)` | TraceAxes contracts axes a and b of a single tensor (a partial trace), summing over the diagonal where the two indices are equal and removing both… |
| `func (t *Tensor) Transpose() *Tensor` | Transpose returns the tensor with its axes reversed. |
| `func (t *Tensor) Unfold(mode int) (*Matrix, error)` | Unfold returns the mode-n matricization of t: a matrix whose rows are indexed by axis mode and whose columns are indexed by the remaining axes taken… |

</details>

<details>
<summary><code>TensorTrain</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func TTSVD(t *Tensor, eps float64) (*TensorTrain, error)` | TTSVD computes a tensor-train decomposition of t using successive truncated SVDs. |
| `func TTSVDRank(t *Tensor, maxRanks []int) (*TensorTrain, error)` | TTSVDRank computes a tensor-train decomposition of t whose internal ranks are capped by maxRanks (length t.Rank()-1). |
| `func (tt *TensorTrain) Compression() float64` | Compression returns the ratio of the number of dense entries to the number of parameters stored by the tensor train, a measure of how much the train… |
| `func (tt *TensorTrain) MaxRank() int` | MaxRank returns the largest internal rank of the tensor train. |
| `func (tt *TensorTrain) Reconstruct() *Tensor` | Reconstruct returns the dense tensor represented by the tensor train. |
| `func (tt *TensorTrain) RelError(t *Tensor) (float64, error)` | RelError returns the relative Frobenius reconstruction error of the tensor train against t. |

</details>

<details>
<summary><code>TuckerDecomposition</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func HOOI(t *Tensor, ranks []int, maxIter int) (*TuckerDecomposition, error)` | HOOI refines a Tucker decomposition of t by higher-order orthogonal iteration, starting from the HOSVD and alternately re-estimating each factor. |
| `func HOSVD(t *Tensor, ranks []int) (*TuckerDecomposition, error)` | HOSVD computes a truncated Tucker decomposition of t by the higher-order singular value decomposition. |
| `func (tu *TuckerDecomposition) Reconstruct() *Tensor` | Reconstruct returns the dense tensor represented by the Tucker decomposition. |
| `func (tu *TuckerDecomposition) RelError(t *Tensor) (float64, error)` | RelError returns the relative Frobenius reconstruction error of the Tucker decomposition against t. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/tensornetwork).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
