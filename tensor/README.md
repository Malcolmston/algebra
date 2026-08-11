# tensor — implements dense, row-major N-dimensional real tensors and the core operations

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/tensor.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/tensor)

Package tensor implements dense, row-major N-dimensional real tensors and the
core operations of multilinear algebra on them.

The central type is `Tensor`, a rank-r array of float64 values stored in a
single contiguous slice using C (row-major) ordering. Rank 0 tensors are
scalars, rank 1 are vectors and rank 2 are matrices; the same code paths
handle arbitrary rank.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/tensor"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Dot(a, b *Tensor) (float64, error)` | Dot returns the inner product of two rank-1 tensors of equal length. |
| `func KroneckerDelta(i, j int) float64` | KroneckerDelta returns the Kronecker delta δ_ij: 1 when i == j and 0 otherwise. |
| `func LeviCivita(indices ...int) int` | LeviCivita returns the value of the Levi-Civita symbol ε for the given indices. |
| `func RavelIndex(shape, idx []int) (int, error)` | RavelIndex converts a multi-index into a flat row-major offset for the given shape. |
| `func UnravelIndex(shape []int, flat int) ([]int, error)` | UnravelIndex converts a flat row-major offset into the corresponding multi-index for the given shape. |

### Types

| Type | What it is |
| --- | --- |
| `Tensor` | Tensor is a dense, row-major, rank-r array of float64 values. |

<details>
<summary><code>Tensor</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Concatenate(axis int, tensors ...*Tensor) (*Tensor, error)` | Concatenate joins tensors along an existing axis. |
| `func Einsum(spec string, tensors ...*Tensor) (*Tensor, error)` | Einsum evaluates an Einstein-summation expression over the given tensors. |
| `func EuclideanMetric(n int) *Tensor` | EuclideanMetric returns the n×n Euclidean metric tensor, i.e. |
| `func FromMatrix(m [][]float64) (*Tensor, error)` | FromMatrix returns a rank-2 tensor holding a copy of the row-major matrix m. |
| `func FromScalar(v float64) *Tensor` | FromScalar returns a rank-0 tensor holding the single value v. |
| `func FromVector(v []float64) *Tensor` | FromVector returns a rank-1 tensor holding a copy of v. |
| `func Full(value float64, shape ...int) *Tensor` | Full returns a tensor of the given shape with every element equal to value. |
| `func Identity(n int) *Tensor` | Identity returns the n×n identity matrix as a rank-2 tensor (the components of the Kronecker delta δ_ij). |
| `func Kronecker(a, b *Tensor) (*Tensor, error)` | Kronecker returns the Kronecker product of two tensors of equal rank. |
| `func KroneckerDeltaTensor(n int) *Tensor` | KroneckerDeltaTensor is a synonym for `Identity`: the n×n matrix of components δ_ij. |
| `func LeviCivitaTensor(n int) *Tensor` | LeviCivitaTensor returns the rank-n Levi-Civita (permutation) tensor in n dimensions: a tensor of shape [n,n,...,n] (n axes) whose component at a… |
| `func LowerIndex(t, g *Tensor, axis int) (*Tensor, error)` | LowerIndex lowers the index on the given axis of t using the metric g, returning a tensor of the same shape whose axis component is t'_{...a...} =… |
| `func MatMul(a, b *Tensor) (*Tensor, error)` | MatMul returns the matrix product of two rank-2 tensors a (m×k) and b (k×n), giving an m×n tensor. |
| `func MinkowskiMetric(dim int) *Tensor` | MinkowskiMetric returns the dim×dim Minkowski metric tensor with the mostly-plus signature (-,+,+,...): component (0,0) is -1 and the remaining… |
| `func New(shape ...int) *Tensor` | New returns a new zero-filled tensor with the given shape. |
| `func NewWithData(shape []int, data []float64) (*Tensor, error)` | NewWithData returns a tensor with the given shape whose elements are the row-major contents of data. |
| `func Ones(shape ...int) *Tensor` | Ones returns a tensor of the given shape with every element equal to 1. |
| `func Outer(a, b *Tensor) *Tensor` | Outer returns the tensor (outer) product a⊗b. |
| `func RaiseIndex(t, gInv *Tensor, axis int) (*Tensor, error)` | RaiseIndex raises the index on the given axis of t using the inverse metric gInv, returning a tensor of the same shape whose axis component is… |
| `func Stack(axis int, tensors ...*Tensor) (*Tensor, error)` | Stack joins tensors along a new axis inserted at position axis, increasing the rank by one. |
| `func TensorDot(a, b *Tensor, axesA, axesB []int) (*Tensor, error)` | TensorDot contracts a and b over the paired axes in axesA and axesB. For each i the axis axesA`i` of a is summed against axis axesB`i` of b, which… |
| `func TensorProduct(a, b *Tensor) *Tensor` | TensorProduct is a synonym for `Outer`, using the name common in multilinear-algebra texts for the product a⊗b. |
| `func Zeros(shape ...int) *Tensor` | Zeros returns a tensor of the given shape filled with zeros. |
| `func (t *Tensor) Abs() *Tensor` | Abs returns a tensor holding the absolute value of every element of t. |
| `func (t *Tensor) Add(other *Tensor) (*Tensor, error)` | Add returns the elementwise sum t+other. |
| `func (t *Tensor) AddScalar(s float64) *Tensor` | AddScalar returns a tensor with s added to every element of t. |
| `func (t *Tensor) AlmostEqual(other *Tensor, tol float64) bool` | AlmostEqual reports whether t and other have the same shape and every pair of corresponding elements differs by at most tol in absolute value. |
| `func (t *Tensor) Apply(f func(float64) float64) *Tensor` | Apply returns a tensor obtained by applying f to every element of t. |
| `func (t *Tensor) At(idx ...int) float64` | At returns the element of t at the given multi-index. |
| `func (t *Tensor) Clone() *Tensor` | Clone returns an independent deep copy of t. |
| `func (t *Tensor) Contract(axisA, axisB int) (*Tensor, error)` | Contract returns the tensor obtained by summing t over a pair of its own axes axisA and axisB, which must be distinct and have equal length. |
| `func (t *Tensor) Data() []float64` | Data returns a copy of t's underlying elements in row-major order. |
| `func (t *Tensor) Dim(axis int) int` | Dim returns the length of the given axis. |
| `func (t *Tensor) Div(other *Tensor) (*Tensor, error)` | Div returns the elementwise quotient t/other. |
| `func (t *Tensor) Equal(other *Tensor) bool` | Equal reports whether t and other have the same shape and bit-for-bit equal elements. |
| `func (t *Tensor) ExpandDims(axis int) (*Tensor, error)` | ExpandDims returns a new tensor with a length-1 axis inserted at the given position, increasing the rank by one. |
| `func (t *Tensor) Flatten() []float64` | Flatten returns a copy of t's elements in row-major order as a plain slice. |
| `func (t *Tensor) IsScalar() bool` | IsScalar reports whether t is a rank-0 tensor. |
| `func (t *Tensor) Max() float64` | Max returns the largest element of t. |
| `func (t *Tensor) Mean() float64` | Mean returns the arithmetic mean of all elements of t. |
| `func (t *Tensor) Min() float64` | Min returns the smallest element of t. |
| `func (t *Tensor) MoveAxis(src, dst int) (*Tensor, error)` | MoveAxis returns a new tensor with the axis at position src moved to position dst, the remaining axes keeping their relative order. |
| `func (t *Tensor) Mul(other *Tensor) (*Tensor, error)` | Mul returns the elementwise (Hadamard) product t*other. |
| `func (t *Tensor) Neg() *Tensor` | Neg returns a tensor with every element of t negated. |
| `func (t *Tensor) Norm() float64` | Norm returns the Frobenius (Euclidean) norm of t, the square root of the sum of squared elements. |
| `func (t *Tensor) Permute(perm ...int) (*Tensor, error)` | Permute returns a new tensor whose axes are a reordering of t's axes: axis i of the result is axis perm`i` of t. |
| `func (t *Tensor) Product() float64` | Product returns the product of all elements of t. |
| `func (t *Tensor) Rank() int` | Rank returns the number of axes of t (0 for a scalar). |
| `func (t *Tensor) Ravel() *Tensor` | Ravel returns a rank-1 tensor holding t's elements in row-major order. |
| `func (t *Tensor) Reshape(shape ...int) (*Tensor, error)` | Reshape returns a new tensor with the given shape and the same elements in row-major order. |
| `func (t *Tensor) ScalarValue() (float64, error)` | ScalarValue returns the single element of a rank-0 tensor. |
| `func (t *Tensor) Scale(s float64) *Tensor` | Scale returns a tensor with every element of t multiplied by s. |
| `func (t *Tensor) Set(value float64, idx ...int)` | Set stores value at the given multi-index of t, following the same index rules as `Tensor.At`. |
| `func (t *Tensor) Shape() []int` | Shape returns a copy of t's dimensions, one entry per axis. |
| `func (t *Tensor) ShapeEqual(other *Tensor) bool` | ShapeEqual reports whether t and other have identical shapes. |
| `func (t *Tensor) Size() int` | Size returns the total number of elements in t (1 for a scalar). |
| `func (t *Tensor) Squeeze() *Tensor` | Squeeze returns a new tensor with every length-1 axis removed. |
| `func (t *Tensor) Strides() []int` | Strides returns a copy of t's row-major strides, one entry per axis. |
| `func (t *Tensor) String() string` | String returns a compact, deterministic textual description of t giving its shape and, for small tensors, its elements. |
| `func (t *Tensor) Sub(other *Tensor) (*Tensor, error)` | Sub returns the elementwise difference t-other. |
| `func (t *Tensor) Sum() float64` | Sum returns the sum of all elements of t. |
| `func (t *Tensor) SumAxis(axis int) (*Tensor, error)` | SumAxis returns a tensor with the given axis summed out, reducing the rank by one. |
| `func (t *Tensor) SwapAxes(a, b int) (*Tensor, error)` | SwapAxes returns a new tensor identical to t but with axes a and b exchanged. |
| `func (t *Tensor) ToMatrix() ([][]float64, error)` | ToMatrix returns the elements of a rank-2 tensor as a freshly allocated row-major [][]float64. |
| `func (t *Tensor) ToVector() ([]float64, error)` | ToVector returns a copy of a rank-1 tensor's elements. |
| `func (t *Tensor) Trace() (float64, error)` | Trace returns the sum of the diagonal of a rank-2 tensor. |
| `func (t *Tensor) Transpose() *Tensor` | Transpose returns a new tensor with all axes reversed. |

</details>

### Variables

`ErrShape`, `ErrIndex`, `ErrAxis`, `ErrRank`, `ErrDataLength`, `ErrSpec`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/tensor).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
