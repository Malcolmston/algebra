# matrix — provides exact symbolic and numeric linear algebra built on top

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/matrix.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/matrix)

Package matrix provides exact symbolic and numeric linear algebra built on top
of the github.com/malcolmston/algebra expression engine.

Every entry of a `Matrix` (and every component of a `Vector`) is an
algebra.Expr, so a matrix can hold arbitrary-precision integers, exact
rationals, symbols and full symbolic expressions side by side. All arithmetic
keeps results exact by canonicalizing every entry through algebra.Simplify,
which means operations such as `Matrix.Det`, `Matrix.Inverse`, `Solve` and
`Matrix.CharPoly` return closed-form, deterministic answers rather than
floating-point approximations.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/matrix"
```

## Usage

This is the package's own `ExampleMatrix_Eigenvalues`, so it compiles and its output is
asserted on every `go test ./matrix/`.

```go
m := matrix.FromInts([][]int64{{2, 0}, {0, 3}})
	ev, _ := m.Eigenvalues()
	for _, e := range ev {
		fmt.Println(e)
	}
```

```
2
3
```

## Exported surface

### Types

| Type | What it is |
| --- | --- |
| `Dense` | Dense is a performance-focused, flat row-major matrix of float64 values. |
| `Matrix` | Matrix is a two-dimensional, rows×cols array of algebra.Expr entries. |
| `Vector` | Vector is a one-dimensional slice of algebra.Expr components. |

<details>
<summary><code>Dense</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewDense(rows, cols int) *Dense` | NewDense returns a rows×cols dense matrix with every element initialized to zero. |
| `func (d *Dense) AddScaled(o *Dense, s float64) (*Dense, error)` | AddScaled performs the in-place AXPY update d += s*o and returns the receiver. |
| `func (d *Dense) At(i, j int) float64` | At returns the element at row i, column j (both 0-based). |
| `func (d *Dense) Cols() int` | Cols returns the number of columns. |
| `func (d *Dense) MulInto(a, b *Dense) (*Dense, error)` | MulInto computes the matrix product a·b and writes it into the receiver, returning the receiver for convenience. |
| `func (d *Dense) Rows() int` | Rows returns the number of rows. |
| `func (d *Dense) Scale(s float64) *Dense` | Scale multiplies every element of the receiver by s in place and returns the receiver. |
| `func (d *Dense) Set(i, j int, v float64)` | Set stores v at row i, column j (both 0-based). |
| `func (d *Dense) ToMatrix() *Matrix` | ToMatrix converts the dense matrix back into a symbolic `Matrix` whose entries are inexact algebra.Flt literals, via `FromFloats`. |
| `func (d *Dense) TransposeInto(dst *Dense) *Dense` | TransposeInto writes the transpose of the receiver into dst and returns dst. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Diag(diag ...algebra.Expr) *Matrix` | Diag returns the square diagonal matrix whose main diagonal holds the given expressions and whose off-diagonal entries are 0. |
| `func FromExpr(vals [][]algebra.Expr) *Matrix` | FromExpr builds a matrix from a rectangular slice of expressions. |
| `func FromFloats(vals [][]float64) *Matrix` | FromFloats builds a matrix from a rectangular slice of float64 values, storing each as an inexact algebra.Flt literal. |
| `func FromInts(vals [][]int64) *Matrix` | FromInts builds a matrix from a rectangular slice of int64 values. |
| `func Identity(n int) *Matrix` | Identity returns the n×n identity matrix. |
| `func New(rows, cols int) *Matrix` | New returns a rows×cols matrix with every entry initialized to 0. |
| `func Ones(rows, cols int) *Matrix` | Ones returns a rows×cols matrix with every entry equal to 1. |
| `func Zeros(rows, cols int) *Matrix` | Zeros returns a rows×cols matrix of zeros. |
| `func (m *Matrix) Add(n *Matrix) (*Matrix, error)` | Add returns the entrywise sum m+n, each entry simplified. |
| `func (m *Matrix) Adjugate() (*Matrix, error)` | Adjugate returns the adjugate (classical adjoint): the transpose of the cofactor matrix. |
| `func (m *Matrix) At(i, j int) algebra.Expr` | At returns the entry at row i, column j (both 0-based). |
| `func (m *Matrix) CharPoly() (algebra.Expr, error)` | CharPoly returns the characteristic polynomial det(A - λI) as an algebra.Expr in the symbol `Lambda`. |
| `func (m *Matrix) Cholesky() (*Matrix, error)` | Cholesky returns the lower-triangular factor L of the (numeric) symmetric positive-definite matrix m, such that A = L·Lᵀ. |
| `func (m *Matrix) Clone() *Matrix` | Clone returns a deep copy of the matrix. |
| `func (m *Matrix) Cofactor(i, j int) (algebra.Expr, error)` | Cofactor returns the (i,j) cofactor: (-1)^(i+j) times the (i,j) minor. |
| `func (m *Matrix) Col(j int) *Vector` | Col returns a copy of column j as a `Vector`. |
| `func (m *Matrix) Cols() int` | Cols returns the number of columns. |
| `func (m *Matrix) ColumnSpaceExact() ([]*Vector, error)` | ColumnSpaceExact returns an exact basis of the column space (image) of A. |
| `func (m *Matrix) Cond2() (float64, error)` | Cond2 returns the spectral (2-norm) condition number of a numeric matrix, the ratio s[0]/s`k-1` of its largest to smallest singular value with k =… |
| `func (m *Matrix) CondP(p int) (float64, error)` | CondP returns the induced condition number of m for the norm selected by p, defined as ‖A‖ₚ · ‖A⁻¹‖ₚ. |
| `func (m *Matrix) DeleteRowCol(i, j int) (*Matrix, error)` | DeleteRowCol returns the matrix with row i and column j removed, the equivalent of SymPy's Matrix.minor_submatrix(). |
| `func (m *Matrix) Det() (algebra.Expr, error)` | Det returns the exact determinant computed by cofactor (Laplace) expansion, which performs no division and therefore stays exact for symbolic entries. |
| `func (m *Matrix) DetLU() (algebra.Expr, error)` | DetLU returns the determinant of the (numeric) square matrix m as an algebra.Flt, computed as sign·∏diag(U) from the partial-pivot LU… |
| `func (m *Matrix) Diagonal() *Vector` | Diagonal returns the main diagonal of m as a `Vector`. |
| `func (m *Matrix) EigSym() (values []float64, vectors *Matrix, err error)` | EigSym computes the eigenvalues and eigenvectors of a real symmetric matrix using the cyclic Jacobi eigenvalue algorithm. |
| `func (m *Matrix) EigSymValues() ([]float64, error)` | EigSymValues computes only the eigenvalues of a real symmetric matrix, in ascending order, using the cyclic Jacobi eigenvalue algorithm. |
| `func (m *Matrix) Eigenvalues() ([]algebra.Expr, error)` | Eigenvalues returns the eigenvalues of a square matrix. |
| `func (m *Matrix) EigenvaluesNumeric() ([]complex128, error)` | EigenvaluesNumeric computes all n eigenvalues of a general real square matrix. |
| `func (m *Matrix) Equal(n *Matrix) bool` | Equal reports whether m and n have the same shape and structurally equal entries. |
| `func (m *Matrix) Exp() (*Matrix, error)` | Exp computes the matrix exponential e^A of a real square matrix by the scaling-and-squaring method with a diagonal Padé approximant. |
| `func (m *Matrix) ExpScaled(t float64) (*Matrix, error)` | ExpScaled returns e^{tA} for the scalar t, the flow map of the linear ODE x' = Ax evaluated at time t. |
| `func (m *Matrix) Floats() ([][]float64, error)` | Floats converts every entry to a float64, returning `ErrUnsupported` wrapping the underlying evaluation error if any entry contains a free symbol or… |
| `func (m *Matrix) HStack(others ...*Matrix) (*Matrix, error)` | HStack returns the horizontal concatenation [m others...], mirroring SymPy's Matrix.hstack. |
| `func (m *Matrix) Hadamard(n *Matrix) (*Matrix, error)` | Hadamard returns the entrywise (Schur) product m∘n, whose (i,j) entry is the product of the (i,j) entries of m and n, each simplified. |
| `func (m *Matrix) Inverse() (*Matrix, error)` | Inverse returns the exact inverse computed as adjugate(m)/det(m). |
| `func (m *Matrix) IsAntiSymmetric() bool` | IsAntiSymmetric reports whether m is square and equal to the negation of its own transpose (m`i``j` == -m`j``i`, which forces a zero diagonal). |
| `func (m *Matrix) IsDiagonal() bool` | IsDiagonal reports whether every off-diagonal entry simplifies to 0. |
| `func (m *Matrix) IsIdentity() bool` | IsIdentity reports whether m is square, has 1 on every main-diagonal entry and 0 everywhere else. |
| `func (m *Matrix) IsLowerTriangular() bool` | IsLowerTriangular reports whether every entry strictly above the main diagonal simplifies to 0. |
| `func (m *Matrix) IsSquare() bool` | IsSquare reports whether the matrix has equal row and column counts. |
| `func (m *Matrix) IsSymmetric() bool` | IsSymmetric reports whether m is square and equal to its own transpose. |
| `func (m *Matrix) IsUpperTriangular() bool` | IsUpperTriangular reports whether every entry strictly below the main diagonal simplifies to 0. |
| `func (m *Matrix) IsZero() bool` | IsZero reports whether every entry of m simplifies to 0. |
| `func (m *Matrix) Kron(n *Matrix) *Matrix` | Kron returns the Kronecker product m⊗n, a block matrix of shape (m.Rows()·n.Rows())×(m.Cols()·n.Cols()) whose (i,j) block is m[i,j]·n. |
| `func (m *Matrix) LU() (l, u, p *Matrix, sign int, err error)` | LU computes the partial-pivot LU factorization P*A = L*U of the (numeric) square matrix m. |
| `func (m *Matrix) Map(f func(algebra.Expr) algebra.Expr) *Matrix` | Map returns a new matrix whose (i,j) entry is the simplified result of f(m`i``j`), mirroring SymPy's Matrix.applyfunc. |
| `func (m *Matrix) Minor(i, j int) (algebra.Expr, error)` | Minor returns the (i,j) minor: the determinant of the submatrix obtained by deleting row i and column j. |
| `func (m *Matrix) Mul(n *Matrix) (*Matrix, error)` | Mul returns the matrix product m·n, each entry simplified. |
| `func (m *Matrix) MulVec(v *Vector) (*Vector, error)` | MulVec returns the matrix-vector product m·v as a `Vector`, treating v as a column vector. |
| `func (m *Matrix) Neg() *Matrix` | Neg returns -m (every entry negated). |
| `func (m *Matrix) Norm1() (float64, error)` | Norm1 returns the induced 1-norm of m, the maximum absolute column sum. |
| `func (m *Matrix) NormFro() algebra.Expr` | NormFro returns the exact Frobenius norm of m, the square root of the sum of the squares of every entry. |
| `func (m *Matrix) NormInf() (float64, error)` | NormInf returns the induced infinity-norm of m, the maximum absolute row sum. |
| `func (m *Matrix) NormMax() (float64, error)` | NormMax returns the max-norm of m, the largest absolute value among all entries. |
| `func (m *Matrix) Nullspace(tol float64) (*Matrix, error)` | Nullspace returns an orthonormal basis of the kernel of a numeric matrix as the columns of the returned matrix. |
| `func (m *Matrix) NullspaceExact() ([]*Vector, error)` | NullspaceExact returns an exact basis of the kernel ker(A) = {x : A·x = 0}. |
| `func (m *Matrix) Pinv(tol float64) (*Matrix, error)` | Pinv returns the Moore-Penrose pseudoinverse V·diag(1/s)·Uᵀ of a numeric matrix, an n×m matrix. |
| `func (m *Matrix) Pow(p int) (*Matrix, error)` | Pow returns m raised to the non-negative integer power p. |
| `func (m *Matrix) QR() (q, r *Matrix, err error)` | QR computes the Householder QR factorization A = Q·R of the (numeric) matrix m, which must have at least as many rows as columns (m ≥ n). |
| `func (m *Matrix) RREF() (*Matrix, int)` | RREF returns the reduced row-echelon form of m together with the rank (the number of nonzero pivot rows). |
| `func (m *Matrix) Rank() int` | Rank returns the rank of m, the number of linearly independent rows, determined from its reduced row-echelon form. |
| `func (m *Matrix) RankNumeric(tol float64) (int, error)` | RankNumeric returns the numerical rank of a numeric matrix: the number of singular values strictly greater than tol. |
| `func (m *Matrix) Row(i int) *Vector` | Row returns a copy of row i as a `Vector`. |
| `func (m *Matrix) RowSpaceExact() ([]*Vector, error)` | RowSpaceExact returns an exact basis of the row space of A. |
| `func (m *Matrix) Rows() int` | Rows returns the number of rows. |
| `func (m *Matrix) SVD() (u *Matrix, s []float64, v *Matrix, err error)` | SVD returns the thin (economy) singular value decomposition A = U·diag(s)·Vᵀ of a numeric matrix. |
| `func (m *Matrix) ScalarMul(s algebra.Expr) *Matrix` | ScalarMul returns m with every entry multiplied by the expression s. |
| `func (m *Matrix) Set(i, j int, e algebra.Expr)` | Set stores the simplified value of e at row i, column j (both 0-based). |
| `func (m *Matrix) SingularValues() ([]float64, error)` | SingularValues returns the singular values of a numeric matrix in non-increasing order. |
| `func (m *Matrix) SolveLU(b *Vector) (*Vector, error)` | SolveLU solves the square linear system m·x = b for the (numeric) matrix m using its partial-pivot LU factorization followed by forward and back… |
| `func (m *Matrix) String() string` | String renders the matrix as an aligned, bracketed grid. |
| `func (m *Matrix) Sub(n *Matrix) (*Matrix, error)` | Sub returns the entrywise difference m-n, each entry simplified. |
| `func (m *Matrix) SubMatrix(r0, r1, c0, c1 int) (*Matrix, error)` | SubMatrix returns the rectangular block of m spanning rows [r0, r1) and columns [c0, c1), the equivalent of SymPy's A[r0:r1, c0:c1]. |
| `func (m *Matrix) Subs(sym algebra.Expr, val algebra.Expr) *Matrix` | Subs substitutes the expression val for the symbol sym in every entry, simplifying the result. |
| `func (m *Matrix) SwapCols(i, j int) (*Matrix, error)` | SwapCols returns a copy of m with columns i and j exchanged. |
| `func (m *Matrix) SwapRows(i, j int) (*Matrix, error)` | SwapRows returns a copy of m with rows i and j exchanged. |
| `func (m *Matrix) ToDense() (*Dense, error)` | ToDense converts the symbolic matrix to a numeric `Dense` using the numeric fast-path of `Matrix.Floats`. |
| `func (m *Matrix) Trace() (algebra.Expr, error)` | Trace returns the sum of the main-diagonal entries. |
| `func (m *Matrix) Transpose() *Matrix` | Transpose returns the transpose of m. |
| `func (m *Matrix) VStack(others ...*Matrix) (*Matrix, error)` | VStack returns the vertical concatenation of m above others, mirroring SymPy's Matrix.vstack. |

</details>

<details>
<summary><code>Vector</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LeastSquares(a *Matrix, b *Vector) (*Vector, error)` | LeastSquares solves the m×n linear system A·x = b in the least-squares sense and returns the solution vector x. |
| `func NewVector(components ...algebra.Expr) *Vector` | NewVector returns a vector holding the simplified given components. |
| `func Solve(a *Matrix, b *Vector) (*Vector, error)` | Solve returns the unique solution x of the linear system a·x = b, where a is an m×n coefficient matrix and b is a length-m vector. |
| `func VectorFromInts(vals ...int64) *Vector` | VectorFromInts returns a vector built from int64 components. |
| `func (v *Vector) Add(w *Vector) (*Vector, error)` | Add returns the component-wise sum v+w. |
| `func (v *Vector) Angle(w *Vector) (float64, error)` | Angle returns the angle between v and w in radians as a numeric float64, computed as arccos(v·w / (‖v‖‖w‖)). |
| `func (v *Vector) At(i int) algebra.Expr` | At returns the i-th component (0-based). |
| `func (v *Vector) Clone() *Vector` | Clone returns a copy of the vector. |
| `func (v *Vector) ColMatrix() *Matrix` | ColMatrix returns the vector as an n×1 column matrix. |
| `func (v *Vector) Cross(w *Vector) (*Vector, error)` | Cross returns the 3-dimensional cross product v×w. |
| `func (v *Vector) Dot(w *Vector) (algebra.Expr, error)` | Dot returns the dot product v·w as a single simplified expression. |
| `func (v *Vector) Equal(w *Vector) bool` | Equal reports whether v and w have the same length and structurally equal components (compared after simplification). |
| `func (v *Vector) Len() int` | Len returns the number of components. |
| `func (v *Vector) Neg() *Vector` | Neg returns -v. |
| `func (v *Vector) Norm() algebra.Expr` | Norm returns the Euclidean norm sqrt(v·v) as an exact symbolic expression. |
| `func (v *Vector) Normalize() (*Vector, error)` | Normalize returns the unit vector v/‖v‖ with each component simplified. |
| `func (v *Vector) RowMatrix() *Matrix` | RowMatrix returns the vector as a 1×n row matrix. |
| `func (v *Vector) ScalarMul(s algebra.Expr) *Vector` | ScalarMul returns the vector scaled by the expression s. |
| `func (v *Vector) Set(i int, e algebra.Expr)` | Set stores the simplified value e at index i (0-based). |
| `func (v *Vector) String() string` | String renders the vector as a bracketed, comma-separated list. |
| `func (v *Vector) Sub(w *Vector) (*Vector, error)` | Sub returns the component-wise difference v-w. |

</details>

### Constants

`NormInfinity`

### Variables

`ErrDimension`, `ErrNotSquare`, `ErrSingular`, `ErrInconsistent`, `ErrUnderdetermined`, `ErrUnsupported`, `ErrNotPositiveDefinite`, `Lambda`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/matrix).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
