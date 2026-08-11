# lattices — implements algorithms for integer and real lattices

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/lattices.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/lattices)

Package lattices implements algorithms for integer and real lattices.

- Gram matrices, Gram determinants and lattice covolume (determinant). -
Gram-Schmidt orthogonalization with the mu coefficient matrix. - LLL
(Lenstra-Lenstra-Lovasz) basis reduction and standalone size reduction. -
Shortest-vector (SVP) and closest-vector (CVP) search by Fincke-Pohst
enumeration. - Babai's rounding and nearest-plane approximate CVP solvers. -
Minkowski and Hermite bounds and successive-minima estimates. - The dual
(reciprocal) lattice. - Hermite normal form of integer matrices.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/lattices"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BallVolume(n int, r float64) float64` | BallVolume returns the volume of the n-dimensional Euclidean ball of radius r, namely pi^(n/2)/Gamma(n/2+1) * r^n. |
| `func GaussianHeuristic(n int, det float64) float64` | GaussianHeuristic returns the Gaussian heuristic estimate of the first minimum, (det / V_n)^(1/n), the expected shortest length for a random lattice… |
| `func HermiteBound(n int, det float64) float64` | HermiteBound returns the Hermite bound sqrt(gamma_n) * det^(1/n) on the first minimum lambda_1 of a rank-n lattice with the given determinant. |
| `func HermiteConstant(n int) float64` | HermiteConstant returns Hermite's constant gamma_n. |
| `func HermiteConstantLowerBound(n int) float64` | HermiteConstantLowerBound returns a Minkowski-Hlawka style lower bound on Hermite's constant, gamma_n >= (n/(2*pi*e)) for large n; it is only a rough… |
| `func HermiteConstantUpperBound(n int) float64` | HermiteConstantUpperBound returns Minkowski's upper bound on Hermite's constant, gamma_n <= (2/pi) * Gamma(2 + n/2)^(2/n). |
| `func MinkowskiBound(n int, det float64) float64` | MinkowskiBound returns the upper bound on the first minimum lambda_1 given by Minkowski's convex body theorem for a lattice of the given rank and… |
| `func MinkowskiSecondLower(n int, det float64) float64` | MinkowskiSecondLower returns the lower bound of Minkowski's second theorem on the product of the successive minima: prod_i lambda_i >= 2^n * det /… |
| `func MinkowskiSecondUpper(n int, det float64) float64` | MinkowskiSecondUpper returns the upper bound of Minkowski's second theorem on the product of the successive minima: prod_i lambda_i <= 2^n * det /… |
| `func RankinBound(n, m int) float64` | RankinBound returns Rankin's bound gamma_{n,m} <= gamma_n^m used to bound the minimum determinant of an m-dimensional sublattice; here it returns… |
| `func SameLattice(a, b Basis) (bool, error)` | SameLattice reports whether two integer bases of the same ambient dimension generate the same lattice, by comparing their Hermite normal forms. |
| `func UnitBallVolume(n int) float64` | UnitBallVolume returns the volume of the n-dimensional unit ball, pi^(n/2)/Gamma(n/2+1). |

### Types

| Type | What it is |
| --- | --- |
| `Basis` | Basis is a lattice basis: an ordered list of row vectors that generate the lattice by integer linear combination. |
| `EnumResult` | EnumResult holds the outcome of a lattice enumeration search: the lattice Vector found, its integer coordinate vector Coeffs in the given basis, and… |
| `GramSchmidt` | GramSchmidt holds the result of Gram-Schmidt orthogonalization of a basis: the orthogonal vectors Star (the b*_i), the strictly lower-triangular… |
| `IntMatrix` | IntMatrix is a dense matrix of exact integers (*big.Int) stored row-major. |
| `Matrix` | Matrix is a dense real matrix stored row-major with float64 entries. |
| `RatMatrix` | RatMatrix is a dense matrix of exact rationals stored row-major. |
| `RatVec` | RatVec is an exact rational vector represented as a slice of *big.Rat coordinates. |
| `Vec` | Vec is a real vector represented as a slice of float64 coordinates. |

<details>
<summary><code>Basis</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BasisFromMatrix(m Matrix) Basis` | BasisFromMatrix builds a Basis whose vectors are the rows of m. |
| `func BasisFromRows(rows [][]float64) Basis` | BasisFromRows builds a Basis from a slice of float64 row slices. |
| `func IdentityBasis(n int) Basis` | IdentityBasis returns the standard basis of Z^n (the n rows of the identity). |
| `func NewBasis(rows ...Vec) Basis` | NewBasis returns a Basis containing copies of the given row vectors. |
| `func (b Basis) BabaiCoeffs(target Vec) []int64` | BabaiCoeffs returns the integer coordinate vector produced by Babai rounding of target, so that b.Point(result) equals BabaiRound(target). |
| `func (b Basis) BabaiNearestPlane(target Vec) Vec` | BabaiNearestPlane returns an approximate closest lattice vector to target using Babai's nearest-plane algorithm, which walks the Gram-Schmidt frame… |
| `func (b Basis) BabaiRound(target Vec) Vec` | BabaiRound returns an approximate closest lattice vector to target using Babai's rounding technique: express the target in basis coordinates (in a… |
| `func (b Basis) Biorthogonality() (Matrix, error)` | Biorthogonality returns the matrix M`i``j` = <d_i, b_j> where d is the dual basis; it should equal the identity matrix up to numerical error. |
| `func (b Basis) Clone() Basis` | Clone returns an independent deep copy of the basis. |
| `func (b Basis) ClosestVector(target Vec) (EnumResult, error)` | ClosestVector returns the lattice vector nearest to target using enumeration seeded by Babai's nearest-plane approximation for the search radius. |
| `func (b Basis) Coordinates(v Vec) (Vec, error)` | Coordinates returns the coefficient vector x such that b represents v, that is sum_i x`i`*b_i = v, by solving the (possibly overdetermined) system in… |
| `func (b Basis) Determinant() float64` | Determinant returns the lattice determinant (covolume) sqrt(det G). |
| `func (b Basis) Dim() int` | Dim returns the ambient dimension (length of each vector), or 0 for an empty basis. |
| `func (b Basis) Dual() (Basis, error)` | Dual returns the dual (reciprocal) lattice basis d_0, ..., d_{n-1}, the unique vectors in the span of the basis satisfying <d_i, b_j> = 1 if i == j… |
| `func (b Basis) DualDeterminant() float64` | DualDeterminant returns the covolume of the dual lattice, which equals the reciprocal of the covolume of the primal lattice. |
| `func (b Basis) DualGramDeterminantRat() (*big.Rat, error)` | DualGramDeterminantRat returns the exact determinant of the dual Gram matrix, which equals the reciprocal of the primal Gram determinant. |
| `func (b Basis) DualRat() ([]RatVec, error)` | DualRat returns the exact rational dual basis, computed with exact rational arithmetic. |
| `func (b Basis) EnumerateVectors(radius float64, limit int) ([]EnumResult, bool, error)` | EnumerateVectors returns all nonzero lattice vectors whose Euclidean norm is at most radius, each as an EnumResult, together with a completeness flag… |
| `func (b Basis) FirstMinimum() (float64, error)` | FirstMinimum returns lambda_1, the length of the shortest nonzero lattice vector, computed exactly by enumeration. |
| `func (b Basis) GaussianHeuristicBasis() float64` | GaussianHeuristicBasis returns the Gaussian heuristic evaluated for the basis. |
| `func (b Basis) Gram() Matrix` | Gram returns the Gram matrix G with G`i``j` = <b_i, b_j>. |
| `func (b Basis) GramDeterminant() float64` | GramDeterminant returns det(G), the determinant of the Gram matrix. |
| `func (b Basis) GramDeterminantRat() *big.Rat` | GramDeterminantRat returns the exact determinant of the Gram matrix. |
| `func (b Basis) GramInverseRat() (RatMatrix, error)` | GramInverseRat returns the exact inverse of the rational Gram matrix, the Gram matrix of the dual basis. |
| `func (b Basis) GramRat() RatMatrix` | GramRat returns the exact rational Gram matrix of the basis. |
| `func (b Basis) GramSchmidtNorms() []float64` | GramSchmidtNorms returns the squared norms of the Gram-Schmidt vectors b*_i. |
| `func (b Basis) HadamardRatio() float64` | HadamardRatio returns covolume / prod_i \|b_i\|, a value in (0, 1]. |
| `func (b Basis) HermiteBoundBasis() float64` | HermiteBoundBasis returns HermiteBound evaluated at the rank and determinant of the basis. |
| `func (b Basis) HermiteNormalForm() (IntMatrix, error)` | HermiteNormalForm returns the row Hermite normal form of an integer basis. |
| `func (b Basis) IsDualTo(other Basis, tol float64) bool` | IsDualTo reports whether other is (numerically, within tol) the dual basis of b, i.e. |
| `func (b Basis) IsFullRank(tol float64) bool` | IsFullRank reports whether the basis vectors are linearly independent, i.e. |
| `func (b Basis) IsLLLReduced(delta, tol float64) bool` | IsLLLReduced reports whether the basis is LLL reduced for parameter delta within tolerance tol: it must be size reduced and satisfy the Lovasz… |
| `func (b Basis) IsSizeReduced(tol float64) bool` | IsSizeReduced reports whether every Gram-Schmidt coefficient Mu`i``j` with j < i has absolute value at most 1/2 + tol, the size-reduced condition. |
| `func (b Basis) IsSquare() bool` | IsSquare reports whether the number of vectors equals the ambient dimension. |
| `func (b Basis) KissingNumber() (int, error)` | KissingNumber returns the number of shortest nonzero lattice vectors, that is the number of lattice vectors of length exactly lambda_1 (which always… |
| `func (b Basis) LLL(delta float64) (Basis, error)` | LLL returns an LLL-reduced basis for the same lattice using the given Lovasz parameter delta, which must lie in (1/4, 1]. |
| `func (b Basis) LLLDefault() Basis` | LLLDefault returns an LLL-reduced basis using DefaultDelta (3/4). |
| `func (b Basis) LLLGramDefect() float64` | LLLGramDefect returns the orthogonality defect of the LLL-reduced basis, a convenient quality measure after reduction. |
| `func (b Basis) Matrix() Matrix` | Matrix returns the basis as a Matrix whose rows are the basis vectors. |
| `func (b Basis) MaxNorm() float64` | MaxNorm returns the largest Euclidean norm among the basis vectors. |
| `func (b Basis) MinkowskiBoundBasis() float64` | MinkowskiBoundBasis returns MinkowskiBound evaluated at the rank and determinant of the basis. |
| `func (b Basis) MuMatrix() Matrix` | MuMatrix returns the strictly lower-triangular matrix of Gram-Schmidt coefficients (with ones on the diagonal) as a Matrix. |
| `func (b Basis) OrthogonalNormal() Basis` | OrthogonalNormal returns the Gram-Schmidt vectors normalized to unit length, forming an orthonormal basis of the span. |
| `func (b Basis) OrthogonalityDefect() float64` | OrthogonalityDefect returns the ratio prod_i \|b_i\| / covolume, a measure of how far the basis is from orthogonal. |
| `func (b Basis) Orthogonalize() GramSchmidt` | Orthogonalize computes the Gram-Schmidt orthogonalization of the basis without normalizing the resulting vectors. |
| `func (b Basis) Point(coeffs []int64) Vec` | Point returns the lattice point sum_i coeffs`i`*b_i for integer coefficients. |
| `func (b Basis) PointFloat(coeffs []float64) Vec` | PointFloat returns sum_i coeffs`i`*b_i for real coefficients (a point in the real span of the lattice). |
| `func (b Basis) ProductOfStarNorms() float64` | ProductOfStarNorms returns prod_i \|b*_i\|, which for a full-rank basis equals the lattice covolume exactly (Gram-Schmidt preserves volume). |
| `func (b Basis) Rank() int` | Rank returns the number of vectors in the basis. |
| `func (b Basis) Row(i int) Vec` | Row returns the i-th basis vector (the stored slice, not a copy). |
| `func (b Basis) ShortestBasisVector() (int, float64)` | ShortestBasisVector returns the index and Euclidean norm of the shortest nonzero vector currently in the basis (not the lattice minimum). |
| `func (b Basis) ShortestVector() (EnumResult, error)` | ShortestVector returns the shortest nonzero vector of the lattice by Fincke-Pohst enumeration. |
| `func (b Basis) SizeReduced() Basis` | SizeReduced returns a copy of the basis that is size reduced: every Gram-Schmidt coefficient mu`i``j` with j < i is made to satisfy \|mu`i``j`\| <=… |
| `func (b Basis) StarBasis() Basis` | StarBasis returns the orthogonalized vectors b*_i as a Basis. |
| `func (b Basis) SuccessiveMinima() ([]float64, error)` | SuccessiveMinima returns the successive minima lambda_1 <= ... |
| `func (b Basis) Swap(i, j int)` | Swap exchanges basis vectors i and j in place. |
| `func (b Basis) Volume() float64` | Volume is an alias for Determinant: the covolume of the lattice. |

</details>

<details>
<summary><code>EnumResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (r EnumResult) Norm() float64` | Norm returns the Euclidean norm (or distance) sqrt(Norm2) of the result. |

</details>

<details>
<summary><code>IntMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IdentityIntMatrix(n int) IntMatrix` | IdentityIntMatrix returns the n-by-n integer identity matrix. |
| `func IntMatrixFromBasis(b Basis) (IntMatrix, error)` | IntMatrixFromBasis converts a basis with integer entries to an IntMatrix. |
| `func NewIntMatrix(rows [][]int64) IntMatrix` | NewIntMatrix builds an integer matrix from int64 rows. |
| `func ZeroIntMatrix(r, c int) IntMatrix` | ZeroIntMatrix returns an r-by-c integer matrix of zeros. |
| `func (m IntMatrix) At(i, j int) *big.Int` | At returns a copy of the entry in row i and column j. |
| `func (m IntMatrix) Clone() IntMatrix` | Clone returns an independent deep copy of m. |
| `func (m IntMatrix) Cols() int` | Cols returns the number of columns of m. |
| `func (m IntMatrix) Det() (*big.Int, error)` | Det returns the exact integer determinant of a square matrix using the fraction-free Bareiss algorithm. |
| `func (m IntMatrix) Equal(n IntMatrix) bool` | Equal reports whether m and n have the same shape and equal entries. |
| `func (m IntMatrix) HNFRank() int` | HNFRank returns the number of nonzero pivot rows in the Hermite normal form, which equals the rank of the integer matrix. |
| `func (m IntMatrix) HermiteNormalForm() IntMatrix` | HermiteNormalForm returns the row Hermite normal form of the integer matrix: an upper-triangular row-echelon matrix H obtained from m by unimodular… |
| `func (m IntMatrix) HermiteNormalFormWithTransform() (h IntMatrix, u IntMatrix)` | HermiteNormalFormWithTransform returns the Hermite normal form H together with a unimodular matrix U such that U*m = H. |
| `func (m IntMatrix) IsHermiteNormalForm() bool` | IsHermiteNormalForm reports whether m is already in row Hermite normal form: pivots occur in strictly increasing columns, every pivot is positive,… |
| `func (m IntMatrix) IsSquare() bool` | IsSquare reports whether m has equal numbers of rows and columns. |
| `func (m IntMatrix) Mul(n IntMatrix) IntMatrix` | Mul returns the exact product m*n. |
| `func (m IntMatrix) Rat() RatMatrix` | Rat returns the exact rational version of m. |
| `func (m IntMatrix) Rows() int` | Rows returns the number of rows of m. |
| `func (m IntMatrix) Set(i, j int, x *big.Int)` | Set stores a copy of x in row i and column j. |
| `func (m IntMatrix) SetInt64(i, j int, x int64)` | SetInt64 stores the int64 value x in row i and column j. |
| `func (m IntMatrix) String() string` | String renders m with one row per line. |
| `func (m IntMatrix) Transpose() IntMatrix` | Transpose returns the transpose of m. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DiagMatrix(diag ...float64) Matrix` | DiagMatrix returns the square matrix with the given diagonal entries. |
| `func IdentityMatrix(n int) Matrix` | IdentityMatrix returns the n-by-n identity matrix. |
| `func NewMatrix(rows [][]float64) Matrix` | NewMatrix builds a rows-by-cols matrix from the given row slices. |
| `func ZeroMatrix(r, c int) Matrix` | ZeroMatrix returns an r-by-c matrix of zeros. |
| `func (m Matrix) Add(n Matrix) Matrix` | Add returns m+n. |
| `func (m Matrix) ApproxEqual(n Matrix, tol float64) bool` | ApproxEqual reports whether m and n have the same shape and every entry agrees within the absolute tolerance tol. |
| `func (m Matrix) At(i, j int) float64` | At returns the entry in row i and column j. |
| `func (m Matrix) Clone() Matrix` | Clone returns an independent copy of m. |
| `func (m Matrix) Col(j int) Vec` | Col returns a copy of column j as a Vec. |
| `func (m Matrix) Cols() int` | Cols returns the number of columns of m. |
| `func (m Matrix) Det() (float64, error)` | Det returns the determinant of a square matrix via LU decomposition with partial pivoting. |
| `func (m Matrix) Equal(n Matrix) bool` | Equal reports whether m and n have the same shape and entries. |
| `func (m Matrix) Inverse() (Matrix, error)` | Inverse returns the inverse of a square matrix by Gauss-Jordan elimination with partial pivoting. |
| `func (m Matrix) IsSquare() bool` | IsSquare reports whether m has equal numbers of rows and columns. |
| `func (m Matrix) IsSymmetric(tol float64) bool` | IsSymmetric reports whether m is square and equal to its transpose within the tolerance tol. |
| `func (m Matrix) Mul(n Matrix) Matrix` | Mul returns the matrix product m*n. |
| `func (m Matrix) MulVec(v Vec) Vec` | MulVec returns the matrix-vector product m*v. |
| `func (m Matrix) Rank(tol float64) int` | Rank returns the numerical rank of m computed by Gaussian elimination, using tol as the threshold below which a pivot is treated as zero. |
| `func (m Matrix) Row(i int) Vec` | Row returns a copy of row i as a Vec. |
| `func (m Matrix) Rows() int` | Rows returns the number of rows of m. |
| `func (m Matrix) Scale(s float64) Matrix` | Scale returns the scalar multiple s*m. |
| `func (m Matrix) Set(i, j int, x float64)` | Set stores x in row i and column j. |
| `func (m Matrix) Solve(b Vec) (Vec, error)` | Solve returns the solution x of m*x = b for a square matrix m using Gaussian elimination with partial pivoting. |
| `func (m Matrix) String() string` | String renders m with one row per line. |
| `func (m Matrix) Sub(n Matrix) Matrix` | Sub returns m-n. |
| `func (m Matrix) Trace() (float64, error)` | Trace returns the sum of the diagonal entries of a square matrix. |
| `func (m Matrix) Transpose() Matrix` | Transpose returns the transpose of m. |

</details>

<details>
<summary><code>RatMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IdentityRatMatrix(n int) RatMatrix` | IdentityRatMatrix returns the n-by-n rational identity matrix. |
| `func NewRatMatrix(rows [][]int64) RatMatrix` | NewRatMatrix builds a rows-by-cols rational matrix from integer rows. |
| `func RatMatrixFromMatrix(m Matrix) RatMatrix` | RatMatrixFromMatrix converts a float64 Matrix to an exact RatMatrix. |
| `func ZeroRatMatrix(r, c int) RatMatrix` | ZeroRatMatrix returns an r-by-c rational matrix of zeros. |
| `func (m RatMatrix) Add(n RatMatrix) RatMatrix` | Add returns the exact sum m+n. |
| `func (m RatMatrix) At(i, j int) *big.Rat` | At returns a copy of the entry in row i and column j. |
| `func (m RatMatrix) Clone() RatMatrix` | Clone returns an independent deep copy of m. |
| `func (m RatMatrix) Cols() int` | Cols returns the number of columns of m. |
| `func (m RatMatrix) Det() (*big.Rat, error)` | Det returns the exact determinant of a square rational matrix by fraction (exact) Gaussian elimination. |
| `func (m RatMatrix) Equal(n RatMatrix) bool` | Equal reports whether m and n have the same shape and exactly equal entries. |
| `func (m RatMatrix) Float() Matrix` | Float returns the float64 approximation of m as a Matrix. |
| `func (m RatMatrix) Inverse() (RatMatrix, error)` | Inverse returns the exact inverse of a square rational matrix by Gauss-Jordan elimination. |
| `func (m RatMatrix) IsSquare() bool` | IsSquare reports whether m has equal numbers of rows and columns. |
| `func (m RatMatrix) Mul(n RatMatrix) RatMatrix` | Mul returns the exact product m*n. |
| `func (m RatMatrix) Rows() int` | Rows returns the number of rows of m. |
| `func (m RatMatrix) Set(i, j int, x *big.Rat)` | Set stores a copy of x in row i and column j. |
| `func (m RatMatrix) String() string` | String renders m with one row per line using rational entries. |
| `func (m RatMatrix) Sub(n RatMatrix) RatMatrix` | Sub returns the exact difference m-n. |
| `func (m RatMatrix) Transpose() RatMatrix` | Transpose returns the transpose of m. |

</details>

<details>
<summary><code>RatVec</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewRatVec(n int) RatVec` | NewRatVec returns an n-dimensional rational zero vector with freshly allocated entries. |
| `func RatVecFromFloats(xs ...float64) RatVec` | RatVecFromFloats builds a RatVec from float64 coordinates, converting each exactly to a rational. |
| `func RatVecFromInts(xs ...int64) RatVec` | RatVecFromInts builds a RatVec from integer coordinates. |
| `func (v RatVec) Add(w RatVec) RatVec` | Add returns the exact sum v+w. |
| `func (v RatVec) AddScaled(w RatVec, s *big.Rat) RatVec` | AddScaled returns the exact combination v + s*w. |
| `func (v RatVec) Clone() RatVec` | Clone returns an independent deep copy of v. |
| `func (v RatVec) Dim() int` | Dim returns the number of coordinates in v. |
| `func (v RatVec) Dot(w RatVec) *big.Rat` | Dot returns the exact inner product <v, w>. |
| `func (v RatVec) Equal(w RatVec) bool` | Equal reports whether v and w have identical dimension and coordinates. |
| `func (v RatVec) Float() Vec` | Float returns the float64 approximation of v as a Vec. |
| `func (v RatVec) IsZero() bool` | IsZero reports whether every coordinate of v is exactly zero. |
| `func (v RatVec) Neg() RatVec` | Neg returns -v. |
| `func (v RatVec) Norm2() *big.Rat` | Norm2 returns the exact squared Euclidean norm <v, v>. |
| `func (v RatVec) Scale(s *big.Rat) RatVec` | Scale returns the exact scalar multiple s*v. |
| `func (v RatVec) String() string` | String renders v as a bracketed, space-separated list of rationals. |
| `func (v RatVec) Sub(w RatVec) RatVec` | Sub returns the exact difference v-w. |

</details>

<details>
<summary><code>Vec</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Combine(coeffs []float64, vecs []Vec) Vec` | Combine returns the linear combination sum_i coeffs`i`*vecs`i`. |
| `func IntCombine(coeffs []int64, vecs []Vec) Vec` | IntCombine returns the integer linear combination sum_i coeffs`i`*vecs`i`. |
| `func NewVec(xs ...float64) Vec` | NewVec returns a Vec with the given coordinates. |
| `func UnitVec(n, i int) Vec` | UnitVec returns the n-dimensional standard basis vector e_i (a 1 in position i and zeros elsewhere). |
| `func VecFromInts(xs ...int64) Vec` | VecFromInts builds a Vec from integer coordinates. |
| `func ZeroVec(n int) Vec` | ZeroVec returns the n-dimensional zero vector. |
| `func (v Vec) Add(w Vec) Vec` | Add returns the sum v+w. |
| `func (v Vec) AddScaled(w Vec, s float64) Vec` | AddScaled returns v + s*w. |
| `func (v Vec) Angle(w Vec) float64` | Angle returns the angle in radians between v and w in the range [0, pi]. |
| `func (v Vec) ApproxEqual(w Vec, tol float64) bool` | ApproxEqual reports whether v and w have the same dimension and every coordinate agrees within the absolute tolerance tol. |
| `func (v Vec) Clone() Vec` | Clone returns an independent copy of v. |
| `func (v Vec) Cosine(w Vec) float64` | Cosine returns the cosine of the angle between v and w. |
| `func (v Vec) Dim() int` | Dim returns the number of coordinates in v. |
| `func (v Vec) Dist(w Vec) float64` | Dist returns the Euclidean distance between v and w. |
| `func (v Vec) Dist2(w Vec) float64` | Dist2 returns the squared Euclidean distance between v and w. |
| `func (v Vec) Dot(w Vec) float64` | Dot returns the Euclidean inner product <v, w>. |
| `func (v Vec) Equal(w Vec) bool` | Equal reports whether v and w have identical dimensions and coordinates. |
| `func (v Vec) IsOrthogonal(w Vec, tol float64) bool` | IsOrthogonal reports whether <v, w> is zero within the tolerance tol. |
| `func (v Vec) IsZero() bool` | IsZero reports whether every coordinate of v is exactly zero. |
| `func (v Vec) L1() float64` | L1 returns the L1 (taxicab) norm of v, the sum of absolute values. |
| `func (v Vec) LInf() float64` | LInf returns the L-infinity (maximum absolute value) norm of v. |
| `func (v Vec) Neg() Vec` | Neg returns -v. |
| `func (v Vec) Norm() float64` | Norm returns the Euclidean (L2) norm of v. |
| `func (v Vec) Norm2() float64` | Norm2 returns the squared Euclidean norm <v, v>. |
| `func (v Vec) Normalize() Vec` | Normalize returns the unit vector in the direction of v. |
| `func (v Vec) Proj(w Vec) Vec` | Proj returns the orthogonal projection of v onto the line spanned by w. |
| `func (v Vec) ProjCoeff(w Vec) float64` | ProjCoeff returns the scalar <v, w>/<w, w>, the coefficient of the orthogonal projection of v onto w. |
| `func (v Vec) Reject(w Vec) Vec` | Reject returns the component of v orthogonal to w, that is v minus its projection onto w. |
| `func (v Vec) Round() Vec` | Round returns the vector obtained by rounding each coordinate of v to the nearest integer (halves away from zero). |
| `func (v Vec) Scale(s float64) Vec` | Scale returns the scalar multiple s*v. |
| `func (v Vec) String() string` | String renders v as a bracketed, space-separated list of coordinates. |
| `func (v Vec) Sub(w Vec) Vec` | Sub returns the difference v-w. |
| `func (v Vec) Sum() float64` | Sum returns the sum of the coordinates of v. |

</details>

### Constants

`DefaultDelta`

### Variables

`ErrBadParameter`, `ErrDimMismatch`, `ErrEmpty`, `ErrNoSolution`, `ErrNotFullRank`, `ErrNotSquare`, `ErrRankMismatch`, `ErrSingular`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/lattices).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
