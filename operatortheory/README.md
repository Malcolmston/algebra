# operatortheory — implements finite-dimensional linear operator and spectral theory using only the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/operatortheory.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/operatortheory)

Package operatortheory implements finite-dimensional linear operator and
spectral theory using only the Go standard library.

- Classification predicates: `Matrix.IsHermitian`, `Matrix.IsSkewHermitian`,
`Matrix.IsUnitary`, `Matrix.IsNormal`, `Matrix.IsProjection`,
`Matrix.IsPositiveDefinite`, `Matrix.IsContraction`, `Matrix.IsIsometry` and
many more. - Norms: the Frobenius, induced 1-, 2- (spectral) and infinity
norms, the nuclear (trace) norm and the numerical radius. - Spectral theory:
`Matrix.Eigenvalues`, `Matrix.HermitianEigen`, `Matrix.Spectrum`,
`Matrix.SpectralRadius`, the `Matrix.Resolvent`, the `Matrix.NumericalRange`
and `Matrix.NumericalAbscissa`. - Decompositions: `Matrix.QR`,
`Matrix.Hessenberg`, `Matrix.SVD`, `Matrix.PolarDecomposition` and the
Hermitian spectral decomposition `Matrix.SpectralDecomposition`. - Functional
calculus on normal operators: `Matrix.ApplyFunction` together with the
ready-made `Matrix.Exp`, `Matrix.Log`, `Matrix.Sqrt`, `Matrix.Sign` and
trigonometric functions. - Pseudospectra estimates: `Matrix.ResolventNorm`,
`Matrix.PseudospectralAbscissa` and `Matrix.PseudospectrumGrid`.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/operatortheory"
```

## Exported surface

### Types

| Type | What it is |
| --- | --- |
| `Eigenpair` | Eigenpair bundles an eigenvalue with a corresponding unit eigenvector. |
| `Matrix` | Matrix is a dense complex matrix stored in row-major order. |
| `ScalarFunction` | ScalarFunction is a complex-valued function of a complex variable, the kind of analytic function to which the functional calculus applies. |
| `SpectralComponent` | SpectralComponent bundles a (real) eigenvalue of a Hermitian operator with the orthogonal projector onto its eigenspace. |
| `Vector` | Vector is a finite-dimensional complex vector, an element of the Hilbert space C^n with inner product <x,y> = sum conj(x_i) y_i. |

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Companion(coeffs []complex128) (*Matrix, error)` | Companion returns the companion matrix of the monic polynomial whose coefficients (excluding the leading 1) are given from the constant term up to… |
| `func Diagonal(d []complex128) *Matrix` | Diagonal returns the square diagonal matrix with the given diagonal entries. |
| `func FromReal(rows, cols int, data []float64) (*Matrix, error)` | FromReal returns a complex matrix whose entries are the given real numbers (in row-major order) with zero imaginary part. |
| `func FromRows(rows [][]complex128) (*Matrix, error)` | FromRows builds a matrix from a slice of rows. |
| `func Identity(n int) *Matrix` | Identity returns the n-by-n identity operator. |
| `func JordanBlock(lambda complex128, n int) *Matrix` | JordanBlock returns the n-by-n Jordan block with eigenvalue lambda: lambda on the diagonal and 1 on the first superdiagonal. |
| `func NewMatrix(rows, cols int) *Matrix` | NewMatrix returns a rows-by-cols zero matrix. |
| `func NewMatrixFromData(rows, cols int, data []complex128) (*Matrix, error)` | NewMatrixFromData returns a rows-by-cols matrix filled from data in row-major order. |
| `func OrthogonalProjector(vectors []Vector) (*Matrix, error)` | OrthogonalProjector returns the orthogonal projection operator onto the subspace spanned by the given vectors of common length n. |
| `func RandomHermitian(n int, seed int64) *Matrix` | RandomHermitian returns a random n-by-n Hermitian matrix built as (A + A^H)/2 from a random matrix A seeded by seed. |
| `func RandomMatrix(rows, cols int, seed int64) *Matrix` | RandomMatrix returns an m-by-n matrix whose entries have real and imaginary parts drawn independently from the standard normal distribution, using a… |
| `func RandomUnitary(n int, seed int64) *Matrix` | RandomUnitary returns a random n-by-n unitary matrix obtained from the QR factorisation of a random complex matrix seeded by seed. |
| `func RealDiagonal(d []float64) *Matrix` | RealDiagonal returns the square diagonal matrix with the given real diagonal entries. |
| `func Reflection(v Vector) (*Matrix, error)` | Reflection returns the Householder reflection I - 2 P, where P is the orthogonal projector onto the line spanned by the nonzero unit-normalisable… |
| `func Zero(rows, cols int) *Matrix` | Zero returns the rows-by-cols zero matrix. |
| `func (m *Matrix) Add(b *Matrix) (*Matrix, error)` | Add returns m + b. |
| `func (m *Matrix) Adjoint() *Matrix` | Adjoint returns the conjugate transpose (Hermitian adjoint) m^H. |
| `func (m *Matrix) AntiCommutator(b *Matrix) (*Matrix, error)` | AntiCommutator returns the anticommutator {m,b} = m*b + b*m. |
| `func (m *Matrix) Apply(v Vector) Vector` | Apply returns the image m*v of the vector v under the operator, panicking on a dimension mismatch. |
| `func (m *Matrix) ApplyFunction(f ScalarFunction) (*Matrix, error)` | ApplyFunction returns f(m) defined by the holomorphic functional calculus for a diagonalisable (in particular, normal) operator: if m = V… |
| `func (m *Matrix) At(i, j int) complex128` | At returns the entry in row i and column j. |
| `func (m *Matrix) CayleyTransform() (*Matrix, error)` | CayleyTransform returns (m - i*I)(m + i*I)^{-1}, which maps a Hermitian operator to a unitary one. |
| `func (m *Matrix) CharacteristicPolynomial() ([]complex128, error)` | CharacteristicPolynomial returns the coefficients of the characteristic polynomial det(x*I - m) in ascending order, so the result has length n+1 with… |
| `func (m *Matrix) Clone() *Matrix` | Clone returns an independent copy of the matrix. |
| `func (m *Matrix) Col(j int) Vector` | Col returns a copy of column j as a Vector. |
| `func (m *Matrix) Cols() int` | Cols returns the number of columns. |
| `func (m *Matrix) Commutator(b *Matrix) (*Matrix, error)` | Commutator returns the commutator [m,b] = m*b - b*m. |
| `func (m *Matrix) ConditionNumber() float64` | ConditionNumber returns the 2-norm condition number, the ratio of the largest to the smallest singular value. |
| `func (m *Matrix) Conjugate() *Matrix` | Conjugate returns the entrywise complex conjugate of m. |
| `func (m *Matrix) Cos() (*Matrix, error)` | Cos returns the matrix cosine cos(m) via the functional calculus. |
| `func (m *Matrix) Cosh() (*Matrix, error)` | Cosh returns the matrix hyperbolic cosine cosh(m) via the functional calculus. |
| `func (m *Matrix) Determinant() (complex128, error)` | Determinant returns the determinant of a square matrix. |
| `func (m *Matrix) Diag() Vector` | Diag returns the main diagonal of the matrix as a Vector of length min(rows,cols). |
| `func (m *Matrix) DiagonalizeHermitian() (u, d *Matrix, err error)` | DiagonalizeHermitian returns a unitary U and a real diagonal matrix D such that m = U*D*U^H, for a Hermitian matrix. |
| `func (m *Matrix) Dims() (int, int)` | Dims returns the number of rows and columns. |
| `func (m *Matrix) DirectSum(b *Matrix) *Matrix` | DirectSum returns the block-diagonal matrix diag(m, b). |
| `func (m *Matrix) DistanceFrobenius(b *Matrix) (float64, error)` | DistanceFrobenius returns the Frobenius norm of m - b, a metric on the space of matrices of equal shape. |
| `func (m *Matrix) DistanceToSingularity() float64` | DistanceToSingularity returns the 2-norm distance from a square matrix to the nearest singular matrix, which equals the smallest singular value. |
| `func (m *Matrix) Eigen() ([]Eigenpair, error)` | Eigen returns eigenvalue/eigenvector pairs for a square matrix. |
| `func (m *Matrix) EigenvalueMultiplicities(tol float64) (values []complex128, mult []int, err error)` | EigenvalueMultiplicities groups the eigenvalues of a square matrix into clusters that agree to within tol and returns the distinct representatives… |
| `func (m *Matrix) Eigenvalues() ([]complex128, error)` | Eigenvalues returns all eigenvalues of a square matrix, counted with multiplicity, computed with the shifted QR algorithm. |
| `func (m *Matrix) Equal(other *Matrix, tol float64) bool` | Equal reports whether m and other have the same shape and agree entrywise to within tol. |
| `func (m *Matrix) Exp() (*Matrix, error)` | Exp returns the matrix exponential exp(m) computed with the scaling and squaring method combined with a truncated Taylor series. |
| `func (m *Matrix) FrobeniusNorm() float64` | FrobeniusNorm returns the Frobenius norm, the square root of the sum of the squared moduli of the entries. |
| `func (m *Matrix) HadamardProduct(b *Matrix) (*Matrix, error)` | HadamardProduct returns the entrywise (Schur) product of m and b. |
| `func (m *Matrix) HermitianEigen() (values []float64, vectors *Matrix, err error)` | HermitianEigen returns the real eigenvalues (in ascending order) and the orthonormal eigenvectors (as the columns of the returned matrix) of a… |
| `func (m *Matrix) HermitianPart() *Matrix` | HermitianPart returns (m + m^H)/2, the Hermitian part of a square matrix. |
| `func (m *Matrix) Hessenberg() (h, q *Matrix, err error)` | Hessenberg returns the upper-Hessenberg reduction H together with the unitary Q such that m = Q*H*Q^H. It returns ErrNotSquare for a non-square… |
| `func (m *Matrix) ImagPart() *Matrix` | ImagPart returns the entrywise imaginary part of m as a real-valued complex matrix. |
| `func (m *Matrix) InPseudospectrum(z complex128, eps float64) bool` | InPseudospectrum reports whether z lies in the eps-pseudospectrum of m, i.e. |
| `func (m *Matrix) Inertia(tol float64) (neg, zero, pos int)` | Inertia returns the numbers of negative, zero and positive eigenvalues of a Hermitian matrix, using tol to decide which eigenvalues count as zero. |
| `func (m *Matrix) InfNorm() float64` | InfNorm returns the induced infinity-norm, the maximum absolute row sum. |
| `func (m *Matrix) Inverse() (*Matrix, error)` | Inverse returns the inverse of a square matrix. |
| `func (m *Matrix) IsContraction(tol float64) bool` | IsContraction reports whether the operator norm of m is at most 1 + tol. |
| `func (m *Matrix) IsDiagonal(tol float64) bool` | IsDiagonal reports whether all off-diagonal entries vanish to within tol. |
| `func (m *Matrix) IsHermitian(tol float64) bool` | IsHermitian reports whether the matrix equals its conjugate transpose to within tol. |
| `func (m *Matrix) IsIdempotent(tol float64) bool` | IsIdempotent is a synonym for IsProjection. |
| `func (m *Matrix) IsInvolution(tol float64) bool` | IsInvolution reports whether m^2 = I to within tol. |
| `func (m *Matrix) IsIsometry(tol float64) bool` | IsIsometry reports whether m preserves norms, i.e. |
| `func (m *Matrix) IsLowerTriangular(tol float64) bool` | IsLowerTriangular reports whether all entries above the main diagonal vanish to within tol. |
| `func (m *Matrix) IsNilpotent(tol float64) bool` | IsNilpotent reports whether some power m^k (k <= n) is the zero matrix to within tol. |
| `func (m *Matrix) IsNormal(tol float64) bool` | IsNormal reports whether m commutes with its adjoint, m^H m = m m^H, to within tol. |
| `func (m *Matrix) IsOrthogonal(tol float64) bool` | IsOrthogonal reports whether a real matrix satisfies m^T m = I to within tol. |
| `func (m *Matrix) IsOrthogonalProjection(tol float64) bool` | IsOrthogonalProjection reports whether m is both Hermitian and idempotent, so that it projects orthogonally onto its range. |
| `func (m *Matrix) IsPartialIsometry(tol float64) bool` | IsPartialIsometry reports whether m^H m is an orthogonal projection, the defining property of a partial isometry. |
| `func (m *Matrix) IsPositiveDefinite(tol float64) bool` | IsPositiveDefinite reports whether m is Hermitian with strictly positive eigenvalues, all exceeding tol. |
| `func (m *Matrix) IsPositiveSemidefinite(tol float64) bool` | IsPositiveSemidefinite reports whether m is Hermitian with eigenvalues no smaller than -tol. |
| `func (m *Matrix) IsProjection(tol float64) bool` | IsProjection reports whether m is idempotent, m^2 = m, to within tol. |
| `func (m *Matrix) IsSchurStable(tol float64) bool` | IsSchurStable reports whether every eigenvalue of m lies strictly inside the unit disc (a discrete-time stable operator). |
| `func (m *Matrix) IsSelfAdjoint(tol float64) bool` | IsSelfAdjoint is a synonym for IsHermitian. |
| `func (m *Matrix) IsSkewHermitian(tol float64) bool` | IsSkewHermitian reports whether m^H = -m to within tol. |
| `func (m *Matrix) IsSquare() bool` | IsSquare reports whether the matrix has equal numbers of rows and columns. |
| `func (m *Matrix) IsStable(tol float64) bool` | IsStable reports whether every eigenvalue of m has strictly negative real part (a continuous-time stable, or Hurwitz, operator). |
| `func (m *Matrix) IsStrictContraction(tol float64) bool` | IsStrictContraction reports whether the operator norm of m is strictly less than 1 - tol. |
| `func (m *Matrix) IsSymmetric(tol float64) bool` | IsSymmetric reports whether m equals its transpose (no conjugation) to within tol. |
| `func (m *Matrix) IsUnitary(tol float64) bool` | IsUnitary reports whether m^H m equals the identity to within tol. |
| `func (m *Matrix) IsUpperTriangular(tol float64) bool` | IsUpperTriangular reports whether all entries below the main diagonal vanish to within tol. |
| `func (m *Matrix) KernelBasis(tol float64) []Vector` | KernelBasis returns an orthonormal basis for the null space (kernel) of a square matrix, obtained from the right singular vectors with negligible… |
| `func (m *Matrix) Kron(b *Matrix) *Matrix` | Kron returns the Kronecker product m (x) b. |
| `func (m *Matrix) KyFanNorm(k int) float64` | KyFanNorm returns the Ky Fan k-norm, the sum of the k largest singular values. |
| `func (m *Matrix) Log() (*Matrix, error)` | Log returns a matrix logarithm of m via the functional calculus, so that exp(Log(m)) = m for a diagonalisable m with eigenvalues off the negative… |
| `func (m *Matrix) MaxAbs() float64` | MaxAbs returns the largest modulus among the entries of m. |
| `func (m *Matrix) MaxEigenvalue() float64` | MaxEigenvalue returns the largest eigenvalue of a Hermitian matrix. |
| `func (m *Matrix) MaxNorm() float64` | MaxNorm returns the largest modulus among the entries (the entrywise max-norm). |
| `func (m *Matrix) MinEigenvalue() float64` | MinEigenvalue returns the smallest eigenvalue of a Hermitian matrix. |
| `func (m *Matrix) Mul(b *Matrix) (*Matrix, error)` | Mul returns the matrix product m*b. |
| `func (m *Matrix) MulVec(v Vector) (Vector, error)` | MulVec returns the matrix-vector product m*v. |
| `func (m *Matrix) Neg() *Matrix` | Neg returns -m. |
| `func (m *Matrix) NuclearNorm() float64` | NuclearNorm returns the nuclear (trace) norm, the sum of the singular values. |
| `func (m *Matrix) NumericalAbscissa() float64` | NumericalAbscissa returns the maximum real part of the numerical range, equal to the largest eigenvalue of the Hermitian part (m + m^H)/2. |
| `func (m *Matrix) NumericalRadius(samples int) float64` | NumericalRadius returns the numerical radius w(m) = max { \|z\| : z in W(m) }, estimated by sampling the field of values at the given number of… |
| `func (m *Matrix) NumericalRange(samples int) []complex128` | NumericalRange returns a convex polygon approximating the boundary of the numerical range (field of values) W(m) = { <x, m x> : \|\|x\|\| = 1 },… |
| `func (m *Matrix) OneNorm() float64` | OneNorm returns the induced 1-norm, the maximum absolute column sum. |
| `func (m *Matrix) OperatorNorm() float64` | OperatorNorm returns the induced 2-norm (spectral norm), the largest singular value. |
| `func (m *Matrix) PolarDecomposition() (unitary, positive *Matrix, err error)` | PolarDecomposition returns the polar factors of a square matrix, U unitary (or a partial isometry when singular) and P Hermitian positive… |
| `func (m *Matrix) Power(k int) *Matrix` | Power returns m raised to the non-negative integer power k using binary exponentiation. |
| `func (m *Matrix) PowerReal(p float64) (*Matrix, error)` | PowerReal returns m raised to a real power p via the functional calculus, exp(p*Log(m)). |
| `func (m *Matrix) PseudospectralAbscissa(eps float64, samples int) (float64, error)` | PseudospectralAbscissa estimates the eps-pseudospectral abscissa, the maximum real part of any point in the eps-pseudospectrum. |
| `func (m *Matrix) PseudospectralRadius(eps float64, samples int) (float64, error)` | PseudospectralRadius estimates the eps-pseudospectral radius, the largest modulus of any point in the eps-pseudospectrum, by scanning a grid… |
| `func (m *Matrix) PseudospectrumGrid(reMin, reMax, imMin, imMax float64, nx, ny int) ([][]float64, error)` | PseudospectrumGrid evaluates the resolvent norm on a rectangular grid in the complex plane spanning [reMin,reMax] x [imMin,imMax] with nx columns and… |
| `func (m *Matrix) QR() (q, r *Matrix)` | QR returns a reduced QR factorisation m = Q*R computed with complex Householder reflections. |
| `func (m *Matrix) QuadraticForm(x, y Vector) complex128` | QuadraticForm returns the sesquilinear form <x, m y>. |
| `func (m *Matrix) RangeProjector() *Matrix` | RangeProjector returns the orthogonal projection onto the column space (range) of m. |
| `func (m *Matrix) Rank(tol float64) int` | Rank returns the numerical rank of the matrix, the number of singular values exceeding tol * (largest singular value). |
| `func (m *Matrix) RayleighQuotient(v Vector) complex128` | RayleighQuotient returns <v, m v> / <v, v> for a nonzero vector v. |
| `func (m *Matrix) RealPart() *Matrix` | RealPart returns the entrywise real part of m as a complex matrix with zero imaginary part. |
| `func (m *Matrix) Resolvent(z complex128) (*Matrix, error)` | Resolvent returns the resolvent operator (z*I - m)^{-1} at the complex point z. |
| `func (m *Matrix) ResolventNorm(z complex128) float64` | ResolventNorm returns the spectral norm of the resolvent at z, \|\|(z*I - m)^{-1}\|\|_2 = 1 / sigma_min(z*I - m). |
| `func (m *Matrix) Row(i int) Vector` | Row returns a copy of row i as a Vector. |
| `func (m *Matrix) Rows() int` | Rows returns the number of rows. |
| `func (m *Matrix) SVD() (u *Matrix, s []float64, v *Matrix)` | SVD returns a reduced singular value decomposition m = U * diag(s) * V^H, where U has orthonormal columns (m-by-k), s holds the k = min(rows,cols)… |
| `func (m *Matrix) Scale(s complex128) *Matrix` | Scale returns the matrix s*m. |
| `func (m *Matrix) ScaleReal(s float64) *Matrix` | ScaleReal returns the matrix s*m for a real scalar s. |
| `func (m *Matrix) SchattenNorm(p float64) float64` | SchattenNorm returns the Schatten p-norm, the l^p norm of the vector of singular values, for p >= 1. |
| `func (m *Matrix) Semigroup(t float64) (*Matrix, error)` | Semigroup returns exp(t*m), the value at time t of the one-parameter operator semigroup generated by m. |
| `func (m *Matrix) Set(i, j int, v complex128)` | Set stores v in row i and column j. |
| `func (m *Matrix) SetCol(j int, v Vector)` | SetCol overwrites column j with the entries of v. |
| `func (m *Matrix) SetRow(i int, v Vector)` | SetRow overwrites row i with the entries of v. |
| `func (m *Matrix) Sign() (*Matrix, error)` | Sign returns the matrix sign function, mapping each eigenvalue lambda to sign(Re lambda) (and 0 when Re lambda is 0). |
| `func (m *Matrix) Sin() (*Matrix, error)` | Sin returns the matrix sine sin(m) via the functional calculus. |
| `func (m *Matrix) SingularValues() []float64` | SingularValues returns the singular values of the matrix in descending order. |
| `func (m *Matrix) Sinh() (*Matrix, error)` | Sinh returns the matrix hyperbolic sine sinh(m) via the functional calculus. |
| `func (m *Matrix) SkewHermitianPart() *Matrix` | SkewHermitianPart returns (m - m^H)/2, the skew-Hermitian part of a square matrix. |
| `func (m *Matrix) SmallestSingularValue() float64` | SmallestSingularValue returns the smallest singular value of the matrix. |
| `func (m *Matrix) Solve(b *Matrix) (*Matrix, error)` | Solve returns the solution X of the linear system m*X = b for a matrix right-hand side. |
| `func (m *Matrix) SolveVec(b Vector) (Vector, error)` | SolveVec returns the solution x of the linear system m*x = b for a vector right-hand side. |
| `func (m *Matrix) SpectralAbscissa() float64` | SpectralAbscissa returns the maximum real part among the eigenvalues, which governs the growth rate of the operator semigroup exp(tA). |
| `func (m *Matrix) SpectralDecomposition(tol float64) ([]SpectralComponent, error)` | SpectralDecomposition returns the spectral resolution of a Hermitian matrix: one SpectralComponent per distinct eigenvalue (eigenvalues within tol… |
| `func (m *Matrix) SpectralGap() float64` | SpectralGap returns the difference between the two smallest distinct eigenvalues of a Hermitian matrix (the gap above the ground state). |
| `func (m *Matrix) SpectralNorm() float64` | SpectralNorm is an alias for OperatorNorm. |
| `func (m *Matrix) SpectralProjector(target, tol float64) *Matrix` | SpectralProjector returns the orthogonal projector onto the eigenspace of a Hermitian matrix associated with the eigenvalue closest to target,… |
| `func (m *Matrix) SpectralRadius() float64` | SpectralRadius returns the largest modulus among the eigenvalues. |
| `func (m *Matrix) Spectrum() ([]complex128, error)` | Spectrum returns the eigenvalues of a square matrix sorted lexicographically by real part and then imaginary part. |
| `func (m *Matrix) Sqrt() (*Matrix, error)` | Sqrt returns the principal square root of m via the functional calculus, so that Sqrt(m)^2 = m for a diagonalisable m. |
| `func (m *Matrix) String() string` | String renders the matrix with each entry formatted to three decimal places. |
| `func (m *Matrix) Sub(b *Matrix) (*Matrix, error)` | Sub returns m - b. |
| `func (m *Matrix) Submatrix(r1, r2, c1, c2 int) (*Matrix, error)` | Submatrix returns the r1..r2-1 by c1..c2-1 block of m. |
| `func (m *Matrix) Trace() complex128` | Trace returns the sum of the diagonal entries. |
| `func (m *Matrix) TraceNorm() float64` | TraceNorm is an alias for NuclearNorm. |
| `func (m *Matrix) Transpose() *Matrix` | Transpose returns the transpose of m (no conjugation). |

</details>

<details>
<summary><code>Vector</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BasisVector(n, i int) Vector` | BasisVector returns the i-th standard basis vector of C^n. |
| `func GramSchmidt(vectors []Vector, tol float64) []Vector` | GramSchmidt orthonormalises the supplied vectors using the modified Gram-Schmidt process with the Hermitian inner product. |
| `func NewVector(n int) Vector` | NewVector returns a zero vector of length n. |
| `func VectorFromReal(data []float64) Vector` | VectorFromReal returns a complex vector whose entries are the given real numbers with zero imaginary part. |
| `func VectorOf(entries ...complex128) Vector` | VectorOf returns a Vector containing the supplied entries. |
| `func (v Vector) Add(w Vector) Vector` | Add returns v + w. |
| `func (v Vector) Angle(w Vector) float64` | Angle returns the angle in radians between the real directions of v and w, computed from \|<v,w>\| / (\|\|v\|\| \|\|w\|\|). |
| `func (v Vector) Clone() Vector` | Clone returns an independent copy of v. |
| `func (v Vector) Conjugate() Vector` | Conjugate returns the entrywise complex conjugate of v. |
| `func (v Vector) Dot(w Vector) complex128` | Dot returns the Hermitian inner product <v,w> = sum conj(v_i) w_i. |
| `func (v Vector) Equal(w Vector, tol float64) bool` | Equal reports whether v and w have the same length and agree entrywise to within tol. |
| `func (v Vector) InfNorm() float64` | InfNorm returns the l-infinity norm, the largest modulus among the entries. |
| `func (v Vector) IsOrthogonal(w Vector, tol float64) bool` | IsOrthogonal reports whether v and w are orthogonal to within tol, i.e. |
| `func (v Vector) IsZero(tol float64) bool` | IsZero reports whether every entry of v has modulus at most tol. |
| `func (v Vector) Len() int` | Len returns the dimension of the vector. |
| `func (v Vector) Neg() Vector` | Neg returns -v. |
| `func (v Vector) Norm() float64` | Norm returns the Euclidean (l2) norm of v. |
| `func (v Vector) NormSquared() float64` | NormSquared returns the squared Euclidean norm <v,v>. |
| `func (v Vector) Normalize() (Vector, float64)` | Normalize returns a unit vector in the direction of v together with the original norm. |
| `func (v Vector) OneNorm() float64` | OneNorm returns the l1 norm, the sum of the moduli of the entries. |
| `func (v Vector) OuterProduct(w Vector) *Matrix` | OuterProduct returns the rank-one operator v w^H, whose (i,j) entry is v_i conj(w_j). |
| `func (v Vector) ProjectOnto(u Vector) Vector` | ProjectOnto returns the orthogonal projection of v onto the line spanned by the nonzero vector u, namely (<u,v>/<u,u>) u. |
| `func (v Vector) Scale(s complex128) Vector` | Scale returns the vector s*v. |
| `func (v Vector) Sub(w Vector) Vector` | Sub returns v - w. |

</details>

### Variables

`ErrDimensionMismatch`, `ErrNotSquare`, `ErrNotHermitian`, `ErrNotNormal`, `ErrSingular`, `ErrEmpty`, `ErrOutOfRange`, `ErrInvalidArgument`, `ErrNoConvergence`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/operatortheory).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
