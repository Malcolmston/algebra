# liealgebra — implements computational tools for Lie algebras and root systems using only the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/liealgebra.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/liealgebra)

Package liealgebra implements computational tools for Lie algebras and root
systems using only the Go standard library.

- Dense real and complex matrices (`Matrix`, `CMatrix`) with the arithmetic
needed for matrix Lie algebras: addition, multiplication, transpose, conjugate
transpose, trace, Frobenius norm, LU based determinant, inverse and linear
solve.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/liealgebra"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AdjointDimension(family string, rank int) (int, error)` | AdjointDimension returns the dimension of the adjoint representation, which equals the dimension of the Lie algebra itself. |
| `func AllRoots(family string, rank int) ([][]float64, error)` | AllRoots returns all roots (positive together with their negatives) of a classical type as Euclidean coordinate vectors. |
| `func CTraceForm(a, b *CMatrix) (complex128, error)` | CTraceForm returns the complex trace form tr(AB) of two complex matrices. |
| `func CartanInteger(alpha, beta []float64) (float64, error)` | CartanInteger returns the Cartan integer ⟨α,β∨⟩ = 2(α,β)/(β,β) of two Euclidean root vectors. |
| `func CartanMatrixDeterminant(family string, rank int) (float64, error)` | CartanMatrixDeterminant returns the determinant of the Cartan matrix, an important invariant equal to the order of the center of the simply connected… |
| `func CasimirEigenvalue(family string, rank int, dynkin []int) (float64, error)` | CasimirEigenvalue returns the eigenvalue of the quadratic Casimir operator on the irreducible representation with the given highest weight (Dynkin… |
| `func CasimirSU2(j float64) float64` | CasimirSU2 returns the physics-convention quadratic Casimir eigenvalue j(j+1) of the spin-j irreducible representation of su(2). |
| `func CasimirTraceNormalization(basis []*Matrix) (float64, error)` | CasimirTraceNormalization returns the constant c such that the Killing form equals c times the trace form on the defining representation, given a… |
| `func Coroot(root []float64) ([]float64, error)` | Coroot returns the coroot α∨ = 2α/(α,α) of a Euclidean root vector. |
| `func CoxeterNumber(family string, rank int) (int, error)` | CoxeterNumber returns the Coxeter number h of the given Dynkin type, equal to (number of roots)/rank and to one plus the height of the highest root. |
| `func Det(m *Matrix) (float64, error)` | Det returns the determinant of a square matrix. |
| `func DimensionSLn(n int) (int, error)` | DimensionSLn returns the dimension n^2-1 of the Lie algebra sl(n) = A_{n-1}. |
| `func DimensionSOn(n int) (int, error)` | DimensionSOn returns the dimension n(n-1)/2 of the Lie algebra so(n). |
| `func DimensionSPn(n int) (int, error)` | DimensionSPn returns the dimension n(2n+1) of the Lie algebra sp(2n) = C_n. |
| `func DimensionSU2(j float64) (int, error)` | DimensionSU2 returns the dimension 2j+1 of the spin-j representation of su(2). |
| `func DimensionSUn(n int) (int, error)` | DimensionSUn returns the dimension n^2-1 of the Lie algebra su(n). |
| `func DualCoxeterFromRoots(family string, rank int) (int, error)` | DualCoxeterFromRoots returns the dual Coxeter number computed as 1 + (ρ, θ∨) where θ is the highest root; provided as an independent cross-check… |
| `func DualCoxeterNumber(family string, rank int) (int, error)` | DualCoxeterNumber returns the dual Coxeter number h∨ of the given Dynkin type. |
| `func FundamentalWeightsRootBasis(family string, rank int) ([][]float64, error)` | FundamentalWeightsRootBasis returns the fundamental weights expressed in the basis of simple roots; row i is ω_i. |
| `func HighestRootLabel(family string, rank int) ([]int, error)` | HighestRootLabel returns the coefficient vector (in the simple-root basis) of the highest root, the unique positive root of maximal height. |
| `func IsAbelian(basis []*Matrix, tol float64) bool` | IsAbelian reports whether every pair of basis elements commutes to within tolerance tol, i.e. |
| `func IsClosedUnderBracket(basis []*Matrix, tol float64) bool` | IsClosedUnderBracket reports whether the bracket of every pair of basis elements lies in the span of the basis to within tolerance tol, i.e. |
| `func IsNilpotentMatrix(m *Matrix, tol float64) (bool, error)` | IsNilpotentMatrix reports whether the square matrix m is nilpotent, i.e. |
| `func IsSemisimple(basis []*Matrix, tol float64) (bool, error)` | IsSemisimple reports whether the Killing form of the given basis is nondegenerate (Cartan's criterion), tested by \|det K\| > tol. |
| `func IsSimplyLaced(family string, rank int) (bool, error)` | IsSimplyLaced reports whether the Dynkin type is simply laced, i.e. |
| `func JacobiResidualConstants(c [][][]float64) float64` | JacobiResidualConstants returns the maximum absolute value of the Jacobi identity applied to structure constants: Σ_m ( c`m``i``j` c`n``m``k` +… |
| `func KillingFormRank(basis []*Matrix, tol float64) (int, error)` | KillingFormRank returns the rank of the Killing form of a basis. |
| `func KillingFormValue(basis []*Matrix, x, y []float64) (float64, error)` | KillingFormValue returns the Killing form value K(X,Y) = tr(ad_X ad_Y) for two elements X and Y given by their coordinate vectors in the supplied… |
| `func LieAlgebraDimension(family string, rank int) (int, error)` | LieAlgebraDimension returns the dimension of the semisimple Lie algebra of the given type, equal to rank + number of roots. |
| `func LieAlgebraName(family string, rank int) (string, error)` | LieAlgebraName returns a human-readable classical name for the given Dynkin type, such as "sl(4)" for A_3 or "so(7)" for B_3. |
| `func NumPositiveRoots(family string, rank int) (int, error)` | NumPositiveRoots returns the number of positive roots of the given type. |
| `func NumRoots(family string, rank int) (int, error)` | NumRoots returns the total number of roots (positive and negative). |
| `func PositiveRootLabels(family string, rank int) ([][]int, error)` | PositiveRootLabels returns every positive root of the given type expressed as an integer coefficient vector in the basis of simple roots. |
| `func PositiveRoots(family string, rank int) ([][]float64, error)` | PositiveRoots returns the positive roots of a classical type as Euclidean coordinate vectors. |
| `func Rank(m *Matrix, tol float64) int` | Rank returns the numerical rank of a matrix using Gaussian elimination with the given tolerance for treating a pivot as zero. |
| `func RootInnerProduct(b *Matrix, x, y []float64) (float64, error)` | RootInnerProduct returns the inner product of two vectors expressed in the simple-root basis using the bilinear form B (see `RootBilinearForm`). |
| `func SO3ToVector(m *Matrix) ([]float64, error)` | SO3ToVector maps a 3x3 antisymmetric matrix to its axial vector (w) such that the matrix acts as w×. |
| `func SatisfiesJacobi(a, b, c *Matrix, tol float64) bool` | SatisfiesJacobi reports whether the Jacobi identity holds for A, B, C to within tolerance tol (measured by the max absolute entry of the residual). |
| `func SimpleRoots(family string, rank int) ([][]float64, error)` | SimpleRoots returns the simple roots of a classical type (A, B, C or D) as Euclidean coordinate vectors. |
| `func Solve(a *Matrix, b []float64) ([]float64, error)` | Solve solves the linear system A x = b for a square A. It returns `ErrSingular` if A is not invertible and `ErrDim` on a length mismatch. |
| `func SolveLeastSquares(a *Matrix, b []float64) ([]float64, error)` | SolveLeastSquares solves the (possibly overdetermined) system A x = b in the least-squares sense via the normal equations AᵀA x = Aᵀb. |
| `func SpectralRadiusBound(m *Matrix) float64` | SpectralRadiusBound returns an upper bound on the spectral radius of a real matrix via the maximum absolute row sum (an induced ∞-norm). |
| `func StructureConstants(basis []*Matrix) ([][][]float64, error)` | StructureConstants computes the structure constants c`k``i``j` of a matrix Lie algebra with the given basis, defined by [e_i, e_j] = Σ_k c`k``i``j`… |
| `func TraceForm(a, b *Matrix) (float64, error)` | TraceForm returns the trace form tr(AB) of two square matrices. |
| `func VecAdd(a, b []float64) ([]float64, error)` | VecAdd returns a+b, or `ErrDim` on a length mismatch. |
| `func VecDot(a, b []float64) (float64, error)` | VecDot returns the Euclidean dot product of two equal-length vectors, or `ErrDim` on a length mismatch. |
| `func VecEqual(a, b []float64, tol float64) bool` | VecEqual reports whether two vectors are equal to within tolerance tol. |
| `func VecNorm(a []float64) float64` | VecNorm returns the Euclidean norm of a vector. |
| `func VecNormSquared(a []float64) float64` | VecNormSquared returns the squared Euclidean norm of a vector. |
| `func VecScale(a []float64, s float64) []float64` | VecScale returns the vector scaled by s. |
| `func VecSub(a, b []float64) ([]float64, error)` | VecSub returns a-b, or `ErrDim` on a length mismatch. |
| `func WeylDimension(family string, rank int, dynkin []int) (int, error)` | WeylDimension returns the dimension of the irreducible representation of the given Dynkin type with highest weight specified by nonnegative Dynkin… |
| `func WeylGroupOrder(family string, rank int) (*big.Int, error)` | WeylGroupOrder returns the order of the Weyl group of the given Dynkin type as an arbitrary-precision integer. |
| `func WeylOrbit(simple [][]float64, v []float64, tol float64) ([][]float64, error)` | WeylOrbit returns the orbit of a Euclidean vector v under the Weyl group generated by reflections in the given simple roots. |
| `func WeylReflection(root, v []float64) ([]float64, error)` | WeylReflection returns the reflection of the vector v across the hyperplane orthogonal to root, s_root(v) = v - 2(v,root)/(root,root) root, in… |
| `func WeylVectorRootBasis(family string, rank int) ([]float64, error)` | WeylVectorRootBasis returns the Weyl vector ρ (the sum of the fundamental weights, equivalently half the sum of the positive roots) expressed in the… |

### Types

| Type | What it is |
| --- | --- |
| `CMatrix` | CMatrix is a dense complex matrix stored in row-major order. |
| `Matrix` | Matrix is a dense real matrix stored in row-major order. |

<details>
<summary><code>CMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CAnticommutator(a, b *CMatrix) (*CMatrix, error)` | CAnticommutator returns {A,B} = AB + BA of two complex matrices. |
| `func CBracket(a, b *CMatrix) (*CMatrix, error)` | CBracket returns the Lie bracket [A,B] = AB - BA of two complex matrices. |
| `func CMatExp(a *CMatrix) (*CMatrix, error)` | CMatExp returns the matrix exponential exp(A) of a square complex matrix by scaling and squaring with a truncated Taylor series. |
| `func CMatPow(a *CMatrix, p int) (*CMatrix, error)` | CMatPow returns a complex matrix raised to a nonnegative integer power. |
| `func DiagCMatrix(d []complex128) *CMatrix` | DiagCMatrix returns the square complex matrix with the given diagonal. |
| `func GellMannMatrices() []*CMatrix` | GellMannMatrices returns all eight Gell-Mann matrices λ_1..λ_8. |
| `func GellMannMatrix(i int) *CMatrix` | GellMannMatrix returns the i-th Gell-Mann matrix, i in 1..8. |
| `func IdentityCMatrix(n int) *CMatrix` | IdentityCMatrix returns the n-by-n complex identity matrix. |
| `func NewCMatrix(rows, cols int) *CMatrix` | NewCMatrix returns a rows-by-cols zero complex matrix. |
| `func NewCMatrixFromRows(rows [][]complex128) (*CMatrix, error)` | NewCMatrixFromRows builds a complex matrix from equal-length rows, returning `ErrDim` if the rows are ragged. |
| `func PauliMatrices() (sx, sy, sz *CMatrix)` | PauliMatrices returns the three Pauli matrices (σ_x, σ_y, σ_z). |
| `func PauliX() *CMatrix` | PauliX returns the Pauli matrix σ_x = [[0,1],[1,0]]. |
| `func PauliY() *CMatrix` | PauliY returns the Pauli matrix σ_y = [[0,-i],[i,0]]. |
| `func PauliZ() *CMatrix` | PauliZ returns the Pauli matrix σ_z = [[1,0],[0,-1]]. |
| `func RealToComplex(m *Matrix) *CMatrix` | RealToComplex promotes a real matrix to a complex matrix. |
| `func SU2Generators() (t1, t2, t3 *CMatrix)` | SU2Generators returns the anti-Hermitian generators of su(2) given by T_i = -i σ_i / 2, satisfying [T_i, T_j] = ε_ijk T_k. |
| `func SU2SpinMatrices() (sx, sy, sz *CMatrix)` | SU2SpinMatrices returns the Hermitian spin-1/2 operators (S_x, S_y, S_z) = (σ_x, σ_y, σ_z)/2 satisfying [S_i,S_j] = i ε_ijk S_k. |
| `func SU3Generators() []*CMatrix` | SU3Generators returns the anti-Hermitian su(3) generators T_a = -i λ_a / 2 for a in 1..8. |
| `func SpinJz(j float64) (*CMatrix, error)` | SpinJz returns the diagonal Jz operator for spin j. |
| `func SpinLowering(j float64) (*CMatrix, error)` | SpinLowering returns the lowering operator J- for spin j. |
| `func SpinMatrices(j float64) (jx, jy, jz *CMatrix, err error)` | SpinMatrices returns the (2j+1)-dimensional Hermitian angular-momentum operators (Jx, Jy, Jz) for spin/quantum number j (a nonnegative half-integer… |
| `func SpinRaising(j float64) (*CMatrix, error)` | SpinRaising returns the raising operator J+ for spin j. |
| `func (m *CMatrix) Add(b *CMatrix) (*CMatrix, error)` | Add returns m+b or `ErrDim`. |
| `func (m *CMatrix) AntiHermitianPart() (*CMatrix, error)` | AntiHermitianPart returns (m-mᴴ)/2, the anti-Hermitian part. |
| `func (m *CMatrix) ApproxEqual(b *CMatrix, tol float64) bool` | ApproxEqual reports shape equality with every entry agreeing to within tol in complex modulus. |
| `func (m *CMatrix) At(i, j int) complex128` | At returns element (i,j), panicking on an out-of-range index. |
| `func (m *CMatrix) Clone() *CMatrix` | Clone returns a deep copy. |
| `func (m *CMatrix) Conjugate() *CMatrix` | Conjugate returns the entrywise complex conjugate. |
| `func (m *CMatrix) ConjugateTranspose() *CMatrix` | ConjugateTranspose returns the Hermitian conjugate (dagger) mᴴ. |
| `func (m *CMatrix) Dagger() *CMatrix` | Dagger is an alias for `CMatrix.ConjugateTranspose`. |
| `func (m *CMatrix) Dims() (int, int)` | Dims returns the number of rows and columns. |
| `func (m *CMatrix) Equal(b *CMatrix) bool` | Equal reports exact equality of shape and entries. |
| `func (m *CMatrix) FrobeniusNorm() float64` | FrobeniusNorm returns the square root of the sum of squared moduli. |
| `func (m *CMatrix) HermitianPart() (*CMatrix, error)` | HermitianPart returns (m+mᴴ)/2, the Hermitian part of a square complex matrix. |
| `func (m *CMatrix) Imag() *Matrix` | Imag returns the imaginary part of each entry as a real matrix. |
| `func (m *CMatrix) IsAntiHermitian(tol float64) bool` | IsAntiHermitian reports whether m equals the negative of its conjugate transpose within tol (a skew-Hermitian matrix, as in su(n)). |
| `func (m *CMatrix) IsHermitian(tol float64) bool` | IsHermitian reports whether m equals its conjugate transpose within tol. |
| `func (m *CMatrix) IsSquare() bool` | IsSquare reports whether the matrix is square. |
| `func (m *CMatrix) IsTraceless(tol float64) bool` | IsTraceless reports whether the trace is zero within tol. |
| `func (m *CMatrix) IsUnitary(tol float64) bool` | IsUnitary reports whether mᴴm is the identity within tol. |
| `func (m *CMatrix) Kronecker(b *CMatrix) *CMatrix` | Kronecker returns the Kronecker (tensor) product m ⊗ b of complex matrices. |
| `func (m *CMatrix) MaxAbs() float64` | MaxAbs returns the largest entry modulus. |
| `func (m *CMatrix) Mul(b *CMatrix) (*CMatrix, error)` | Mul returns the matrix product m*b or `ErrDim`. |
| `func (m *CMatrix) Neg() *CMatrix` | Neg returns the additive inverse. |
| `func (m *CMatrix) Real() *Matrix` | Real returns the real part of each entry as a real matrix. |
| `func (m *CMatrix) Scale(s complex128) *CMatrix` | Scale returns the matrix times the complex scalar s. |
| `func (m *CMatrix) Set(i, j int, v complex128)` | Set assigns v to element (i,j), panicking on an out-of-range index. |
| `func (m *CMatrix) Sub(b *CMatrix) (*CMatrix, error)` | Sub returns m-b or `ErrDim`. |
| `func (m *CMatrix) Trace() (complex128, error)` | Trace returns the sum of diagonal entries, or `ErrNotSquare`. |
| `func (m *CMatrix) Transpose() *CMatrix` | Transpose returns the (unconjugated) transpose. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AdjointAction(x, y *Matrix) (*Matrix, error)` | AdjointAction returns ad_X(Y) = [X,Y], the action of X in the adjoint representation on the matrix Y. |
| `func AdjointMatrix(basis []*Matrix, index int) (*Matrix, error)` | AdjointMatrix returns the matrix of ad_{e_index} in the given basis: the d-by-d matrix M with M`k``j` = c`k``index``j` where c are the structure… |
| `func Anticommutator(a, b *Matrix) (*Matrix, error)` | Anticommutator returns {A,B} = AB + BA of two square real matrices. |
| `func BCHApprox(x, y *Matrix) (*Matrix, error)` | BCHApprox returns the second-order Baker-Campbell-Hausdorff approximation to log(exp(X)exp(Y)): |
| `func Bracket(a, b *Matrix) (*Matrix, error)` | Bracket returns the Lie bracket (commutator) [A,B] = AB - BA of two square real matrices. |
| `func CartanMatrix(family string, rank int) (*Matrix, error)` | CartanMatrix returns the Cartan matrix of the Lie algebra of the given Dynkin type. |
| `func CartanMatrixA(n int) (*Matrix, error)` | CartanMatrixA returns the Cartan matrix of type A_n (the Lie algebra sl(n+1)). |
| `func CartanMatrixB(n int) (*Matrix, error)` | CartanMatrixB returns the Cartan matrix of type B_n (the Lie algebra so(2n+1)). |
| `func CartanMatrixC(n int) (*Matrix, error)` | CartanMatrixC returns the Cartan matrix of type C_n (the Lie algebra sp(2n)). |
| `func CartanMatrixD(n int) (*Matrix, error)` | CartanMatrixD returns the Cartan matrix of type D_n (the Lie algebra so(2n)). |
| `func CartanMatrixFromRoots(simple [][]float64) (*Matrix, error)` | CartanMatrixFromRoots builds the Cartan matrix A`i``j` = 2(α_i,α_j)/(α_j,α_j) from a list of simple roots given as Euclidean coordinate vectors. |
| `func DiagMatrix(d []float64) *Matrix` | DiagMatrix returns the square matrix with the given diagonal entries. |
| `func DynkinBondMatrix(family string, rank int) (*Matrix, error)` | DynkinBondMatrix returns the matrix B`i``j` = A`i``j`*A`j``i`, whose value (0,1,2,3) counts the number of bonds between nodes i and j in the Dynkin… |
| `func DynkinDiagramAdjacency(family string, rank int) (*Matrix, error)` | DynkinDiagramAdjacency returns the 0/1 adjacency matrix of the Dynkin diagram: nodes i and j are adjacent when the Cartan entry A`i``j` is nonzero. |
| `func ExpBracketSeries(x, y *Matrix, terms int) (*Matrix, error)` | ExpBracketSeries returns the truncated series for Ad_{exp(X)}(Y) = exp(ad_X)(Y) = Σ_{k≥0} (1/k!) ad_X^k(Y), summed to the given number of terms… |
| `func ExpMap(a *Matrix) (*Matrix, error)` | ExpMap is an alias for `MatExp`; for a Lie-algebra element X it returns the corresponding group element exp(X). |
| `func HeisenbergGenerators() (x, y, z *Matrix)` | HeisenbergGenerators returns the three matrices (X, Y, Z) generating the 3-dimensional Heisenberg Lie algebra of strictly upper-triangular 3x3 real… |
| `func IdentityMatrix(n int) *Matrix` | IdentityMatrix returns the n-by-n identity matrix. |
| `func Inverse(m *Matrix) (*Matrix, error)` | Inverse returns the matrix inverse, or `ErrSingular` if the matrix is not invertible. |
| `func InverseCartanMatrix(family string, rank int) (*Matrix, error)` | InverseCartanMatrix returns the inverse of the Cartan matrix. |
| `func JacobiResidual(a, b, c *Matrix) (*Matrix, error)` | JacobiResidual returns [[A,B],C] + [[B,C],A] + [[C,A],B], which is identically the zero matrix for any associative product and hence for matrix Lie… |
| `func KillingForm(basis []*Matrix) (*Matrix, error)` | KillingForm returns the Killing form matrix K of a Lie algebra given by a basis, with entries K`i``j` = tr(ad_{e_i} ad_{e_j}). |
| `func MatExp(a *Matrix) (*Matrix, error)` | MatExp returns the matrix exponential exp(A) of a square real matrix, computed by scaling and squaring with a truncated Taylor series. |
| `func MatPow(a *Matrix, p int) (*Matrix, error)` | MatPow returns A raised to a nonnegative integer power using exponentiation by squaring; A^0 is the identity. |
| `func NestedBracket(xs ...*Matrix) (*Matrix, error)` | NestedBracket returns the left-nested bracket [x0,[x1,[...,xk]...]] for a list of matrices. |
| `func NewMatrix(rows, cols int) *Matrix` | NewMatrix returns a rows-by-cols zero matrix. |
| `func NewMatrixFromRows(rows [][]float64) (*Matrix, error)` | NewMatrixFromRows builds a matrix from a slice of equal-length rows. |
| `func RootBilinearForm(family string, rank int) (*Matrix, error)` | RootBilinearForm returns the symmetric matrix B with B`i``j` = (α_i,α_j), the inner products of the simple roots normalised so that the longest… |
| `func SL2CartanH() *Matrix` | SL2CartanH returns the Cartan generator H of sl(2). |
| `func SL2Generators() (e, f, h *Matrix)` | SL2Generators returns the standard basis (E, F, H) of the Lie algebra sl(2) of traceless 2x2 real matrices, satisfying [H,E]=2E, [H,F]=-2F and… |
| `func SL2LoweringF() *Matrix` | SL2LoweringF returns the lowering generator F of sl(2). |
| `func SL2RaisingE() *Matrix` | SL2RaisingE returns the raising generator E of sl(2). |
| `func SO3FromVector(w []float64) (*Matrix, error)` | SO3FromVector returns the 3x3 antisymmetric matrix `w`_× representing the cross product w×(·). |
| `func SO3Generators() (lx, ly, lz *Matrix)` | SO3Generators returns the three real antisymmetric generators (Lx, Ly, Lz) of the rotation algebra so(3), satisfying [Lx,Ly]=Lz and cyclic… |
| `func WeylReflectionMatrix(root []float64) (*Matrix, error)` | WeylReflectionMatrix returns the orthogonal matrix of the reflection across the hyperplane orthogonal to the given root in Euclidean coordinates. |
| `func ZeroMatrix(n int) *Matrix` | ZeroMatrix returns an n-by-n zero matrix. |
| `func (m *Matrix) Add(b *Matrix) (*Matrix, error)` | Add returns m+b or `ErrDim` on a shape mismatch. |
| `func (m *Matrix) AntiSymmetrize() (*Matrix, error)` | AntiSymmetrize returns (m-mᵀ)/2, the antisymmetric part of a square matrix. |
| `func (m *Matrix) ApproxEqual(b *Matrix, tol float64) bool` | ApproxEqual reports whether a and b have the same shape and every entry agrees to within absolute tolerance tol. |
| `func (m *Matrix) At(i, j int) float64` | At returns element (i,j). |
| `func (m *Matrix) Clone() *Matrix` | Clone returns a deep copy of the matrix. |
| `func (m *Matrix) Col(j int) []float64` | Col returns a copy of column j. |
| `func (m *Matrix) Dims() (int, int)` | Dims returns the number of rows and columns. |
| `func (m *Matrix) Equal(b *Matrix) bool` | Equal reports whether a and b have the same shape and identical entries. |
| `func (m *Matrix) FrobeniusNorm() float64` | FrobeniusNorm returns the square root of the sum of squared entries. |
| `func (m *Matrix) IsAntisymmetric(tol float64) bool` | IsAntisymmetric reports whether m equals the negative of its transpose within tolerance tol (a real antisymmetric/skew matrix). |
| `func (m *Matrix) IsDiagonal(tol float64) bool` | IsDiagonal reports whether all off-diagonal entries vanish within tol. |
| `func (m *Matrix) IsSquare() bool` | IsSquare reports whether the matrix has equal row and column counts. |
| `func (m *Matrix) IsSymmetric(tol float64) bool` | IsSymmetric reports whether m equals its transpose within tolerance tol. |
| `func (m *Matrix) IsTraceless(tol float64) bool` | IsTraceless reports whether the trace is zero within tolerance tol. |
| `func (m *Matrix) Kronecker(b *Matrix) *Matrix` | Kronecker returns the Kronecker (tensor) product m ⊗ b. |
| `func (m *Matrix) MatVec(x []float64) ([]float64, error)` | MatVec returns the matrix-vector product m*x or `ErrDim` on a mismatch. |
| `func (m *Matrix) MaxAbs() float64` | MaxAbs returns the largest absolute entry, or 0 for an empty matrix. |
| `func (m *Matrix) Mul(b *Matrix) (*Matrix, error)` | Mul returns the matrix product m*b or `ErrDim` if the inner dimensions disagree. |
| `func (m *Matrix) Neg() *Matrix` | Neg returns the additive inverse of the matrix. |
| `func (m *Matrix) Row(i int) []float64` | Row returns a copy of row i. |
| `func (m *Matrix) Scale(s float64) *Matrix` | Scale returns the matrix multiplied by the scalar s. |
| `func (m *Matrix) Set(i, j int, v float64)` | Set assigns v to element (i,j). |
| `func (m *Matrix) Sub(b *Matrix) (*Matrix, error)` | Sub returns m-b or `ErrDim` on a shape mismatch. |
| `func (m *Matrix) Symmetrize() (*Matrix, error)` | Symmetrize returns (m+mᵀ)/2, the symmetric part of a square matrix. |
| `func (m *Matrix) Trace() (float64, error)` | Trace returns the sum of the diagonal entries. |
| `func (m *Matrix) Transpose() *Matrix` | Transpose returns the transpose of the matrix. |
| `func (m *Matrix) Vec() []float64` | Vec flattens the matrix into a single row-major slice (a fresh copy). |

</details>

### Variables

`ErrDim`, `ErrNotSquare`, `ErrSingular`, `ErrRank`, `ErrType`, `ErrRange`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/liealgebra).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
