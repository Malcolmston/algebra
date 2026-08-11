# tropical — implements algebra over the tropical (min-plus and max-plus) semirings

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/tropical.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/tropical)

Package tropical implements algebra over the tropical (min-plus and max-plus)
semirings and the numerical objects built on top of them.

A semiring replaces ordinary addition and multiplication with two new
operations. In the min-plus semiring the "sum" of two numbers is their minimum
and the "product" is their ordinary sum; its neutral element for addition (the
tropical zero) is +Inf and its neutral element for multiplication (the
tropical one) is 0. The max-plus semiring is the order-dual: addition is the
maximum, the tropical zero is -Inf and the tropical one is still 0. Both are
idempotent, commutative, complete dioids and the two are exchanged by negating
every entry.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/tropical"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func IsInf(x float64) bool` | IsInf reports whether x is either positive or negative infinity, i.e. |
| `func MaxPlusAdd(a, b float64) float64` | MaxPlusAdd returns the max-plus tropical sum of a and b, that is max(a, b). |
| `func MaxPlusDiv(a, b float64) float64` | MaxPlusDiv returns the max-plus tropical quotient a/b (a-b on finite values). |
| `func MaxPlusMul(a, b float64) float64` | MaxPlusMul returns the max-plus tropical product of a and b (a+b, with -Inf absorbing). |
| `func MaxPlusOne() float64` | MaxPlusOne returns the tropical one of the max-plus semiring, 0. |
| `func MaxPlusPow(a float64, n int) float64` | MaxPlusPow returns the max-plus tropical power a^n (n*a on finite values). |
| `func MaxPlusProd(xs ...float64) float64` | MaxPlusProd returns the max-plus tropical product (ordinary sum) of the arguments, or 0 when none are given. |
| `func MaxPlusStar(a float64) float64` | MaxPlusStar returns the max-plus scalar Kleene star: 0 for a <= 0 and +Inf otherwise. |
| `func MaxPlusSum(xs ...float64) float64` | MaxPlusSum returns the max-plus tropical sum (maximum) of the arguments, or -Inf when none are given. |
| `func MaxPlusZero() float64` | MaxPlusZero returns the tropical zero of the max-plus semiring, -Inf. |
| `func MinPlusAdd(a, b float64) float64` | MinPlusAdd returns the min-plus tropical sum of a and b, that is min(a, b). |
| `func MinPlusDiv(a, b float64) float64` | MinPlusDiv returns the min-plus tropical quotient a/b (a-b on finite values). |
| `func MinPlusMul(a, b float64) float64` | MinPlusMul returns the min-plus tropical product of a and b (a+b, with +Inf absorbing). |
| `func MinPlusOne() float64` | MinPlusOne returns the tropical one of the min-plus semiring, 0. |
| `func MinPlusPow(a float64, n int) float64` | MinPlusPow returns the min-plus tropical power a^n (n*a on finite values). |
| `func MinPlusProd(xs ...float64) float64` | MinPlusProd returns the min-plus tropical product (ordinary sum) of the arguments, or 0 when none are given. |
| `func MinPlusStar(a float64) float64` | MinPlusStar returns the min-plus scalar Kleene star: 0 for a >= 0 and -Inf otherwise. |
| `func MinPlusSum(xs ...float64) float64` | MinPlusSum returns the min-plus tropical sum (minimum) of the arguments, or +Inf when none are given. |
| `func MinPlusZero() float64` | MinPlusZero returns the tropical zero of the min-plus semiring, +Inf. |

### Types

| Type | What it is |
| --- | --- |
| `Kind` | Kind enumerates the two tropical semirings supported by this package. |
| `Matrix` | Matrix is a dense tropical matrix stored in row-major order together with the semiring under which its entries are interpreted. |
| `Poly` | Poly is a univariate tropical polynomial. |
| `PolygonVertex` | PolygonVertex is a vertex of a Newton polygon: a coefficient index paired with its (finite) coefficient value. |
| `Root` | Root is a tropical root together with its multiplicity (the horizontal width of the Newton-polygon edge that produced it). |
| `Semiring` | Semiring is a value describing which tropical semiring a computation uses. |
| `Vector` | Vector is a dense tropical vector: a slice of scalar entries together with the semiring under which it is interpreted. |

<details>
<summary><code>Kind</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (k Kind) String() string` | String returns the human-readable name of the semiring kind. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Constant(sr Semiring, r, c int, value float64) Matrix` | Constant returns an r-by-c matrix with every entry equal to value. |
| `func Diag(sr Semiring, values []float64) Matrix` | Diag returns a square matrix with the given values on the diagonal and the tropical zero elsewhere. |
| `func Identity(sr Semiring, n int) Matrix` | Identity returns the n-by-n tropical identity matrix: the tropical one (0) on the diagonal and the tropical zero elsewhere. |
| `func MaxPlusIdentity(n int) Matrix` | MaxPlusIdentity returns the n-by-n max-plus identity matrix. |
| `func MaxPlusMatrix(data [][]float64) Matrix` | MaxPlusMatrix returns a max-plus Matrix holding a copy of data. |
| `func MaxPlusZeros(r, c int) Matrix` | MaxPlusZeros returns an r-by-c max-plus matrix of tropical zeros (-Inf). |
| `func MinPlusIdentity(n int) Matrix` | MinPlusIdentity returns the n-by-n min-plus identity matrix. |
| `func MinPlusMatrix(data [][]float64) Matrix` | MinPlusMatrix returns a min-plus Matrix holding a copy of data. |
| `func MinPlusZeros(r, c int) Matrix` | MinPlusZeros returns an r-by-c min-plus matrix of tropical zeros (+Inf). |
| `func NewMatrix(sr Semiring, data [][]float64) Matrix` | NewMatrix returns a Matrix over sr holding a deep copy of the rectangular slice data. |
| `func Zeros(sr Semiring, r, c int) Matrix` | Zeros returns an r-by-c matrix whose entries are all the tropical zero of the semiring. |
| `func (m Matrix) Add(n Matrix) (Matrix, error)` | Add returns the elementwise tropical sum of m and n. |
| `func (m Matrix) AssignmentValue(perm []int) float64` | AssignmentValue returns the tropical product of the entries m`i`[perm`i`] for a square matrix, the weight of the permutation as an assignment. |
| `func (m Matrix) At(i, j int) float64` | At returns the entry in row i and column j. |
| `func (m Matrix) Clone() Matrix` | Clone returns an independent deep copy of the matrix. |
| `func (m Matrix) Cofactor(i, j int) (float64, error)` | Cofactor returns the tropical cofactor of entry (i,j): the tropical permanent of the minor obtained by deleting row i and column j. |
| `func (m Matrix) Col(j int) []float64` | Col returns a fresh copy of column j. |
| `func (m Matrix) ColVector(j int) Vector` | ColVector returns column j as a Vector. |
| `func (m Matrix) Cols() int` | Cols returns the number of columns. |
| `func (m Matrix) CountFinite() int` | CountFinite returns the number of finite entries. |
| `func (m Matrix) CriticalNodes() ([]int, bool, error)` | CriticalNodes returns the indices of the nodes that lie on a cycle attaining the tropical eigenvalue (the critical cycles). |
| `func (m Matrix) Determinant() (float64, error)` | Determinant returns the tropical determinant of a square matrix. |
| `func (m Matrix) Diagonal() []float64` | Diagonal returns a fresh copy of the main diagonal. |
| `func (m Matrix) DirectSum(n Matrix) (Matrix, error)` | DirectSum returns the block-diagonal matrix with m in the top-left and n in the bottom-right, the off-diagonal blocks filled with tropical zeros. |
| `func (m Matrix) Dual() Matrix` | Dual returns the same matrix reinterpreted in the dual semiring, with every entry negated. |
| `func (m Matrix) Eigenvalue() (float64, bool, error)` | Eigenvalue returns the tropical eigenvalue of an irreducible square matrix: the maximum cycle mean for max-plus and the minimum cycle mean for… |
| `func (m Matrix) Eigenvector() (Vector, float64, bool, error)` | Eigenvector returns a tropical eigenvector v of an irreducible square matrix together with its eigenvalue lambda, so that A (*) v = lambda (*) v. |
| `func (m Matrix) ElementwiseMul(n Matrix) (Matrix, error)` | ElementwiseMul returns the elementwise (Hadamard) tropical product of m and n, whose entries are the ordinary sums of the corresponding entries. |
| `func (m Matrix) Equal(n Matrix) bool` | Equal reports whether m and n have the same semiring, shape and identical entries. |
| `func (m Matrix) EqualTol(n Matrix, tol float64) bool` | EqualTol reports whether m and n have the same semiring and shape and every pair of entries agrees to within tol (infinities must match exactly). |
| `func (m Matrix) GreatestSubSolution(b Vector) (Vector, error)` | GreatestSubSolution returns the greatest x with A (*) x <= b componentwise. |
| `func (m Matrix) HStack(n Matrix) (Matrix, error)` | HStack returns the horizontal concatenation [m \| n]. |
| `func (m Matrix) HasBadCycle() (bool, error)` | HasBadCycle reports whether the matrix has a cycle that makes its tropical closure diverge: a negative-weight cycle for min-plus or a positive-weight… |
| `func (m Matrix) IsDiagonal() bool` | IsDiagonal reports whether every off-diagonal entry is the tropical zero. |
| `func (m Matrix) IsFinite() bool` | IsFinite reports whether every entry is finite (no tropical zeros or other infinities). |
| `func (m Matrix) IsIdentity(tol float64) bool` | IsIdentity reports whether the matrix is square and equals the tropical identity to within tol. |
| `func (m Matrix) IsSquare() bool` | IsSquare reports whether the matrix has equally many rows and columns. |
| `func (m Matrix) IsSymmetric() bool` | IsSymmetric reports whether the matrix is square and equal to its transpose. |
| `func (m Matrix) IsTropicallySingular(tol float64) (bool, error)` | IsTropicallySingular reports whether a square matrix is tropically singular: its optimal assignment value is the tropical zero (no finite matching)… |
| `func (m Matrix) IsZero() bool` | IsZero reports whether every entry is the tropical zero. |
| `func (m Matrix) KroneckerProduct(n Matrix) (Matrix, error)` | KroneckerProduct returns the tropical Kronecker product of m and n, the (m.rows*n.rows)-by-(m.cols*n.cols) matrix whose block (i,j) is m`i``j`… |
| `func (m Matrix) LeftResidual(b Matrix) (Matrix, error)` | LeftResidual returns the greatest matrix X with A (*) X <= B, written A \ B. The (k,j) entry is the dual tropical sum over i of B`i``j` (/) A`i``k`. |
| `func (m Matrix) LongestPaths() (Matrix, error)` | LongestPaths returns the all-pairs longest-path matrix of a max-plus weighted adjacency matrix, treating the tropical zero (-Inf) as the absence of… |
| `func (m Matrix) Map(f func(float64) float64) Matrix` | Map returns a new matrix obtained by applying f to every entry. |
| `func (m Matrix) MaxCycleMean() (float64, bool, error)` | MaxCycleMean returns the maximum cycle mean of a square matrix interpreted as a weighted digraph (entry (i,j) is the weight of the edge i->j; the… |
| `func (m Matrix) MaxEntry() float64` | MaxEntry returns the largest entry in the usual numeric order, or -Inf for an empty matrix. |
| `func (m Matrix) MinCycleMean() (float64, bool, error)` | MinCycleMean returns the minimum cycle mean of a square matrix interpreted as a weighted digraph, computed with Karp's algorithm; it equals the… |
| `func (m Matrix) MinEntry() float64` | MinEntry returns the smallest entry in the usual numeric order, or +Inf for an empty matrix. |
| `func (m Matrix) Minor(p, q int) Matrix` | Minor returns the matrix with row p and column q removed. |
| `func (m Matrix) Mul(n Matrix) (Matrix, error)` | Mul returns the tropical matrix product m (*) n, where the (i,j) entry is the tropical sum over k of m`i``k` (*) n`k``j`. |
| `func (m Matrix) MulVec(v Vector) (Vector, error)` | MulVec returns the tropical matrix-vector product m (*) v. |
| `func (m Matrix) OptimalAssignment() (float64, []int, error)` | OptimalAssignment returns the tropical permanent value and an optimal permutation for a square matrix. |
| `func (m Matrix) PathWeight(nodes []int) float64` | PathWeight returns the tropical product of the edge weights along the walk visiting the given node sequence, that is the sum of… |
| `func (m Matrix) Permanent() (float64, error)` | Permanent returns the tropical permanent of a square matrix, the tropical sum over all permutations sigma of the tropical product of the entries A[i,… |
| `func (m Matrix) PermanentBrute() (float64, error)` | PermanentBrute returns the tropical permanent computed by explicitly enumerating every permutation. |
| `func (m Matrix) PlusClosure() (Matrix, error)` | PlusClosure returns the tropical Kleene plus A^+ = A (+) A^2 (+) A^3 (+) ..., the tropical sum of all positive powers of a square matrix, computed… |
| `func (m Matrix) Pow(n int) (Matrix, error)` | Pow returns the tropical matrix power m^n by repeated squaring. |
| `func (m Matrix) PrincipalSolution(b Vector) (Vector, error)` | PrincipalSolution returns the principal solution of the tropical linear system A (*) x = b obtained by residuation. |
| `func (m Matrix) Raw() [][]float64` | Raw returns a fresh deep copy of the entries as a [][]float64. |
| `func (m Matrix) Reachability() ([][]bool, error)` | Reachability returns a boolean matrix whose (i,j) entry reports whether j is reachable from i in at least one step, i.e. |
| `func (m Matrix) RightResidual(b Matrix) (Matrix, error)` | RightResidual returns the greatest matrix X with X (*) A <= B, written B / A. The (i,k) entry is the dual tropical sum over j of B`i``j` (/) A`k``j`. |
| `func (m Matrix) Row(i int) []float64` | Row returns a fresh copy of row i. |
| `func (m Matrix) RowVector(i int) Vector` | RowVector returns row i as a Vector. |
| `func (m Matrix) Rows() int` | Rows returns the number of rows. |
| `func (m Matrix) ScalarAdd(c float64) Matrix` | ScalarAdd returns the matrix obtained by tropically adding c to every entry. |
| `func (m Matrix) ScalarMul(c float64) Matrix` | ScalarMul returns the matrix obtained by tropically multiplying every entry by c. |
| `func (m Matrix) Semiring() Semiring` | Semiring returns the semiring under which the matrix is interpreted. |
| `func (m Matrix) Set(i, j int, value float64)` | Set stores value in row i and column j. |
| `func (m Matrix) ShortestPaths() (Matrix, error)` | ShortestPaths returns the all-pairs shortest-path matrix of a min-plus weighted adjacency matrix, treating the tropical zero (+Inf) as the absence of… |
| `func (m Matrix) SingleSourceShortestPaths(src int) (Vector, error)` | SingleSourceShortestPaths returns the vector of best walk weights from src to every node, computed with the Bellman-Ford relaxation: shortest paths… |
| `func (m Matrix) SolveAffine(b Vector) (Vector, error)` | SolveAffine returns the least solution of the affine tropical fixed-point equation x = A (*) x (+) b, which equals A* (*) b where A* is the Kleene… |
| `func (m Matrix) SolveExact(b Vector, tol float64) (Vector, bool, error)` | SolveExact returns the principal solution x of A (*) x = b together with a boolean reporting whether it solves the system exactly to within tol. |
| `func (m Matrix) Star() (Matrix, error)` | Star returns the tropical Kleene star A* = I (+) A (+) A^2 (+) ..., the reflexive-transitive closure of a square matrix. |
| `func (m Matrix) StarSeries(k int) (Matrix, error)` | StarSeries returns the bounded tropical closure I (+) A (+) A^2 (+) ... |
| `func (m Matrix) String() string` | String renders the matrix with one row per line and space-separated entries. |
| `func (m Matrix) SubMatrix(r0, r1, c0, c1 int) Matrix` | SubMatrix returns the submatrix spanning rows [r0, r1) and columns [c0, c1). |
| `func (m Matrix) Trace() (float64, error)` | Trace returns the tropical trace, the tropical sum of the diagonal entries. |
| `func (m Matrix) Transpose() Matrix` | Transpose returns the transpose of the matrix. |
| `func (m Matrix) VStack(n Matrix) (Matrix, error)` | VStack returns the vertical concatenation of m stacked above n. |
| `func (m Matrix) VecMul(v Vector) (Vector, error)` | VecMul returns the tropical row-vector-matrix product v^T (*) m. |

</details>

<details>
<summary><code>Poly</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FromRoots(sr Semiring, roots []float64) Poly` | FromRoots returns the monic tropical polynomial over sr whose tropical roots are exactly the given values, that is the tropical product of the linear… |
| `func MaxPlusFromRoots(roots []float64) Poly` | MaxPlusFromRoots returns the monic max-plus polynomial with the given roots. |
| `func MaxPlusPoly(coeffs []float64) Poly` | MaxPlusPoly returns a max-plus tropical polynomial with the given coefficients. |
| `func MinPlusFromRoots(roots []float64) Poly` | MinPlusFromRoots returns the monic min-plus polynomial with the given roots. |
| `func MinPlusPoly(coeffs []float64) Poly` | MinPlusPoly returns a min-plus tropical polynomial with the given coefficients. |
| `func NewPoly(sr Semiring, coeffs []float64) Poly` | NewPoly returns the tropical polynomial over sr with the given coefficients, coeffs`i` belonging to x^i. |
| `func (p Poly) Add(q Poly) Poly` | Add returns the tropical sum of p and q, computed coefficientwise. |
| `func (p Poly) Clone() Poly` | Clone returns an independent copy of the polynomial. |
| `func (p Poly) Coeff(i int) float64` | Coeff returns the coefficient of x^i, or the tropical zero when i is out of range. |
| `func (p Poly) Coeffs() []float64` | Coeffs returns a fresh copy of the coefficient slice, index i holding the coefficient of x^i. |
| `func (p Poly) ConstantCoeff() float64` | ConstantCoeff returns the coefficient of x^0. |
| `func (p Poly) Degree() int` | Degree returns the degree of the polynomial: the highest power with a non-zero coefficient. |
| `func (p Poly) Dual() Poly` | Dual returns the polynomial reinterpreted in the dual semiring with every coefficient negated; its roots are the negatives of the original roots. |
| `func (p Poly) Equal(q Poly) bool` | Equal reports whether p and q are over the same semiring and have identical coefficients. |
| `func (p Poly) Eval(x float64) float64` | Eval returns the value of the tropical polynomial at x. |
| `func (p Poly) EvalAll(xs []float64) []float64` | EvalAll evaluates the polynomial at every x in xs and returns the values. |
| `func (p Poly) EvalArg(x float64) (float64, int)` | EvalArg returns the value of the polynomial at x together with the index of a term attaining it (the active monomial). |
| `func (p Poly) IsMonic() bool` | IsMonic reports whether the leading coefficient equals the tropical one (0). |
| `func (p Poly) IsRoot(x, tol float64) bool` | IsRoot reports whether x is a tropical root of the polynomial, i.e. |
| `func (p Poly) LeadingCoeff() float64` | LeadingCoeff returns the coefficient of the highest power of the polynomial. |
| `func (p Poly) Mul(q Poly) Poly` | Mul returns the tropical product of p and q, computed as a tropical convolution of their coefficients. |
| `func (p Poly) NewtonPolygon() []PolygonVertex` | NewtonPolygon returns the vertices of the Newton polygon of the polynomial: the lower convex hull of the finite coefficient points (i, coeff`i`) for… |
| `func (p Poly) Pow(n int) Poly` | Pow returns the tropical power p^n. |
| `func (p Poly) RootValues() []float64` | RootValues returns the distinct tropical roots repeated according to their multiplicity, ordered by increasing value. |
| `func (p Poly) Roots() []Root` | Roots returns the tropical roots of the polynomial with their multiplicities, read off the Newton polygon and ordered by increasing root value. |
| `func (p Poly) ScalarMul(c float64) Poly` | ScalarMul returns the polynomial with every coefficient tropically multiplied by c. |
| `func (p Poly) Semiring() Semiring` | Semiring returns the semiring under which the polynomial is interpreted. |
| `func (p Poly) ShiftUp(k int) Poly` | ShiftUp returns the polynomial multiplied by x^k, prepending k tropical-zero coefficients. |
| `func (p Poly) String() string` | String renders the polynomial in the form "a0 (+) a1 x (+) a2 x^2 (+) ...", omitting terms with a tropical-zero coefficient. |
| `func (p Poly) TropicallyVanishes(x, tol float64) bool` | TropicallyVanishes reports whether the tropical polynomial vanishes at x, i.e. |

</details>

<details>
<summary><code>Semiring</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MaxPlusSemiring() Semiring` | MaxPlusSemiring returns the max-plus semiring. |
| `func MinPlusSemiring() Semiring` | MinPlusSemiring returns the min-plus semiring. |
| `func NewSemiring(kind Kind) Semiring` | NewSemiring returns the Semiring for the given kind. |
| `func (s Semiring) Add(a, b float64) float64` | Add returns the tropical sum a (+) b: the minimum for min-plus and the maximum for max-plus. |
| `func (s Semiring) AtLeastAsGood(a, b float64) bool` | AtLeastAsGood reports whether a is at least as good as b in the sense of tropical addition, i.e. |
| `func (s Semiring) Better(a, b float64) float64` | Better returns whichever of a and b tropical addition selects, identical to Add but named for readability at call sites that compare candidates. |
| `func (s Semiring) Div(a, b float64) float64` | Div returns the tropical quotient a (/) b, the residual of multiplication: on finite operands it is the ordinary difference a-b. |
| `func (s Semiring) Dual() Semiring` | Dual returns the order-dual semiring: MinPlus for MaxPlus and vice versa. |
| `func (s Semiring) EqualScalar(a, b, tol float64) bool` | EqualScalar reports whether a and b agree to within tol, with infinities required to match exactly and in sign. |
| `func (s Semiring) FormatScalar(a float64) string` | FormatScalar renders a tropical scalar, printing the tropical zero as the mathematical symbol for infinity with the appropriate sign. |
| `func (s Semiring) IsMaxPlus() bool` | IsMaxPlus reports whether s is the max-plus semiring. |
| `func (s Semiring) IsMinPlus() bool` | IsMinPlus reports whether s is the min-plus semiring. |
| `func (s Semiring) IsOne(a float64) bool` | IsOne reports whether a equals the tropical one (0). |
| `func (s Semiring) IsZero(a float64) bool` | IsZero reports whether a is the tropical zero of the semiring. |
| `func (s Semiring) Kind() Kind` | Kind returns the kind of the semiring. |
| `func (s Semiring) Mul(a, b float64) float64` | Mul returns the tropical product a (*) b. |
| `func (s Semiring) One() float64` | One returns the tropical one (multiplicative identity) of the semiring, which is 0 for both min-plus and max-plus. |
| `func (s Semiring) Pow(a float64, n int) float64` | Pow returns the tropical power a raised to the integer exponent n, which equals n*a for finite a. |
| `func (s Semiring) Prod(xs ...float64) float64` | Prod returns the tropical product of all arguments, or the tropical one when no arguments are given. |
| `func (s Semiring) Star(a float64) float64` | Star returns the tropical Kleene star a* = 1 (+) a (+) a^2 (+) ... |
| `func (s Semiring) StarConverges(a float64) bool` | StarConverges reports whether the scalar Kleene star a* is finite: a >= 0 for min-plus and a <= 0 for max-plus. |
| `func (s Semiring) String() string` | String returns the human-readable name of the semiring. |
| `func (s Semiring) Sum(xs ...float64) float64` | Sum returns the tropical sum of all arguments, or the tropical zero when no arguments are given. |
| `func (s Semiring) Zero() float64` | Zero returns the tropical zero (additive identity) of the semiring: +Inf for min-plus and -Inf for max-plus. |

</details>

<details>
<summary><code>Vector</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MaxPlusVector(data []float64) Vector` | MaxPlusVector returns a max-plus Vector holding a copy of data. |
| `func MinPlusVector(data []float64) Vector` | MinPlusVector returns a min-plus Vector holding a copy of data. |
| `func NewVector(sr Semiring, data []float64) Vector` | NewVector returns a Vector over the given semiring holding a copy of data. |
| `func OneVector(sr Semiring, n int) Vector` | OneVector returns a length-n Vector whose entries are all the tropical one (0) of the semiring. |
| `func UnitVector(sr Semiring, n, i int) Vector` | UnitVector returns the length-n tropical unit vector: the tropical one (0) at index i and the tropical zero elsewhere. |
| `func ZeroVector(sr Semiring, n int) Vector` | ZeroVector returns a length-n Vector whose entries are all the tropical zero of the semiring. |
| `func (v Vector) Add(w Vector) Vector` | Add returns the elementwise tropical sum of v and w. |
| `func (v Vector) ArgBest() (int, float64)` | ArgBest returns the index of the entry selected by tropical addition (the index of the minimum for min-plus, of the maximum for max-plus) and its… |
| `func (v Vector) Argmax() int` | Argmax returns the index of the numerically largest entry, or -1 for an empty vector. |
| `func (v Vector) Argmin() int` | Argmin returns the index of the numerically smallest entry, or -1 for an empty vector. |
| `func (v Vector) At(i int) float64` | At returns the i-th entry. |
| `func (v Vector) Clone() Vector` | Clone returns an independent copy of the vector. |
| `func (v Vector) Dot(w Vector) float64` | Dot returns the tropical inner product of v and w: the tropical sum over i of v`i` (*) w`i`. |
| `func (v Vector) Dual() Vector` | Dual returns the same vector reinterpreted in the dual semiring, with every entry negated. |
| `func (v Vector) Equal(w Vector) bool` | Equal reports whether v and w have the same semiring, length and identical entries. |
| `func (v Vector) EqualTol(w Vector, tol float64) bool` | EqualTol reports whether v and w have the same semiring and length and every pair of entries agrees to within tol (with infinities required to match… |
| `func (v Vector) IsZero() bool` | IsZero reports whether every entry equals the tropical zero. |
| `func (v Vector) Len() int` | Len returns the number of entries in the vector. |
| `func (v Vector) Map(f func(float64) float64) Vector` | Map returns a new vector obtained by applying f to every entry. |
| `func (v Vector) Max() float64` | Max returns the largest entry in the usual numeric order. |
| `func (v Vector) Min() float64` | Min returns the smallest entry, ignoring sign of infinities in the usual numeric order. |
| `func (v Vector) Normalize() Vector` | Normalize returns the vector with its tropical-product "length" shifted so that the best entry becomes the tropical one (0): every entry has the best… |
| `func (v Vector) Prod() float64` | Prod returns the tropical product of all entries (their ordinary sum). |
| `func (v Vector) ScalarAdd(c float64) Vector` | ScalarAdd returns the vector obtained by tropically adding c to every entry. |
| `func (v Vector) ScalarMul(c float64) Vector` | ScalarMul returns the vector obtained by tropically multiplying every entry by c. |
| `func (v Vector) Semiring() Semiring` | Semiring returns the semiring under which the vector is interpreted. |
| `func (v Vector) Set(i int, value float64)` | Set stores value at index i. |
| `func (v Vector) Slice() []float64` | Slice returns a fresh copy of the underlying entries. |
| `func (v Vector) String() string` | String renders the vector as a bracketed, space-separated list. |
| `func (v Vector) Sum() float64` | Sum returns the tropical sum of all entries (the minimum for min-plus, the maximum for max-plus). |

</details>

### Variables

`ErrDim`, `ErrDivergent`, `ErrNotSquare`, `ErrTooLarge`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/tropical).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
