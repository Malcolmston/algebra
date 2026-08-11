# grouprep — implements the ordinary representation theory of finite groups over the complex

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/grouprep.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/grouprep)

Package grouprep implements the ordinary representation theory of finite
groups over the complex numbers, built entirely on the Go standard library.

A (matrix) representation of a finite group G is a homomorphism ρ: G →
GL(d, ℂ) assigning to every group element an invertible d×d complex matrix
such that ρ(gh) = ρ(g)ρ(h). This package provides the concrete data
structures and algorithms of the classical theory: complex matrices
(`Matrix`), finite groups given by their Cayley table (`Group`), matrix
representations (`Rep`), characters (`Character`) and character tables
(`CharacterTable`).

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/grouprep"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AbsC(z complex128) float64` | AbsC returns the modulus \|z\| of the complex number z. |
| `func ApproxEqualC(a, b complex128, tol float64) bool` | ApproxEqualC reports whether a and b are within tol of each other in modulus, i.e. |
| `func ApproxZeroC(z complex128, tol float64) bool` | ApproxZeroC reports whether \|z\| <= tol. |
| `func CharacterDimension(chi Character) int` | CharacterDimension returns the dimension of the representation affording chi, the rounded real part of chi on the identity. |
| `func ClassValues(g *Group, chi Character, tol float64) ([]complex128, error)` | ClassValues reduces chi to one value per conjugacy class of g, taking the value on each class representative (the smallest element). |
| `func Conj(z complex128) complex128` | Conj returns the complex conjugate of z. |
| `func Cplx(re, im float64) complex128` | Cplx builds a complex128 from its real and imaginary parts. |
| `func DecomposeCharacter(g *Group, chi Character, irrs []Character) []int` | DecomposeCharacter returns the multiplicities of chi against the given list of irreducible characters, so that chi = Σ mult`i`·irrs`i`. |
| `func Imag(z complex128) float64` | Imag returns the imaginary part of z. |
| `func InnerProduct(g *Group, a, b Character) complex128` | InnerProduct returns the Hermitian inner product of two characters of g, |
| `func IsCharacter(g *Group, chi Character, irrs []Character, tol float64) bool` | IsCharacter reports whether chi is a (non-virtual) character with respect to the given complete list of irreducibles: all multiplicities are… |
| `func IsIntegerC(z complex128, tol float64) bool` | IsIntegerC reports whether z is a (real) integer to within tol. |
| `func IsIrreducible(g *Group, chi Character, tol float64) bool` | IsIrreducible reports whether chi is an irreducible character, i.e. |
| `func IsRealC(z complex128, tol float64) bool` | IsRealC reports whether z is real to within tol, i.e. |
| `func Multiplicity(g *Group, chi, irr Character) int` | Multiplicity returns the multiplicity of the irreducible character irr in the character chi, the rounded value of 〈irr, chi〉. |
| `func NormSquared(g *Group, chi Character) float64` | NormSquared returns the real number 〈chi, chi〉, the squared norm of chi in the character inner product. |
| `func PrimitiveRootOfUnity(n int) complex128` | PrimitiveRootOfUnity returns exp(2πi / n), the standard primitive n-th root of unity. |
| `func Real(z complex128) float64` | Real returns the real part of z. |
| `func RootOfUnity(n, k int) complex128` | RootOfUnity returns the primitive-n-th root exp(2πi k / n), the k-th power of the standard primitive n-th root of unity. |
| `func RoundC(z complex128, decimals int) complex128` | RoundC rounds both the real and imaginary parts of z to the given number of decimal places. |
| `func RoundToInteger(z complex128) complex128` | RoundToInteger returns the nearest Gaussian integer to z. |

### Types

| Type | What it is |
| --- | --- |
| `Character` | Character is a class function on a group, stored as one complex value per element index: chi`i` is the value on element i. |
| `CharacterTable` | CharacterTable holds a complete set of irreducible characters of a group, tabulated against its conjugacy classes. |
| `Group` | Group is a finite group given by an explicit Cayley (multiplication) table. |
| `Matrix` | Matrix is a dense complex matrix stored row-major: m`i``j` is the entry in row i and column j. |
| `Perm` | Perm is a permutation of {0,...,n-1} in one-line notation: p`i` is the image of i. |
| `Rep` | Rep is a finite-dimensional complex matrix representation of a `Group`: a map assigning to every element i an invertible dim×dim matrix mats`i` such… |

<details>
<summary><code>Character</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AddCharacters(a, b Character) Character` | AddCharacters returns the pointwise sum a+b, the character of the direct sum. |
| `func CharacterFromRep(r *Rep) Character` | CharacterFromRep returns the character of the representation r; it is a free function equivalent to r.Character(). |
| `func ConjugateCharacter(chi Character) Character` | ConjugateCharacter returns the complex conjugate of chi, the character of the dual representation. |
| `func InduceCharacter(gG *Group, gH *Group, emb []int, psi Character) (Character, error)` | InduceCharacter induces a character psi of a subgroup H up to the group G, using the embedding emb (emb`h` the index in G of the H-element h). |
| `func RegularCharacter(g *Group) Character` | RegularCharacter returns the character of the regular representation of g: \|G\| on the identity and 0 elsewhere. |
| `func RestrictCharacter(chi Character, emb []int) Character` | RestrictCharacter restricts a character of a group G to a subgroup H given by an embedding emb, where emb`h` is the index in G of the element h of H.… |
| `func ScaleCharacter(z complex128, chi Character) Character` | ScaleCharacter returns z·chi. |
| `func SubCharacters(a, b Character) Character` | SubCharacters returns the pointwise difference a-b (a virtual character). |
| `func TensorCharacters(a, b Character) Character` | TensorCharacters returns the pointwise product a·b, the character of the tensor product representation. |
| `func TrivialCharacter(g *Group) Character` | TrivialCharacter returns the character of the trivial representation of g, identically 1. |
| `func ZeroCharacter(g *Group) Character` | ZeroCharacter returns the identically zero class function on g. |
| `func (chi Character) ApproxEqual(psi Character, tol float64) bool` | ApproxEqual reports whether chi and psi agree entrywise to within tol. |
| `func (chi Character) Clone() Character` | Clone returns an independent copy of chi. |
| `func (chi Character) Degree() float64` | Degree returns the degree (dimension) of the character, its value on the identity element, as a real number. |
| `func (chi Character) IsClassFunction(g *Group, tol float64) bool` | IsClassFunction reports whether chi is constant on the conjugacy classes of g to within tol. |
| `func (chi Character) Value(i int) complex128` | Value returns the value of chi on element i. |

</details>

<details>
<summary><code>CharacterTable</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BuildCharacterTable(g *Group, irrs []Character) *CharacterTable` | BuildCharacterTable assembles a character table for g from a list of irreducible characters (each a full class function over all \|G\| elements). |
| `func CharacterTableCyclic(n int) *CharacterTable` | CharacterTableCyclic returns the character table of the cyclic group C_n. |
| `func CharacterTableDihedral(n int) *CharacterTable` | CharacterTableDihedral returns the character table of the dihedral group D_n. |
| `func CharacterTableKleinFour() *CharacterTable` | CharacterTableKleinFour returns the character table of the Klein four-group. |
| `func CharacterTableQuaternion() *CharacterTable` | CharacterTableQuaternion returns the character table of the quaternion group Q8: four one-dimensional characters and one two-dimensional character… |
| `func CharacterTableSymmetric(n int) *CharacterTable` | CharacterTableSymmetric returns the character table of the symmetric group S_n for 1 <= n <= 4. |
| `func (ct *CharacterTable) ClassElementOrders() []int` | ClassElementOrders returns the order of a representative element of each conjugacy class, in column order. |
| `func (ct *CharacterTable) ClassRepresentatives() []int` | ClassRepresentatives returns the smallest element index of each conjugacy class, in column order. |
| `func (ct *CharacterTable) ClassSizes() []int` | ClassSizes returns the number of elements in each conjugacy class, in column order. |
| `func (ct *CharacterTable) ColumnOrthogonality(tol float64) bool` | ColumnOrthogonality reports whether the columns satisfy the second orthogonality relation: for classes c and d with representatives of centralizer… |
| `func (ct *CharacterTable) Decompose(chi Character) []int` | Decompose returns the multiplicities of an arbitrary character chi against the irreducibles of the table, so that chi = Σ mult`i`·χᵢ. |
| `func (ct *CharacterTable) Dimensions() []int` | Dimensions returns the degrees of the irreducible characters, in row order. |
| `func (ct *CharacterTable) Group() *Group` | Group returns the group of the table. |
| `func (ct *CharacterTable) Irreducible(i int) Character` | Irreducible returns the i-th irreducible character as a full class function over the group. |
| `func (ct *CharacterTable) Irreducibles() []Character` | Irreducibles returns all irreducible characters as full class functions. |
| `func (ct *CharacterTable) IsComplete() bool` | IsComplete reports whether the table has as many irreducibles as classes and satisfies Σ dᵢ² = \|G\|. |
| `func (ct *CharacterTable) NumClasses() int` | NumClasses returns the number of conjugacy classes (columns). |
| `func (ct *CharacterTable) NumIrreducibles() int` | NumIrreducibles returns the number of irreducible characters (rows). |
| `func (ct *CharacterTable) RowOrthogonality(tol float64) bool` | RowOrthogonality reports whether the rows satisfy the first orthogonality relation: 〈χᵢ, χⱼ〉 = δᵢⱼ to within tol, using the weighted… |
| `func (ct *CharacterTable) String() string` | String renders the character table as an aligned text grid: one row per irreducible, one column per class, values rounded to two decimals. |
| `func (ct *CharacterTable) SumOfSquaresOfDims() int` | SumOfSquaresOfDims returns Σ dᵢ², which for a complete table equals the group order \|G\|. |
| `func (ct *CharacterTable) Value(i, c int) complex128` | Value returns the value of irreducible i on conjugacy class c. |

</details>

<details>
<summary><code>Group</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AlternatingGroup(n int) *Group` | AlternatingGroup returns the alternating group A_n of even permutations of n points, of order n!/2 for n >= 2. |
| `func CyclicGroup(n int) *Group` | CyclicGroup returns the cyclic group Z/nZ of order n written additively: element i is i mod n and i·j = (i+j) mod n. |
| `func DihedralGroup(n int) *Group` | DihedralGroup returns the dihedral group D_n of order 2n, the symmetry group of a regular n-gon. |
| `func DirectProduct(g, h *Group) *Group` | DirectProduct returns the external direct product g×h. |
| `func KleinFourGroup() *Group` | KleinFourGroup returns the Klein four-group V = C2×C2 of order 4. |
| `func QuaternionGroup() *Group` | QuaternionGroup returns the quaternion group Q8 = {±1, ±i, ±j, ±k} of order 8. |
| `func SymmetricGroup(n int) *Group` | SymmetricGroup returns the symmetric group S_n of all permutations of n points, of order n!. |
| `func TrivialGroup() *Group` | TrivialGroup returns the group with one element. |
| `func (g *Group) CayleyTable() [][]int` | CayleyTable returns a copy of the multiplication table of g. |
| `func (g *Group) Center() []int` | Center returns the sorted indices of the elements commuting with every element of g. |
| `func (g *Group) Centralizer(x int) []int` | Centralizer returns the sorted indices of the elements commuting with x. |
| `func (g *Group) ClassOf(x int) []int` | ClassOf returns the conjugacy class containing element x. |
| `func (g *Group) Commutator(i, j int) int` | Commutator returns the group commutator i j i⁻¹ j⁻¹. |
| `func (g *Group) ConjugacyClasses() [][]int` | ConjugacyClasses returns the conjugacy classes of g, each a sorted slice of element indices. |
| `func (g *Group) Conjugate(a, x int) int` | Conjugate returns the conjugate g x g⁻¹ (with a the conjugator and x the element), i.e. |
| `func (g *Group) ElementName(i int) string` | ElementName returns a human-readable label for element i. |
| `func (g *Group) ElementOrder(i int) int` | ElementOrder returns the order of element i, the least positive k with i^k = identity. |
| `func (g *Group) Elements() []int` | Elements returns the slice {0,1,...,Order()-1}. |
| `func (g *Group) Exponent() int` | Exponent returns the exponent of g, the least common multiple of all element orders. |
| `func (g *Group) GeneratedBy(gens []int) []int` | GeneratedBy returns the sorted set of elements of the subgroup generated by gens (the closure under multiplication and inverses). |
| `func (g *Group) Identity() int` | Identity returns the index of the identity element, which is always 0. |
| `func (g *Group) Inverse(i int) int` | Inverse returns the inverse of element i. |
| `func (g *Group) IsAbelian() bool` | IsAbelian reports whether g is commutative. |
| `func (g *Group) IsSubgroup(elems []int) bool` | IsSubgroup reports whether the given set of element indices is closed under the group operation and inverses (and hence forms a subgroup). |
| `func (g *Group) IsValid() error` | IsValid re-checks that g really is a group. |
| `func (g *Group) Mul(i, j int) int` | Mul returns the product i·j. |
| `func (g *Group) Name() string` | Name returns the display name of g. |
| `func (g *Group) NumConjugacyClasses() int` | NumConjugacyClasses returns the number of conjugacy classes of g, which equals the number of complex irreducible representations. |
| `func (g *Group) Order() int` | Order returns the number of elements of g. |
| `func (g *Group) Pow(i, k int) int` | Pow returns element i raised to the integer power k (negative k uses the inverse). |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Commutator(a, b Matrix) Matrix` | Commutator returns the matrix commutator ab-ba. |
| `func DiagMatrix(diag []complex128) Matrix` | DiagMatrix returns the square diagonal matrix whose diagonal entries are the elements of diag. |
| `func IdentityMatrix(n int) Matrix` | IdentityMatrix returns the n×n identity matrix. |
| `func MatAdd(a, b Matrix) (Matrix, error)` | MatAdd returns a+b, erroring on a shape mismatch. |
| `func MatDirectSum(a, b Matrix) Matrix` | MatDirectSum returns the block-diagonal direct sum a⊕b. |
| `func MatKronecker(a, b Matrix) Matrix` | MatKronecker returns the Kronecker product a⊗b. |
| `func MatMul(a, b Matrix) (Matrix, error)` | MatMul returns the product a·b. |
| `func MatSub(a, b Matrix) (Matrix, error)` | MatSub returns a-b, erroring on a shape mismatch. |
| `func MatrixFromReal(rows [][]float64) Matrix` | MatrixFromReal builds a complex matrix whose entries are the given real numbers (with zero imaginary part). |
| `func MatrixFromRows(rows [][]complex128) Matrix` | MatrixFromRows builds a matrix from an explicit slice of rows. |
| `func NewMatrix(r, c int) Matrix` | NewMatrix returns an r×c matrix with every entry equal to zero. |
| `func PermutationMatrix(p []int) Matrix` | PermutationMatrix returns the n×n permutation matrix P of the permutation p (in one-line notation, p`i` the image of i), defined by P[p`i`]`i` = 1. |
| `func ReflectionMatrix(theta float64) Matrix` | ReflectionMatrix returns the 2×2 real matrix reflecting the plane across the line through the origin at angle theta to the x-axis. |
| `func RotationMatrix(theta float64) Matrix` | RotationMatrix returns the 2×2 real rotation matrix through angle theta (counter-clockwise). |
| `func ScalarMatrix(n int, z complex128) Matrix` | ScalarMatrix returns the n×n diagonal matrix with z on the diagonal, i.e. |
| `func ZeroMatrix(r, c int) Matrix` | ZeroMatrix is an alias for `NewMatrix` returning the r×c zero matrix. |
| `func (m Matrix) Add(n Matrix) (Matrix, error)` | Add returns m+n. |
| `func (m Matrix) ApproxEqual(n Matrix, tol float64) bool` | ApproxEqual reports whether m and n have the same shape and entries agreeing to within tol in modulus. |
| `func (m Matrix) At(i, j int) complex128` | At returns the entry m`i``j`. |
| `func (m Matrix) Clone() Matrix` | Clone returns an independent deep copy of m. |
| `func (m Matrix) Col(j int) []complex128` | Col returns a copy of column j. |
| `func (m Matrix) Cols() int` | Cols returns the number of columns of m (zero for an empty matrix). |
| `func (m Matrix) Conjugate() Matrix` | Conjugate returns the entrywise complex conjugate of m. |
| `func (m Matrix) ConjugateTranspose() Matrix` | ConjugateTranspose returns the conjugate transpose (Hermitian adjoint) m†. |
| `func (m Matrix) Det() complex128` | Det returns the determinant of the square matrix m, computed by Gaussian elimination with partial pivoting. |
| `func (m Matrix) DirectSum(n Matrix) Matrix` | DirectSum returns the block-diagonal matrix diag(m, n), the direct sum m⊕n. |
| `func (m Matrix) Equal(n Matrix) bool` | Equal reports whether m and n have the same shape and identical entries. |
| `func (m Matrix) FrobeniusNorm() float64` | FrobeniusNorm returns the Frobenius norm of m, the square root of the sum of the squared moduli of all entries. |
| `func (m Matrix) Inverse() (Matrix, error)` | Inverse returns the inverse of the square matrix m via Gauss-Jordan elimination. |
| `func (m Matrix) IsDiagonal(tol float64) bool` | IsDiagonal reports whether all off-diagonal entries of the square matrix m vanish to within tol. |
| `func (m Matrix) IsHermitian(tol float64) bool` | IsHermitian reports whether m equals its conjugate transpose to within tol. |
| `func (m Matrix) IsIdentity(tol float64) bool` | IsIdentity reports whether m is the identity matrix to within tol. |
| `func (m Matrix) IsSquare() bool` | IsSquare reports whether m has equal row and column counts. |
| `func (m Matrix) IsUnitary(tol float64) bool` | IsUnitary reports whether m is unitary, i.e. |
| `func (m Matrix) IsZero(tol float64) bool` | IsZero reports whether every entry of m is zero to within tol. |
| `func (m Matrix) Kronecker(n Matrix) Matrix` | Kronecker returns the Kronecker (tensor) product m⊗n, an (r_m·r_n)×(c_m·c_n) matrix. |
| `func (m Matrix) MaxAbsDiff(n Matrix) float64` | MaxAbsDiff returns the largest modulus \|m`i``j`-n`i``j`\| over all entries. |
| `func (m Matrix) Mul(n Matrix) (Matrix, error)` | Mul returns the matrix product m·n. |
| `func (m Matrix) MulVec(v []complex128) ([]complex128, error)` | MulVec returns the matrix-vector product m·v. |
| `func (m Matrix) Pow(k int) (Matrix, error)` | Pow returns m raised to the non-negative integer power k. |
| `func (m Matrix) Rank(tol float64) int` | Rank returns the numerical rank of m, the number of pivots found by Gaussian elimination with entries below tol treated as zero. |
| `func (m Matrix) Round(decimals int) Matrix` | Round returns a copy of m with every entry rounded to the given number of decimal places. |
| `func (m Matrix) Row(i int) []complex128` | Row returns a copy of row i. |
| `func (m Matrix) Rows() int` | Rows returns the number of rows of m. |
| `func (m Matrix) Scale(z complex128) Matrix` | Scale returns z·m, the matrix with every entry multiplied by z. |
| `func (m Matrix) Set(i, j int, z complex128)` | Set assigns z to entry m`i``j`. |
| `func (m Matrix) String() string` | String renders m with entries rounded to three decimals, one row per line. |
| `func (m Matrix) Sub(n Matrix) (Matrix, error)` | Sub returns m-n. |
| `func (m Matrix) Submatrix(rows, cols []int) Matrix` | Submatrix returns the matrix formed by selecting the given rows and columns in the given order. |
| `func (m Matrix) Trace() complex128` | Trace returns the sum of the diagonal entries of the square matrix m. |
| `func (m Matrix) Transpose() Matrix` | Transpose returns the transpose of m. |

</details>

<details>
<summary><code>Perm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AllPermutations(n int) []Perm` | AllPermutations returns every permutation of {0,...,n-1} in lexicographic order. |
| `func IdentityPerm(n int) Perm` | IdentityPerm returns the identity permutation of degree n. |
| `func (p Perm) Apply(i int) int` | Apply returns the image p`i` of point i. |
| `func (p Perm) Compose(q Perm) Perm` | Compose returns the composition p∘q, applying q first: (p∘q)`i` = p[q`i`]. |
| `func (p Perm) CycleDecomposition() [][]int` | CycleDecomposition returns the disjoint cycles of p as a slice of slices, including fixed points as length-one cycles. |
| `func (p Perm) CycleType() []int` | CycleType returns the multiset of cycle lengths of p in non-increasing order, a partition of the degree. |
| `func (p Perm) Degree() int` | Degree returns the size n of the set {0,...,n-1} that p permutes. |
| `func (p Perm) Equal(q Perm) bool` | Equal reports whether p and q are the same permutation. |
| `func (p Perm) Inverse() Perm` | Inverse returns the inverse permutation of p. |
| `func (p Perm) IsValid() bool` | IsValid reports whether p is a genuine permutation of {0,...,len(p)-1}. |
| `func (p Perm) Order() int` | Order returns the multiplicative order of p, the least positive k with p^k = identity. |
| `func (p Perm) PermMatrix() Matrix` | PermMatrix returns the permutation matrix of p; see `PermutationMatrix`. |
| `func (p Perm) Sign() int` | Sign returns the sign of p: +1 if p is an even permutation and -1 if odd, computed from the parity of the number of transpositions. |

</details>

<details>
<summary><code>Rep</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CyclicRep(n, k int) *Rep` | CyclicRep returns the one-dimensional representation of the cyclic group C_n sending the generator r to the k-th root of unity ω^k = exp(2πik/n). |
| `func DihedralRep2D(n, h int) *Rep` | DihedralRep2D returns the two-dimensional representation of the dihedral group D_n in which the rotation r acts by rotation through 2πh/n and the… |
| `func DirectSumReps(reps ...*Rep) *Rep` | DirectSumReps returns the direct sum of several representations of the same group, in order. |
| `func NaturalRepSymmetric(n int) *Rep` | NaturalRepSymmetric returns the natural n-dimensional permutation representation of the symmetric group S_n, in which each permutation acts by its… |
| `func NewRep(g *Group, mats []Matrix, tol float64) (*Rep, error)` | NewRep builds a representation of g from an explicit matrix per element and verifies that it is a genuine homomorphism to within tol. |
| `func QuaternionRep() *Rep` | QuaternionRep returns the two-dimensional irreducible representation of the quaternion group Q8 by Pauli-like matrices, with i, j, k acting as the… |
| `func RegularRep(g *Group) *Rep` | RegularRep returns the left regular representation of g, of dimension \|G\|: element x acts on the basis {e_y} by e_y ↦ e_{x·y}. |
| `func SignRepSymmetric(n int) *Rep` | SignRepSymmetric returns the one-dimensional sign representation of S_n, sending each permutation to its sign (+1 or -1). |
| `func TrivialRep(g *Group) *Rep` | TrivialRep returns the one-dimensional representation sending every element to [1]. |
| `func (r *Rep) Character() Character` | Character returns the character of r, the class function χ(i) = tr mats`i`. |
| `func (r *Rep) Dim() int` | Dim returns the dimension (degree) of the representation. |
| `func (r *Rep) DirectSum(s *Rep) *Rep` | DirectSum returns the direct sum r⊕s, the block-diagonal representation of dimension Dim(r)+Dim(s). |
| `func (r *Rep) Dual() *Rep` | Dual returns the dual (contragredient) representation, sending g to the inverse transpose of mats`g`. |
| `func (r *Rep) Group() *Group` | Group returns the group being represented. |
| `func (r *Rep) IsHomomorphism(tol float64) bool` | IsHomomorphism reports whether mats[i·j] = mats`i`·mats`j` for all i, j to within tol. |
| `func (r *Rep) IsUnitary(tol float64) bool` | IsUnitary reports whether every representing matrix is unitary to within tol. |
| `func (r *Rep) Matrix(i int) Matrix` | Matrix returns the representing matrix of element i. |
| `func (r *Rep) RepString() string` | RepString renders a short multi-line description of r for debugging. |
| `func (r *Rep) TensorProduct(s *Rep) *Rep` | TensorProduct returns the tensor (Kronecker) product r⊗s, of dimension Dim(r)·Dim(s). |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/grouprep).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
