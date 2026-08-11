# designs — implements combinatorial designs and finite geometry

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/designs.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/designs)

Package designs implements combinatorial designs and finite geometry.

- Balanced incomplete block designs: the parameter quintuple `BIBDParams` (v,
b, r, k, lambda) with the necessary divisibility conditions, Fisher's
inequality, symmetric designs and their complements, derived and residual
designs; a general `Design` type that stores an incidence structure and
recovers its parameters, replication numbers, concurrence (pairwise lambda)
matrix, dual and complement. - Latin squares and mutually orthogonal Latin
squares (`LatinSquare`, `MOLS`), Graeco-Latin squares and orthogonality
testing, built directly and from finite fields. - Steiner triple and quadruple
systems (`SteinerTripleSystem`, `SteinerQuadrupleSystem`) via the Bose
construction, Stinson's hill-climbing algorithm and the Boolean construction
on GF(2)^k. - Resolvable designs, parallel classes, one-factorizations of
complete graphs and round-robin tournament schedules (`OneFactorization`). -
Hadamard matrices (`HadamardMatrix`) by the Sylvester and Paley constructions,
conference matrices, and the symmetric 2-designs they produce. - Difference
sets (`DifferenceSet`) in cyclic groups, the Paley (quadratic-residue) and
Singer planar families, and development into symmetric designs. - Finite
projective and affine planes (`ProjectivePlane`, `AffinePlane`) of prime-power
order with full incidence, joins, meets and parallel classes.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/designs"
```

## Usage

This is the package's own `ExampleGaloisField`, so it compiles and its output is
asserted on every `go test ./designs/`.

```go
f, err := designs.NewGaloisField(8)
	if err != nil {
		panic(err)
	}
	fmt.Println("order:", f.Order(), "char:", f.Characteristic(), "degree:", f.Degree())

	for _, a := range []int{1, 2, 5, 7} {
		fmt.Printf("%d^7=%d;", a, f.Pow(a, 7))
	}
	fmt.Println()

	g, _ := f.PrimitiveElement()
	ord, _ := f.ElementOrder(g)
	fmt.Println("primitive element", g, "has order", ord)
```

```
order: 8 char: 2 degree: 3
1^7=1;2^7=1;5^7=1;7^7=1;
primitive element 2 has order 7
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AreOrthogonal(a, b LatinSquare) (bool, error)` | AreOrthogonal reports whether two Latin squares of the same order are orthogonal: superimposing them yields each of the n^2 ordered symbol pairs… |
| `func Binomial(n, k int) int` | Binomial returns the binomial coefficient C(n,k). |
| `func Divisors(n int) []int` | Divisors returns the positive divisors of n>=1 in increasing order. |
| `func EulerPhi(n int) int` | EulerPhi returns Euler's totient of n>=1, the count of integers in [1,n] coprime to n. |
| `func ExtendedGcd(a, b int) (g, x, y int)` | ExtendedGcd returns g = gcd(a,b) together with Bezout coefficients x, y such that a*x + b*y = g. |
| `func FactorPrimePower(q int) (p, k int, err error)` | FactorPrimePower writes q = p**k with p prime and k>=1, returning p and k. |
| `func Factorial(n int) (int, error)` | Factorial returns n! for 0<=n. |
| `func Factorize(n int) (primes, exponents []int)` | Factorize returns the prime factorization of n>=1 as sorted prime bases with their exponents. |
| `func FisherInequality(d *Design) (bool, error)` | FisherInequality reports whether the design d satisfies Fisher's inequality b>=v. |
| `func Gcd(a, b int) int` | Gcd returns the non-negative greatest common divisor of a and b. |
| `func GraecoLatinSquare(a, b LatinSquare) ([][]int, error)` | GraecoLatinSquare superimposes two orthogonal Latin squares a and b into a single array of ordered pairs encoded as a`i``j`*n + b`i``j`. |
| `func HadamardOrderValid(n int) bool` | HadamardOrderValid reports whether n is a possible Hadamard order: 1, 2, or a positive multiple of 4. |
| `func IntSqrt(n int) int` | IntSqrt returns the integer square root floor(sqrt(n)) for n>=0. |
| `func IrreduciblePolynomial(p, k int) ([]int, error)` | IrreduciblePolynomial returns a monic irreducible polynomial of degree k over GF(p), as coefficients indexed by degree (length k+1, leading… |
| `func IsDifferenceSetSlice(n int, elems []int) bool` | IsDifferenceSetSlice reports whether the given residues form a (n,k,lambda) difference set in Z_n. |
| `func IsIrreducible(f []int, p int) bool` | IsIrreducible reports whether the polynomial f (coefficients indexed by degree) is irreducible over the prime field GF(p), using Rabin's test. |
| `func IsMOLS(squares []LatinSquare) bool` | IsMOLS reports whether the given squares are pairwise mutually orthogonal Latin squares of a common order. |
| `func IsOneFactorization(m int, factors [][][2]int) bool` | IsOneFactorization reports whether the given factors form a one-factorization of K_m: every factor is a perfect matching and every edge of K_m… |
| `func IsParallelClass(d *Design, blockIdx []int) bool` | IsParallelClass reports whether the given block indices of the design form a parallel class: the referenced blocks are pairwise disjoint and together… |
| `func IsPerfectSquare(n int) bool` | IsPerfectSquare reports whether n is a perfect square (n>=0). |
| `func IsPrime(n int) bool` | IsPrime reports whether n is a prime number using trial division. |
| `func IsPrimePower(q int) bool` | IsPrimePower reports whether q equals p**k for a prime p and integer k>=1. |
| `func IsQuadraticResidue(a, p int) bool` | IsQuadraticResidue reports whether a is a non-zero quadratic residue modulo the odd prime p. |
| `func IsResolvableDesign(d *Design, classes [][]int) bool` | IsResolvableDesign reports whether the design admits the supplied resolution into parallel classes. |
| `func IsSteinerQuadrupleSystem(d *Design) bool` | IsSteinerQuadrupleSystem reports whether the design is a Steiner quadruple system: every 3-subset of points lies in exactly one block of size four. |
| `func IsSteinerTripleSystem(d *Design) bool` | IsSteinerTripleSystem reports whether the design is a Steiner triple system: a 2-(v,3,1) design in which every pair of points lies in exactly one… |
| `func JacobsthalMatrix(q int) ([][]int, error)` | JacobsthalMatrix returns the q-by-q Jacobsthal matrix Q of the field GF(q), with Q`i``j` = chi(i-j) where chi is the quadratic character. |
| `func Lcm(a, b int) int` | Lcm returns the least common multiple of a and b. |
| `func LegendreSymbol(a, p int) (int, error)` | LegendreSymbol returns the Legendre symbol (a\|p) for an odd prime p: it is 0 when p divides a, +1 when a is a non-zero quadratic residue modulo p,… |
| `func MaxMOLS(n int) int` | MaxMOLS returns the theoretical upper bound n-1 on the number of mutually orthogonal Latin squares of order n, for n>=2. |
| `func ModExp(base, exp, m int) int` | ModExp returns base**exp modulo m, computed by binary exponentiation, as a value in [0,m). |
| `func ModInverse(a, m int) (int, error)` | ModInverse returns the multiplicative inverse of a modulo m as a value in [0,m). |
| `func MultiplicativeOrder(a, m int) (int, error)` | MultiplicativeOrder returns the multiplicative order of a modulo m, the least e>=1 with a**e == 1 (mod m). |
| `func NearOneFactorization(m int) ([][][2]int, error)` | NearOneFactorization returns a near-one-factorization of the complete graph K_m on an odd number m of vertices: m near-perfect matchings, each… |
| `func NextPrime(n int) int` | NextPrime returns the smallest prime strictly greater than n. |
| `func NumDivisors(n int) int` | NumDivisors returns the number of positive divisors of n>=1. |
| `func OneFactorization(m int) ([][][2]int, error)` | OneFactorization returns a one-factorization of the complete graph K_m on m vertices (m even): a partition of its edges into m-1 perfect matchings… |
| `func OrthogonalArrayFromMOLS(squares []LatinSquare) ([][]int, error)` | OrthogonalArrayFromMOLS converts a set of m mutually orthogonal Latin squares of order n into the rows of the corresponding orthogonal array OA(n^2,… |
| `func PrevPrime(n int) (int, error)` | PrevPrime returns the largest prime strictly less than n, or an error when no such prime exists (n<=2). |
| `func PrimePowerBase(q int) (int, error)` | PrimePowerBase returns the prime p when q = p**k is a prime power, together with an error otherwise. |
| `func PrimePowerExponent(q int) (int, error)` | PrimePowerExponent returns the exponent k when q = p**k is a prime power, together with an error otherwise. |
| `func PrimesUpTo(n int) []int` | PrimesUpTo returns all primes p with p<=n in increasing order, using the sieve of Eratosthenes. |
| `func PrimitiveRoot(p int) (int, error)` | PrimitiveRoot returns the smallest primitive root modulo the prime p (a generator of the multiplicative group). |
| `func QuadraticCharacter(q, x int) (int, error)` | QuadraticCharacter returns the quadratic character chi(x) of the element x in GF(q): 0 when x is zero, +1 when x is a non-zero square, and -1 when x… |
| `func QuadraticNonResidues(p int) []int` | QuadraticNonResidues returns the quadratic non-residues modulo the odd prime p in increasing order. |
| `func QuadraticResidues(p int) []int` | QuadraticResidues returns the non-zero quadratic residues modulo the odd prime p in increasing order. |
| `func RoundRobinSchedule(teams int) ([][][2]int, error)` | RoundRobinSchedule returns a round-robin tournament schedule for the given number of teams: a list of rounds, each round a set of {home,away}… |
| `func SteinerQuadrupleAdmissible(v int) bool` | SteinerQuadrupleAdmissible reports whether v is an admissible order for a Steiner quadruple system S(3,4,v), i.e. |
| `func SteinerQuadrupleParameters(v int) (int, error)` | SteinerQuadrupleParameters returns the number of blocks v(v-1)(v-2)/24 of a Steiner quadruple system on v points. |
| `func SteinerSystemBlockCount(t, k, v int) (int, error)` | SteinerSystemBlockCount returns the number of blocks C(v,t)/C(k,t) of a Steiner system S(t,k,v), the number of k-subsets needed so every t-subset is… |
| `func SteinerTripleAdmissible(v int) bool` | SteinerTripleAdmissible reports whether v is an admissible order for a Steiner triple system, i.e. |
| `func SteinerTripleParameters(v int) (blocks int, err error)` | SteinerTripleParameters reports the number of blocks v(v-1)/6 of a Steiner triple system on v points. |
| `func SymmetricConferenceMatrix(q int) ([][]int, error)` | SymmetricConferenceMatrix returns the symmetric conference matrix of order q+1 built from GF(q) for a prime power q congruent to 1 modulo 4. |

### Types

| Type | What it is |
| --- | --- |
| `AffinePlane` | AffinePlane is the finite affine plane AG(2,q) of prime-power order q, with q^2 points and q^2+q lines partitioned into q+1 parallel classes of q… |
| `BIBDParams` | BIBDParams holds the five parameters of a balanced incomplete block design (BIBD, or 2-design): v points, b blocks, each point in r blocks, each… |
| `Design` | Design is an incidence structure: a finite set of points {0,1,...,Points-1} together with a list of blocks, each block being a sorted set of distinct… |
| `DifferenceSet` | DifferenceSet is a subset D of the cyclic group Z_n. |
| `GaloisField` | GaloisField is a finite field GF(q) with q = p**k elements. |
| `HadamardMatrix` | HadamardMatrix is a square matrix with entries in {-1,+1} whose rows are pairwise orthogonal, i.e. |
| `LatinSquare` | LatinSquare is an n-by-n array of symbols drawn from {0,...,n-1} in which every symbol occurs exactly once in each row and once in each column. |
| `ProjectivePlane` | ProjectivePlane is the finite Desarguesian projective plane PG(2,q) of prime-power order q. |
| `Resolution` | Resolution pairs a design with a partition of its blocks into parallel classes. |

<details>
<summary><code>AffinePlane</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewAffinePlane(q int) (*AffinePlane, error)` | NewAffinePlane constructs AG(2,q) for a prime power q>=2. |
| `func (ap *AffinePlane) AreParallel(l1, l2 int) bool` | AreParallel reports whether two distinct lines belong to the same parallel class (are disjoint). |
| `func (ap *AffinePlane) IncidenceDesign() *Design` | IncidenceDesign returns the affine plane as a resolvable 2-(q^2, q, 1) design whose blocks are the lines. |
| `func (ap *AffinePlane) IsAffinePlane() bool` | IsAffinePlane verifies that the incidence structure is a 2-(q^2,q,1) design with q+1 parallel classes. |
| `func (ap *AffinePlane) IsIncident(p, l int) bool` | IsIncident reports whether point index p lies on line index l. |
| `func (ap *AffinePlane) NumLines() int` | NumLines returns the number of lines q^2+q. |
| `func (ap *AffinePlane) NumParallelClasses() int` | NumParallelClasses returns the number of parallel classes q+1. |
| `func (ap *AffinePlane) NumPoints() int` | NumPoints returns the number of points q^2. |
| `func (ap *AffinePlane) Order() int` | Order returns the order q of the affine plane. |
| `func (ap *AffinePlane) ParallelClasses() [][]int` | ParallelClasses returns the parallel classes, each as a slice of line indices; lines in the same class are pairwise disjoint and together cover every… |
| `func (ap *AffinePlane) PointsOnLine(l int) []int` | PointsOnLine returns the indices of the q points on line l. |

</details>

<details>
<summary><code>BIBDParams</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AffinePlaneParams(n int) (BIBDParams, error)` | AffinePlaneParams returns the BIBD parameters 2-(n^2, n, 1) of an affine plane of order n, for n>=2. |
| `func HadamardDesignParams(t int) (BIBDParams, error)` | HadamardDesignParams returns the parameters 2-(4t-1, 2t-1, t-1) of the symmetric design associated with a Hadamard matrix of order 4t, for t>=2. |
| `func NewBIBDParams(v, k, lambda int) (BIBDParams, error)` | NewBIBDParams constructs a parameter set from v, k and lambda, deriving the replication number r = lambda*(v-1)/(k-1) and block count b = v*r/k. |
| `func ProjectivePlaneParams(n int) (BIBDParams, error)` | ProjectivePlaneParams returns the BIBD parameters 2-(n^2+n+1, n+1, 1) of a projective plane of order n, for n>=2. |
| `func SteinerTripleParams(v int) (BIBDParams, error)` | SteinerTripleParams returns the BIBD parameters 2-(v,3,1) of a Steiner triple system on v points. |
| `func (p BIBDParams) BlockIntersection() (int, error)` | BlockIntersection returns the number of points common to any two distinct blocks of a symmetric design, which is the constant lambda. |
| `func (p BIBDParams) Complement() (BIBDParams, error)` | Complement returns the parameters of the complementary design, obtained by replacing each block with its complement. |
| `func (p BIBDParams) Derived() (BIBDParams, error)` | Derived returns the parameters of the derived design of a symmetric design with respect to a block: it is a (k, k-1... |
| `func (p BIBDParams) FisherInequalityHolds() bool` | FisherInequalityHolds reports whether Fisher's inequality b>=v holds for the parameters. |
| `func (p BIBDParams) IsSymmetric() bool` | IsSymmetric reports whether the design is symmetric, meaning b==v (equivalently r==k). |
| `func (p BIBDParams) Order() int` | Order returns the order n = r - lambda of the design, a quantity that governs much of the theory of symmetric designs. |
| `func (p BIBDParams) Residual() (BIBDParams, error)` | Residual returns the parameters of the residual design of a symmetric 2-(v,k,lambda) design: a 2-(v-k, k-lambda, lambda) design. |
| `func (p BIBDParams) SatisfiesNecessary() bool` | SatisfiesNecessary reports whether the parameters meet the necessary existence conditions: the divisibility identities (Valid) together with Fisher's… |
| `func (p BIBDParams) Valid() bool` | Valid reports whether the parameters satisfy the two fundamental counting identities r*(k-1) = lambda*(v-1) and b*k = v*r, along with the basic range… |

</details>

<details>
<summary><code>Design</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BooleanQuadrupleSystem(k int) (*Design, error)` | BooleanQuadrupleSystem constructs the Steiner quadruple system on the 2**k points of the vector space GF(2)^k, whose blocks are the 4-subsets… |
| `func BoseSteinerTripleSystem(n int) (*Design, error)` | BoseSteinerTripleSystem constructs a Steiner triple system on n = 6t+3 points by the Bose construction, using the idempotent commutative quasigroup… |
| `func HadamardToDesign(h HadamardMatrix) (*Design, error)` | HadamardToDesign converts a normalized Hadamard matrix of order 4t into the symmetric 2-(4t-1, 2t-1, t-1) design obtained by deleting the first row… |
| `func NewDesign(points int, blocks [][]int) (*Design, error)` | NewDesign builds a Design on the given number of points from the supplied blocks. |
| `func SteinerQuadrupleSystem(v int) (*Design, error)` | SteinerQuadrupleSystem constructs a Steiner quadruple system on v points. |
| `func SteinerTripleSystem(v int, rng *rand.Rand) (*Design, error)` | SteinerTripleSystem constructs a Steiner triple system on v points using Stinson's hill-climbing algorithm, driven by the supplied random source so… |
| `func (d *Design) BlockContains(j, i int) bool` | BlockContains reports whether block index j contains point i. |
| `func (d *Design) BlockSizes() []int` | BlockSizes returns the size of every block, indexed by block. |
| `func (d *Design) Complement() *Design` | Complement returns the complementary design on the same points, whose blocks are the set complements of the original blocks. |
| `func (d *Design) Concurrence(i, j int) int` | Concurrence returns the number of blocks that contain both point i and point j. |
| `func (d *Design) ConcurrenceMatrix() [][]int` | ConcurrenceMatrix returns the symmetric points-by-points matrix whose (i,j) entry is the number of blocks containing both points i and j; the… |
| `func (d *Design) Dual() *Design` | Dual returns the dual design, in which the roles of points and blocks are exchanged: the dual has one point per original block, and one block per… |
| `func (d *Design) IncidenceMatrix() [][]int` | IncidenceMatrix returns the points-by-blocks 0/1 incidence matrix M, where M`i``j` is 1 exactly when point i lies in block j. |
| `func (d *Design) IsBIBD() (BIBDParams, bool)` | IsBIBD reports whether the design is a balanced incomplete block design and, when it is, returns its parameters. |
| `func (d *Design) IsPairBalanced() (int, bool)` | IsPairBalanced reports whether every pair of distinct points lies in the same number of blocks lambda, returning lambda and true, or (0,false)… |
| `func (d *Design) IsPointRegular() (int, bool)` | IsPointRegular reports whether every point has the same replication number, returning that common value r and true, or (0,false) otherwise. |
| `func (d *Design) IsSimple() bool` | IsSimple reports whether the design has no repeated blocks. |
| `func (d *Design) IsUniform() (int, bool)` | IsUniform reports whether all blocks have the same size k, returning that common size and true, or (0,false) otherwise. |
| `func (d *Design) NumBlocks() int` | NumBlocks returns the number of blocks b of the design. |
| `func (d *Design) NumPoints() int` | NumPoints returns the number of points v of the design. |
| `func (d *Design) Parameters() (BIBDParams, error)` | Parameters attempts to recover the balanced incomplete block design parameters (v,b,r,k,lambda) of the design. |
| `func (d *Design) PointDegree(i int) int` | PointDegree returns the number of blocks containing point i (its replication number). |
| `func (d *Design) ReplicationNumbers() []int` | ReplicationNumbers returns the replication number of every point, indexed by point. |

</details>

<details>
<summary><code>DifferenceSet</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewDifferenceSet(n int, elems []int) (*DifferenceSet, error)` | NewDifferenceSet builds a subset of Z_n from the given residues, reducing them modulo n, removing duplicates and sorting. |
| `func QuadraticResidueDifferenceSet(p int) (*DifferenceSet, error)` | QuadraticResidueDifferenceSet returns the Paley difference set of quadratic residues in Z_p for a prime p congruent to 3 modulo 4, a (p, (p-1)/2,… |
| `func SingerDifferenceSet(q int) (*DifferenceSet, error)` | SingerDifferenceSet returns the Singer planar difference set of order q, a (q^2+q+1, q+1, 1) difference set constructed from the trace-zero… |
| `func (d *DifferenceSet) Complement() *DifferenceSet` | Complement returns the complement of D within Z_n. |
| `func (d *DifferenceSet) Develop() *Design` | Develop returns the symmetric design obtained by developing D: the point set is Z_n and the blocks are the n translates D+g for g in Z_n. |
| `func (d *DifferenceSet) DifferenceCounts() []int` | DifferenceCounts returns a slice indexed by g in [0,n) giving the number of ordered pairs (x,y) with x,y in D and x-y == g (mod n). |
| `func (d *DifferenceSet) IsDifferenceSet() bool` | IsDifferenceSet reports whether D is a (n,k,lambda) difference set: every non-zero group element occurs equally often as a difference. |
| `func (d *DifferenceSet) IsMultiplier(t int) bool` | IsMultiplier reports whether t is a numerical multiplier of D, meaning t is coprime to n and t*D is a translate D+s of D for some s. |
| `func (d *DifferenceSet) IsPlanar() bool` | IsPlanar reports whether D is a planar (lambda=1) difference set, which develops into a projective plane. |
| `func (d *DifferenceSet) Lambda() (int, bool)` | Lambda returns the common difference multiplicity lambda together with true when D is a difference set, or (0,false) otherwise. |
| `func (d *DifferenceSet) Multipliers() []int` | Multipliers returns all numerical multipliers t in [1,n) of the difference set, i.e. |
| `func (d *DifferenceSet) Multiply(t int) *DifferenceSet` | Multiply returns the dilated set t*D (mod n). |
| `func (d *DifferenceSet) Order() int` | Order returns the order n of the ambient cyclic group. |
| `func (d *DifferenceSet) Parameters() (n, k, lambda int, err error)` | Parameters returns the (n,k,lambda) parameters of the difference set. |
| `func (d *DifferenceSet) Size() int` | Size returns the number of elements k of the set. |
| `func (d *DifferenceSet) Translate(s int) *DifferenceSet` | Translate returns the shifted set D+s (mod n). |

</details>

<details>
<summary><code>GaloisField</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGaloisField(q int) (*GaloisField, error)` | NewGaloisField constructs GF(q) for a prime power q = p**k. |
| `func (f *GaloisField) Add(a, b int) int` | Add returns the field sum a+b. |
| `func (f *GaloisField) Characteristic() int` | Characteristic returns the prime characteristic p of the field. |
| `func (f *GaloisField) Degree() int` | Degree returns the extension degree k, so that q = p**k. |
| `func (f *GaloisField) Div(a, b int) (int, error)` | Div returns a/b for non-zero b. |
| `func (f *GaloisField) ElementOrder(a int) (int, error)` | ElementOrder returns the multiplicative order of a non-zero element a, the least e>=1 with a**e == 1. |
| `func (f *GaloisField) Elements() []int` | Elements returns all q elements of the field, 0,1,...,q-1. |
| `func (f *GaloisField) FromCoefficients(c []int) int` | FromCoefficients returns the field element whose polynomial coefficients (constant term first) are the given values, reduced modulo p. |
| `func (f *GaloisField) Inv(a int) (int, error)` | Inv returns the multiplicative inverse of a non-zero element a. |
| `func (f *GaloisField) IsPrimitiveElement(a int) bool` | IsPrimitiveElement reports whether the non-zero element a generates the multiplicative group, i.e. |
| `func (f *GaloisField) Mul(a, b int) int` | Mul returns the field product a*b. |
| `func (f *GaloisField) Neg(a int) int` | Neg returns the field additive inverse -a. |
| `func (f *GaloisField) Norm(a int) int` | Norm returns the field norm of a down to the prime field, defined as the product a * a**p * ... |
| `func (f *GaloisField) One() int` | One returns the multiplicative identity of the field. |
| `func (f *GaloisField) Order() int` | Order returns the number of elements q = p**k of the field. |
| `func (f *GaloisField) Poly() []int` | Poly returns a copy of the defining monic irreducible polynomial (nil for a prime field's trivial polynomial is instead returned as [0,1]). |
| `func (f *GaloisField) Pow(a, n int) int` | Pow returns a**n for n>=0 in the field. |
| `func (f *GaloisField) PrimitiveElement() (int, error)` | PrimitiveElement returns a generator of the multiplicative group of the field. |
| `func (f *GaloisField) Sub(a, b int) int` | Sub returns the field difference a-b. |
| `func (f *GaloisField) ToCoefficients(e int) []int` | ToCoefficients returns the polynomial coefficients (constant term first, length k) of the field element e. |
| `func (f *GaloisField) Trace(a int) int` | Trace returns the field trace of a down to the prime field, defined as a + a**p + a**(p^2) + ... |
| `func (f *GaloisField) Zero() int` | Zero returns the additive identity of the field. |

</details>

<details>
<summary><code>HadamardMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func KroneckerProduct(a, b HadamardMatrix) HadamardMatrix` | KroneckerProduct returns the Kronecker (tensor) product of two Hadamard matrices, itself a Hadamard matrix of order equal to the product of the two… |
| `func PaleyConstructionI(q int) (HadamardMatrix, error)` | PaleyConstructionI returns the Paley type-I Hadamard matrix of order q+1, defined for a prime power q congruent to 3 modulo 4. |
| `func PaleyConstructionII(q int) (HadamardMatrix, error)` | PaleyConstructionII returns the Paley type-II Hadamard matrix of order 2(q+1), defined for a prime power q congruent to 1 modulo 4. |
| `func SylvesterHadamard(k int) (HadamardMatrix, error)` | SylvesterHadamard returns the Sylvester Hadamard matrix of order 2**k, built recursively from H(1)=[1] by H(2m)=[[H,H],[H,-H]]. |
| `func (h HadamardMatrix) IsHadamard() bool` | IsHadamard reports whether the matrix has entries in {-1,+1} and satisfies H*H^T = n*I, the defining property of a Hadamard matrix. |
| `func (h HadamardMatrix) IsNormalized() bool` | IsNormalized reports whether the first row and first column consist entirely of +1 entries. |
| `func (h HadamardMatrix) Normalize() HadamardMatrix` | Normalize returns an equivalent Hadamard matrix whose first row and first column are all +1, obtained by negating the rows and columns that begin… |
| `func (h HadamardMatrix) Order() int` | Order returns the side length of the matrix. |
| `func (h HadamardMatrix) Row(i int) []int` | Row returns a copy of row i. |
| `func (h HadamardMatrix) Transpose() HadamardMatrix` | Transpose returns the transpose of the matrix, which is again Hadamard. |

</details>

<details>
<summary><code>LatinSquare</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CyclicLatinSquare(n int) (LatinSquare, error)` | CyclicLatinSquare returns the Latin square of order n given by l`i``j` = (i+j) mod n, the Cayley table of the cyclic group Z_n. |
| `func LatinSquareFromField(q, a int) (LatinSquare, error)` | LatinSquareFromField returns the Latin square l`x``y` = a*x + y over GF(q) for a fixed non-zero field element a, with rows and columns indexed by… |
| `func MOLS(q int) ([]LatinSquare, error)` | MOLS returns a complete set of q-1 mutually orthogonal Latin squares of order q, constructed from the finite field GF(q) as L_a`x``y` = a*x + y for… |
| `func NewLatinSquare(rows [][]int) (LatinSquare, error)` | NewLatinSquare validates that the supplied array is a Latin square and returns it. |
| `func (l LatinSquare) At(i, j int) int` | At returns the symbol in row i, column j. |
| `func (l LatinSquare) Col(j int) []int` | Col returns a copy of column j. |
| `func (l LatinSquare) IsIdempotent() bool` | IsIdempotent reports whether l`i``i`==i for every i, i.e. |
| `func (l LatinSquare) IsLatin() bool` | IsLatin reports whether the array is a valid Latin square: square shape, entries in range, and each symbol appearing once per row and once per column. |
| `func (l LatinSquare) IsReduced() bool` | IsReduced reports whether the square is in reduced (standard) form: its first row and first column are both 0,1,...,n-1 in order. |
| `func (l LatinSquare) IsSymmetric() bool` | IsSymmetric reports whether the square equals its transpose (a commutative quasigroup). |
| `func (l LatinSquare) Order() int` | Order returns the side length n of the square. |
| `func (l LatinSquare) Row(i int) []int` | Row returns a copy of row i. |
| `func (l LatinSquare) Transpose() LatinSquare` | Transpose returns the transpose of the square, which is again a Latin square. |

</details>

<details>
<summary><code>ProjectivePlane</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FanoPlane() *ProjectivePlane` | FanoPlane returns the Fano plane PG(2,2), the unique projective plane of order 2 with 7 points and 7 lines. |
| `func NewProjectivePlane(q int) (*ProjectivePlane, error)` | NewProjectivePlane constructs PG(2,q) for a prime power q>=2. |
| `func (pl *ProjectivePlane) IncidenceDesign() *Design` | IncidenceDesign returns the points-and-lines incidence structure as a Design: one point per plane point and one block per line listing the points on… |
| `func (pl *ProjectivePlane) IsIncident(p, l int) bool` | IsIncident reports whether point index p lies on line index l, i.e. |
| `func (pl *ProjectivePlane) IsProjectivePlane() bool` | IsProjectivePlane verifies the incidence axioms: the plane has q^2+q+1 points and lines, every line has q+1 points, and every pair of points lies on… |
| `func (pl *ProjectivePlane) Line(j int) [3]int` | Line returns the homogeneous coordinates of line index j. |
| `func (pl *ProjectivePlane) LineThroughPoints(p1, p2 int) (int, error)` | LineThroughPoints returns the index of the unique line joining two distinct points. |
| `func (pl *ProjectivePlane) LinesThroughPoint(p int) []int` | LinesThroughPoint returns the indices of the q+1 lines incident with point p. |
| `func (pl *ProjectivePlane) MeetOfLines(l1, l2 int) (int, error)` | MeetOfLines returns the index of the unique point common to two distinct lines. |
| `func (pl *ProjectivePlane) NumLines() int` | NumLines returns the number of lines q^2+q+1. |
| `func (pl *ProjectivePlane) NumPoints() int` | NumPoints returns the number of points q^2+q+1. |
| `func (pl *ProjectivePlane) Order() int` | Order returns the order q of the plane. |
| `func (pl *ProjectivePlane) Point(i int) [3]int` | Point returns the homogeneous coordinates of point index i. |
| `func (pl *ProjectivePlane) PointsOnLine(l int) []int` | PointsOnLine returns the indices of the q+1 points incident with line l. |

</details>

<details>
<summary><code>Resolution</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AffineResolution(q int) (*Resolution, error)` | AffineResolution returns the affine plane AG(2,q) as a resolvable design together with its natural resolution into q+1 parallel classes (the parallel… |
| `func NewResolution(d *Design, classes [][]int) (*Resolution, error)` | NewResolution validates that the given classes partition all blocks of the design into parallel classes and returns the resolution. |
| `func (r *Resolution) Class(i int) []int` | Class returns the block indices making up parallel class i. |
| `func (r *Resolution) IsValid() bool` | IsValid reports whether the resolution is well formed: every class is a parallel class and together they partition the blocks. |
| `func (r *Resolution) NumClasses() int` | NumClasses returns the number of parallel classes in the resolution. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/designs).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
