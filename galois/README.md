# galois — implements arithmetic and structure theory for finite fields GF(p^n)

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/galois.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/galois)

Package galois implements arithmetic and structure theory for finite fields
GF(p^n), built entirely on the Go standard library.

A finite field (Galois field) of order q = p^n exists for every prime p and
positive integer n. This package models both the prime field GF(p) and the
extension fields GF(p^n), providing element arithmetic (add, subtract,
multiply, invert, divide, exponentiate) together with the deeper structural
tools that make finite fields useful in coding theory and cryptography.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/galois"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AddMod(a, b, p *big.Int) *big.Int` | AddMod returns (a + b) mod p in the range [0, p). |
| `func BabyStepGiantStep(base, target *Element) (*big.Int, error)` | BabyStepGiantStep is an alias for DiscreteLog. |
| `func CyclotomicCoset(q *big.Int, s, n int) []int64` | CyclotomicCoset returns the q-cyclotomic coset of s modulo q^n - 1, that is the sorted set {s, s·q, s·q², …} reduced modulo q^n - 1. |
| `func CyclotomicCosets(q *big.Int, n int) [][]int64` | CyclotomicCosets returns the partition of {0, 1, …, q^n-2} into q-cyclotomic cosets modulo q^n - 1, each coset sorted and the list ordered by… |
| `func DiscreteLog(base, target *Element) (*big.Int, error)` | DiscreteLog returns the discrete logarithm x of target to the base, the least non-negative x with base^x = target, computed by the baby-step… |
| `func DistinctFactorCount(f *Poly) int` | DistinctFactorCount returns the number of distinct irreducible factors of f. |
| `func DivMod(a, b, p *big.Int) (*big.Int, error)` | DivMod returns (a / b) mod p, that is a times the inverse of b, or an error when b is not invertible. |
| `func DivisorsInt(n *big.Int) []*big.Int` | DivisorsInt returns all positive divisors of \|n\| in increasing order. |
| `func EulerPhi(n *big.Int) *big.Int` | EulerPhi returns Euler's totient φ(n): the count of integers in [1, n] coprime to n. |
| `func ExtendedGcd(a, b *big.Int) (g, x, y *big.Int)` | ExtendedGcd returns g = gcd(a, b) together with Bézout coefficients x and y satisfying a*x + b*y = g. |
| `func FactorCount(f *Poly) int` | FactorCount returns the number of irreducible factors of f counted with multiplicity. |
| `func FieldNorm(a *Element) *big.Int` | FieldNorm returns the norm of a to the prime field. |
| `func FieldTrace(a *Element) *big.Int` | FieldTrace returns the trace of a to the prime field. |
| `func Gcd(a, b *big.Int) *big.Int` | Gcd returns the non-negative greatest common divisor of a and b. |
| `func IntPow(base *big.Int, exp int) *big.Int` | IntPow returns base raised to the power exp (exp >= 0) as an exact integer. |
| `func InvMod(a, p *big.Int) (*big.Int, error)` | InvMod returns the multiplicative inverse of a modulo p, or an error when a is not invertible (gcd(a, p) != 1). |
| `func IsIrreducible(f *Poly) bool` | IsIrreducible reports whether f is irreducible over GF(p) using Rabin's test. |
| `func IsPrimeInt(n *big.Int) bool` | IsPrimeInt reports whether n is a prime integer. |
| `func IsPrimitivePoly(f *Poly) bool` | IsPrimitivePoly reports whether f is a primitive polynomial over GF(p): a monic irreducible polynomial whose root x generates the whole… |
| `func IsPrimitiveRoot(a, p *big.Int) bool` | IsPrimitiveRoot reports whether a is a primitive root modulo the prime p, i.e. |
| `func IsQuadraticResidue(a, p *big.Int) bool` | IsQuadraticResidue reports whether a is a square modulo the prime p. |
| `func IsSquareFree(f *Poly) bool` | IsSquareFree reports whether f has no repeated irreducible factor. |
| `func IsSubfieldDegree(n, d int) bool` | IsSubfieldDegree reports whether GF(p^d) embeds in GF(p^n), i.e. |
| `func Jacobi(a, n *big.Int) int` | Jacobi returns the Jacobi symbol (a/n) for odd n > 0, one of -1, 0, +1. |
| `func Lcm(a, b *big.Int) *big.Int` | Lcm returns the non-negative least common multiple of a and b. |
| `func Legendre(a, p *big.Int) int` | Legendre returns the Legendre symbol (a/p) for an odd prime p: it is 0 when p divides a, +1 when a is a non-zero quadratic residue, and -1 otherwise. |
| `func MobiusMu(n *big.Int) int` | MobiusMu returns the Möbius function μ(n): 0 if n is divisible by a square greater than 1, otherwise (-1)^k where k is the number of distinct… |
| `func Mod(a, p *big.Int) *big.Int` | Mod returns a reduced into the canonical range [0, p) for the modulus p > 0. |
| `func MulMod(a, b, p *big.Int) *big.Int` | MulMod returns (a * b) mod p in the range [0, p). |
| `func MultiplicativeOrder(a, p *big.Int) (*big.Int, error)` | MultiplicativeOrder returns the multiplicative order of a in the group of units modulo the prime p: the smallest positive k with a^k ≡ 1 (mod p). |
| `func NegMod(a, p *big.Int) *big.Int` | NegMod returns (-a) mod p in the range [0, p). |
| `func NextPrime(n *big.Int) *big.Int` | NextPrime returns the smallest prime strictly greater than n. |
| `func NumDivisors(n *big.Int) int` | NumDivisors returns the number of positive divisors of \|n\|. |
| `func NumberOfIrreducibles(p *big.Int, n int) *big.Int` | NumberOfIrreducibles returns the number of monic irreducible polynomials of degree n over GF(p), via the necklace formula (1/n)·Σ_{d\|n}… |
| `func PolyEqual(a, b *Poly) bool` | PolyEqual reports whether two polynomials are equal. |
| `func PowMod(a, e, p *big.Int) (*big.Int, error)` | PowMod returns a raised to the power e modulo p. |
| `func PrimeFactors(n *big.Int) []*big.Int` | PrimeFactors returns the distinct prime factors of \|n\| in increasing order. |
| `func PrimitiveRoot(p *big.Int) (*big.Int, error)` | PrimitiveRoot returns the smallest primitive root modulo the prime p. |
| `func SqrtMod(a, p *big.Int) (*big.Int, error)` | SqrtMod returns a square root r of a modulo the prime p, so that r*r ≡ a (mod p). |
| `func SubMod(a, b, p *big.Int) *big.Int` | SubMod returns (a - b) mod p in the range [0, p). |
| `func SubfieldDegrees(n int) []int` | SubfieldDegrees returns the degrees of the subfields of GF(p^n): the positive divisors of n. |

### Types

| Type | What it is |
| --- | --- |
| `Element` | Element is a member of a Field, represented by its unique residue polynomial of degree less than the field degree. |
| `Field` | Field represents the finite field GF(p^n). |
| `MatModP` | MatModP is a dense matrix with entries in the prime field GF(p). |
| `Poly` | Poly is a dense univariate polynomial over the prime field GF(p). |
| `PolyFactor` | PolyFactor pairs an irreducible (or square-free) monic polynomial factor with the multiplicity at which it occurs. |
| `PrimePower` | PrimePower records a prime together with its exponent in a factorisation, so that n = Prime^Exp times the other factors. |

<details>
<summary><code>Element</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ConjugateElements(a *Element) []*Element` | ConjugateElements returns the Galois conjugates of a. |
| `func (a *Element) Add(b *Element) *Element` | Add returns a + b in the field. |
| `func (a *Element) Coeffs() []*big.Int` | Coeffs returns the little-endian coefficient vector of the element, padded to the field degree n. |
| `func (a *Element) Conjugates() []*Element` | Conjugates returns the Galois conjugates of a, namely a, a^p, a^(p^2), … up to but not including the first repetition. |
| `func (a *Element) Copy() *Element` | Copy returns an independent copy of the element. |
| `func (a *Element) Div(b *Element) (*Element, error)` | Div returns a / b in the field, or an error when b is zero. |
| `func (a *Element) Equal(b *Element) bool` | Equal reports whether a and b are the same element of the same field. |
| `func (a *Element) Frobenius() *Element` | Frobenius returns the image of a under the Frobenius automorphism x ↦ x^p, the canonical generator of the Galois group of GF(p^n) over GF(p). |
| `func (a *Element) FrobeniusPow(k int) *Element` | FrobeniusPow returns the k-fold Frobenius image a^(p^k). |
| `func (a *Element) Inv() (*Element, error)` | Inv returns the multiplicative inverse of a, or an error when a is zero. |
| `func (a *Element) IsOne() bool` | IsOne reports whether a is the multiplicative identity. |
| `func (a *Element) IsPrimitive() bool` | IsPrimitive reports whether a generates the whole multiplicative group of the field, i.e. |
| `func (a *Element) IsZero() bool` | IsZero reports whether a is the additive identity. |
| `func (a *Element) MinimalPoly() *Poly` | MinimalPoly returns the minimal polynomial of a over GF(p): the unique monic polynomial of least degree with a as a root, equal to the product of (x… |
| `func (a *Element) Mul(b *Element) *Element` | Mul returns a * b in the field. |
| `func (a *Element) Neg() *Element` | Neg returns the additive inverse of a. |
| `func (a *Element) Norm() *big.Int` | Norm returns the field norm N(a) = a^((p^n-1)/(p-1)) = a · a^p · … · a^(p^(n-1)) of a down to the prime field GF(p), returned as an integer in… |
| `func (a *Element) Order() (*big.Int, error)` | Order returns the multiplicative order of the non-zero element a: the least positive k with a^k = 1. |
| `func (a *Element) Pow(e *big.Int) (*Element, error)` | Pow returns a raised to the integer power e. |
| `func (a *Element) String() string` | String renders the element via its residue polynomial. |
| `func (a *Element) Sub(b *Element) *Element` | Sub returns a - b in the field. |
| `func (a *Element) ToInt() *big.Int` | ToInt returns the integer index in [0, p^n) obtained by reading the element's coefficient vector as base-p digits. |
| `func (a *Element) Trace() *big.Int` | Trace returns the field trace Tr(a) = a + a^p + … + a^(p^(n-1)) of a down to the prime field GF(p), returned as an integer in [0, p). |

</details>

<details>
<summary><code>Field</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewField(p *big.Int, mod *Poly) (*Field, error)` | NewField returns GF(p^n) defined by the given monic modulus polynomial, which must be irreducible of degree n >= 1 over GF(p). |
| `func NewFieldFromDegree(p *big.Int, n int) (*Field, error)` | NewFieldFromDegree returns GF(p^n) using an automatically discovered monic irreducible modulus of degree n. |
| `func NewPrimeField(p *big.Int) (*Field, error)` | NewPrimeField returns the prime field GF(p). |
| `func (f *Field) Add(a, b *Element) *Element` | Add returns a + b, a convenience wrapper over Element.Add. |
| `func (f *Field) Characteristic() *big.Int` | Characteristic returns the prime characteristic p of the field. |
| `func (f *Field) Degree() int` | Degree returns the extension degree n, so the field is GF(p^n). |
| `func (f *Field) Div(a, b *Element) (*Element, error)` | Div returns a / b. |
| `func (f *Field) Element(coeffs ...int64) *Element` | Element returns the field element for the given little-endian int64 coefficients over GF(p), reduced into the field. |
| `func (f *Field) ElementFromCoeffs(coeffs []*big.Int) *Element` | ElementFromCoeffs returns the field element for the given little-endian big.Int coefficients, reduced into the field. |
| `func (f *Field) ElementFromInt(k *big.Int) *Element` | ElementFromInt maps an integer index k in [0, p^n) to a field element by writing k in base p as the element's coefficient vector. |
| `func (f *Field) Elements() []*Element` | Elements enumerates every element of the field in index order. |
| `func (f *Field) FrobeniusAutomorphism() func(*Element) *Element` | FrobeniusAutomorphism returns the Frobenius map x ↦ x^p of the field as a Go function. |
| `func (f *Field) FromPoly(p *Poly) *Element` | FromPoly returns the field element represented by the polynomial p reduced modulo the field modulus. |
| `func (f *Field) Generator() (*Element, error)` | Generator is a synonym for PrimitiveElement, returning a multiplicative generator of the field. |
| `func (f *Field) Inv(a *Element) (*Element, error)` | Inv returns the inverse of a. |
| `func (f *Field) Modulus() *Poly` | Modulus returns a copy of the field's defining irreducible polynomial. |
| `func (f *Field) Mul(a, b *Element) *Element` | Mul returns a * b. |
| `func (f *Field) Neg(a *Element) *Element` | Neg returns -a. |
| `func (f *Field) NonzeroElements() []*Element` | NonzeroElements enumerates every non-zero element of the field in index order. |
| `func (f *Field) One() *Element` | One returns the multiplicative identity of the field. |
| `func (f *Field) Order() *big.Int` | Order returns the number of elements in the field, p^n. |
| `func (f *Field) Pow(a *Element, e *big.Int) (*Element, error)` | Pow returns a^e. |
| `func (f *Field) PrimitiveElement() (*Element, error)` | PrimitiveElement returns a generator of the multiplicative group of the field: an element of order p^n - 1. |
| `func (f *Field) PrimitiveElementX() *Element` | PrimitiveElementX returns the residue class of x, the canonical generator of GF(p^n) as GF(p)`x`/(mod). |
| `func (f *Field) RandomElement(r *rand.Rand) *Element` | RandomElement returns a uniformly random element drawn from r. |
| `func (f *Field) String() string` | String describes the field as "GF(p^n)". |
| `func (f *Field) Sub(a, b *Element) *Element` | Sub returns a - b. |
| `func (f *Field) SubfieldOrders() []*big.Int` | SubfieldOrders returns the orders p^d of all subfields of the field, in increasing order. |
| `func (f *Field) Subfields() []int` | Subfields returns the degrees d dividing n for which GF(p^d) is a subfield of the field GF(p^n), in increasing order. |
| `func (f *Field) Zero() *Element` | Zero returns the additive identity of the field. |

</details>

<details>
<summary><code>MatModP</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IdentityModP(p *big.Int, n int) *MatModP` | IdentityModP returns the n×n identity matrix over GF(p). |
| `func NewMatModP(p *big.Int, rows, cols int) *MatModP` | NewMatModP returns a rows×cols zero matrix over GF(p). |
| `func (m *MatModP) Clone() *MatModP` | Clone returns an independent copy of the matrix. |
| `func (m *MatModP) Get(i, j int) *big.Int` | Get returns the entry at (i, j). |
| `func (m *MatModP) Mul(other *MatModP) *MatModP` | Mul returns the matrix product m*other over GF(p). |
| `func (m *MatModP) NullSpace() [][]*big.Int` | NullSpace returns a basis for the right null space {x : m·x = 0} over GF(p), each basis vector represented as a length-Cols slice. |
| `func (m *MatModP) Rank() int` | Rank returns the rank of the matrix over GF(p). |
| `func (m *MatModP) RowReduce() (*MatModP, int)` | RowReduce returns the reduced row-echelon form of the matrix together with its rank. |
| `func (m *MatModP) Set(i, j int, v *big.Int)` | Set stores v (reduced mod p) at (i, j). |
| `func (m *MatModP) Sub(other *MatModP) *MatModP` | Sub returns the entrywise difference m-other over GF(p). |
| `func (m *MatModP) Transpose() *MatModP` | Transpose returns the transpose of the matrix. |

</details>

<details>
<summary><code>Poly</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AllIrreduciblePolys(p *big.Int, n int) []*Poly` | AllIrreduciblePolys returns every monic irreducible polynomial of degree n over GF(p). |
| `func CanonicalIrreduciblePoly(p *big.Int, n int) (*Poly, error)` | CanonicalIrreduciblePoly returns the deterministic canonical (lexicographically least by coefficient index) monic irreducible polynomial of degree n… |
| `func ConwayStylePoly(p *big.Int, n int) (*Poly, error)` | ConwayStylePoly returns a deterministic Conway-style defining polynomial for GF(p^n): the lexicographically least primitive polynomial of degree n… |
| `func FindIrreducible(p *big.Int, n int) (*Poly, error)` | FindIrreducible returns the first monic irreducible polynomial of degree n over GF(p), scanning candidates by ascending coefficient index. |
| `func FindPrimitivePoly(p *big.Int, n int) (*Poly, error)` | FindPrimitivePoly returns the first primitive polynomial of degree n over GF(p), scanning candidates by ascending coefficient index. |
| `func MinimalPolynomial(a *Element) *Poly` | MinimalPolynomial returns the minimal polynomial of a over GF(p). |
| `func NewPoly(p *big.Int, coeffs ...int64) *Poly` | NewPoly builds a polynomial over GF(p) from int64 coefficients given little-endian (constant term first). |
| `func NewPolyBig(p *big.Int, coeffs []*big.Int) *Poly` | NewPolyBig builds a polynomial over GF(p) from big.Int coefficients given little-endian. |
| `func PolyConst(p, c *big.Int) *Poly` | PolyConst returns the constant polynomial c over GF(p). |
| `func PolyGcd(a, b *Poly) *Poly` | PolyGcd returns the monic greatest common divisor of a and b. |
| `func PolyLcm(a, b *Poly) *Poly` | PolyLcm returns the monic least common multiple of a and b. |
| `func PolyMonomial(p, coeff *big.Int, deg int) *Poly` | PolyMonomial returns the monomial coeff·x^deg over GF(p). |
| `func PolyMulMod(a, b, m *Poly) *Poly` | PolyMulMod returns (a * b) mod m over GF(p). |
| `func PolyOne(p *big.Int) *Poly` | PolyOne returns the constant polynomial 1 over GF(p). |
| `func PolyPowMod(a *Poly, e *big.Int, m *Poly) (*Poly, error)` | PolyPowMod returns a^e mod m over GF(p) for a non-negative exponent e. |
| `func PolyX(p *big.Int) *Poly` | PolyX returns the monomial x over GF(p). |
| `func PolyZero(p *big.Int) *Poly` | PolyZero returns the zero polynomial over GF(p). |
| `func PthRootPoly(f *Poly) *Poly` | PthRootPoly returns the p-th root of f, assuming f is a perfect p-th power in GF(p)`x` (its only non-zero coefficients occur at exponents divisible… |
| `func RandomIrreducible(p *big.Int, n int, r *rand.Rand) (*Poly, error)` | RandomIrreducible returns a random monic irreducible polynomial of degree n over GF(p), drawing candidate coefficient vectors from r until one is… |
| `func (a *Poly) Add(b *Poly) *Poly` | Add returns the sum a + b over GF(p). |
| `func (a *Poly) Coefficient(i int) *big.Int` | Coefficient returns the coefficient of x^i, or 0 when i is out of range. |
| `func (a *Poly) Compose(b *Poly) *Poly` | Compose returns the composition a(b(x)) over GF(p). |
| `func (a *Poly) Copy() *Poly` | Copy returns an independent deep copy of the polynomial. |
| `func (a *Poly) Degree() int` | Degree returns the degree of the polynomial. |
| `func (a *Poly) Derivative() *Poly` | Derivative returns the formal derivative of the polynomial over GF(p). |
| `func (a *Poly) DivMod(b *Poly) (q, r *Poly, err error)` | DivMod returns the quotient q and remainder r of a divided by b, satisfying a = q*b + r with deg(r) < deg(b). |
| `func (a *Poly) Equal(b *Poly) bool` | Equal reports whether a and b are equal as polynomials over the same field. |
| `func (a *Poly) Eval(x *big.Int) *big.Int` | Eval evaluates the polynomial at the field element x using Horner's rule, returning a value in [0, p). |
| `func (a *Poly) ExtendedGcd(b *Poly) (g, s, t *Poly)` | ExtendedGcd returns g = gcd(a, b) together with cofactors s and t satisfying a*s + b*t = g. |
| `func (a *Poly) Gcd(b *Poly) *Poly` | Gcd returns the monic greatest common divisor of a and b over GF(p). |
| `func (a *Poly) IsConstant() bool` | IsConstant reports whether the polynomial has degree 0 or is zero. |
| `func (a *Poly) IsMonic() bool` | IsMonic reports whether the leading coefficient is 1. |
| `func (a *Poly) IsOne() bool` | IsOne reports whether the polynomial is the constant 1. |
| `func (a *Poly) IsZero() bool` | IsZero reports whether the polynomial is the zero polynomial. |
| `func (a *Poly) LeadingCoeff() *big.Int` | LeadingCoeff returns the coefficient of the highest-degree term, or 0 for the zero polynomial. |
| `func (a *Poly) Monic() *Poly` | Monic returns a scaled so that its leading coefficient is 1. |
| `func (a *Poly) Mul(b *Poly) *Poly` | Mul returns the product a * b over GF(p). |
| `func (a *Poly) Neg() *Poly` | Neg returns the additive inverse -a over GF(p). |
| `func (a *Poly) Pow(e int) *Poly` | Pow returns a raised to the non-negative integer power e. |
| `func (a *Poly) PowMod(e *big.Int, m *Poly) (*Poly, error)` | PowMod returns a^e reduced modulo the polynomial m, using square-and-multiply. |
| `func (a *Poly) Quo(b *Poly) (*Poly, error)` | Quo returns the quotient of a divided by b. |
| `func (a *Poly) Rem(b *Poly) (*Poly, error)` | Rem returns the remainder of a divided by b (a mod b). |
| `func (a *Poly) ScalarMul(s *big.Int) *Poly` | ScalarMul returns the polynomial scaled by the field element s. |
| `func (a *Poly) String() string` | String renders the polynomial in descending-degree infix form such as "x^2 + 3x + 1", or "0" for the zero polynomial. |
| `func (a *Poly) Sub(b *Poly) *Poly` | Sub returns the difference a - b over GF(p). |
| `func (a *Poly) ToIndex() *big.Int` | ToIndex encodes a polynomial's coefficient vector as a big integer for stable ordering: Σ Coeff`i`·p^i. |

</details>

<details>
<summary><code>PolyFactor</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Factor(f *Poly) (lead *big.Int, factors []PolyFactor)` | Factor returns the complete factorisation of f over GF(p): the leading coefficient together with the distinct monic irreducible factors and their… |
| `func SquareFreeFactorization(f *Poly) []PolyFactor` | SquareFreeFactorization returns the square-free decomposition of f: a set of pairwise-coprime square-free monic factors together with their… |

</details>

<details>
<summary><code>PrimePower</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FactorInt(n *big.Int) []PrimePower` | FactorInt returns the prime factorisation of \|n\| as a slice of PrimePower entries sorted by increasing prime. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/galois).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
