# ellipticcurves — implements the arithmetic and number theory of elliptic curves in short Weierstrass

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/ellipticcurves.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/ellipticcurves)

Package ellipticcurves implements the arithmetic and number theory of elliptic
curves in short Weierstrass form y^2 = x^3 + A*x + B, both over a prime finite
field F_p and over the field of rational numbers Q. The focus is mathematical
rather than cryptographic: the package is meant for exploring the group law,
invariants, torsion, point counting and pairings of small curves.

Over F_p the central types are CurveFp and its affine points PointFp. They
support the geometric group law (addition, doubling, negation and scalar
multiplication), the curve invariants discriminant and j-invariant, point
counting (naive enumeration, baby-step/giant-step point order and a
Mestre-style baby-step/giant-step curve order), division polynomials and the
Weil pairing.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/ellipticcurves"
```

## Usage

This is the package's own `ExampleCRT`, so it compiles and its output is
asserted on every `go test ./ellipticcurves/`.

```go
x, m, err := ec.CRT(
		[]*big.Int{big.NewInt(2), big.NewInt(3), big.NewInt(2)},
		[]*big.Int{big.NewInt(3), big.NewInt(5), big.NewInt(7)},
	)
	fmt.Println(x, m, err)
	_, _, err = ec.CRT([]*big.Int{big.NewInt(1)}, []*big.Int{big.NewInt(0)})
	fmt.Println(err)
```

```
23 105 <nil>
ellipticcurves: inconsistent CRT system
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func CRT(residues, moduli []*big.Int) (*big.Int, *big.Int, error)` | CRT solves the system x = residues`i` (mod moduli`i`) for pairwise coprime moduli, returning the least non-negative solution and the product of the… |
| `func Divisors(n *big.Int) []*big.Int` | Divisors returns all positive divisors of \|n\| in ascending order. |
| `func EulerTotient(n *big.Int) *big.Int` | EulerTotient returns Euler's totient phi(\|n\|), the count of integers in [1, n] coprime to n. |
| `func Factorize(n *big.Int) map[*big.Int]int` | Factorize returns the prime factorization of \|n\| as a map from prime to exponent. |
| `func Gcd(a, b *big.Int) *big.Int` | Gcd returns the greatest common divisor of a and b as a non-negative integer. |
| `func IntSqrt(n *big.Int) *big.Int` | IntSqrt returns the integer square root floor(sqrt(n)) for n >= 0. |
| `func IsPerfectSquare(n *big.Int) (*big.Int, bool)` | IsPerfectSquare reports whether n is a perfect square, returning the square root when it is. |
| `func IsPrimitiveRoot(g, p *big.Int) bool` | IsPrimitiveRoot reports whether g is a primitive root modulo the odd prime p, i.e. |
| `func IsProbablePrime(n *big.Int, rounds int) bool` | IsProbablePrime reports whether n is prime according to the Miller-Rabin and Baillie-PSW tests in big.Int.ProbablyPrime with the given number of… |
| `func IsQuadraticResidue(a, p *big.Int) bool` | IsQuadraticResidue reports whether a is a non-zero quadratic residue modulo the odd prime p. |
| `func JacobiSymbol(a, n *big.Int) int` | JacobiSymbol returns the Jacobi symbol (a/n) for a positive odd integer n. |
| `func KroneckerSymbol(a, n *big.Int) int` | KroneckerSymbol returns the Kronecker symbol (a/n), the total multiplicative extension of the Jacobi and Legendre symbols to all integers n, taking… |
| `func Lcm(a, b *big.Int) *big.Int` | Lcm returns the least common multiple of a and b as a non-negative integer. |
| `func LcmSlice(xs []*big.Int) *big.Int` | LcmSlice returns the least common multiple of a slice of integers, or 1 for an empty slice. |
| `func LegendreSymbol(a, p *big.Int) int` | LegendreSymbol returns the Legendre symbol (a/p) for an odd prime p, one of -1, 0 or 1. |
| `func Mod(a, m *big.Int) *big.Int` | Mod returns the canonical representative of a modulo m in the range [0, m). |
| `func ModAdd(a, b, m *big.Int) *big.Int` | ModAdd returns (a + b) mod m reduced into [0, m). |
| `func ModDiv(a, b, m *big.Int) (*big.Int, error)` | ModDiv returns (a * b^-1) mod m, or ErrNotInvertible when b is not invertible modulo m. |
| `func ModDouble(a, m *big.Int) *big.Int` | ModDouble returns (2*a) mod m reduced into [0, m). |
| `func ModExp(a, e, m *big.Int) (*big.Int, error)` | ModExp returns a^e mod m using big.Int.Exp. |
| `func ModInverse(a, m *big.Int) (*big.Int, error)` | ModInverse returns the multiplicative inverse of a modulo m, or ErrNotInvertible when gcd(a, m) != 1. |
| `func ModMul(a, b, m *big.Int) *big.Int` | ModMul returns (a * b) mod m reduced into [0, m). |
| `func ModNeg(a, m *big.Int) *big.Int` | ModNeg returns (-a) mod m reduced into [0, m). |
| `func ModPow(a, e, m *big.Int) *big.Int` | ModPow returns a^e mod m for a non-negative exponent e, a convenience wrapper over big.Int.Exp that reduces a into [0, m) first. |
| `func ModSqrt(a, p *big.Int) (*big.Int, error)` | ModSqrt returns a square root of a modulo the odd prime p, or ErrNoSquareRoot when a is a quadratic non-residue. |
| `func ModSqrtTonelli(a, p *big.Int) (*big.Int, error)` | ModSqrtTonelli returns a square root of a modulo the odd prime p using an explicit Tonelli-Shanks implementation. |
| `func ModSquare(a, m *big.Int) *big.Int` | ModSquare returns (a*a) mod m reduced into [0, m). |
| `func ModSub(a, b, m *big.Int) *big.Int` | ModSub returns (a - b) mod m reduced into [0, m). |
| `func MultiplicativeOrder(a, p *big.Int) *big.Int` | MultiplicativeOrder returns the multiplicative order of a modulo the odd prime p, that is the least k > 0 with a^k = 1 mod p. |
| `func NextPrime(n *big.Int) *big.Int` | NextPrime returns the least prime strictly greater than n. |
| `func NumDivisors(n *big.Int) int` | NumDivisors returns the number of positive divisors of \|n\|. |
| `func OrderModN(a, n *big.Int) *big.Int` | OrderModN returns the multiplicative order of a modulo n for a general modulus n > 1 with gcd(a, n) = 1, computed by iterated multiplication. |
| `func PollardRho(n *big.Int) *big.Int` | PollardRho returns a non-trivial factor of the composite odd integer n using Pollard's rho algorithm, or nil on failure. |
| `func PrimeFactors(n *big.Int) []*big.Int` | PrimeFactors returns the distinct prime factors of \|n\| in ascending order. |
| `func PrimitiveRoot(p *big.Int) *big.Int` | PrimitiveRoot returns a generator of the multiplicative group modulo the odd prime p, i.e. |
| `func RandomFieldElement(p *big.Int, rng *rand.Rand) *big.Int` | RandomFieldElement returns a uniformly random element of F_p using the supplied source. |
| `func RandomNonzeroFieldElement(p *big.Int, rng *rand.Rand) *big.Int` | RandomNonzeroFieldElement returns a uniformly random element of F_p^* using the supplied source. |
| `func SqrtCandidates(a, p *big.Int) (*big.Int, *big.Int, error)` | SqrtCandidates returns both square roots of a modulo the odd prime p as (r, p-r). |
| `func TrialDivision(n, limit *big.Int) *big.Int` | TrialDivision returns the smallest prime factor of n found by trial division up to the given limit, or nil if none is found in range. |

### Types

| Type | What it is |
| --- | --- |
| `CurveFp` | CurveFp is an elliptic curve in short Weierstrass form y^2 = x^3 + A*x + B over the prime field F_p, with p > 3. |
| `CurveQ` | CurveQ is an elliptic curve y^2 = x^3 + A*x + B over the rationals Q, with A and B stored as reduced big.Rat values. |
| `PointFp` | PointFp is an affine point on a CurveFp, or the point at infinity when Infinity is set. |
| `PointQ` | PointQ is an affine rational point on a CurveQ, or the point at infinity when Infinity is set. |

<details>
<summary><code>CurveFp</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCurveFp(a, b, p *big.Int) (*CurveFp, error)` | NewCurveFp constructs the curve y^2 = x^3 + A*x + B over F_p. |
| `func (c *CurveFp) Add(p, q PointFp) PointFp` | Add returns the group sum p + q using the affine chord-and-tangent law. |
| `func (c *CurveFp) C4() *big.Int` | C4 returns the invariant c4 = -48*A reduced modulo p. |
| `func (c *CurveFp) C6() *big.Int` | C6 returns the invariant c6 = -864*B reduced modulo p. |
| `func (c *CurveFp) Cofactor(n *big.Int) *big.Int` | Cofactor returns #E(F_p) / n, the index of the order-n subgroup, or nil when n does not divide the group order. |
| `func (c *CurveFp) Contains(pt PointFp) bool` | Contains reports whether pt lies on the curve; it is an alias for IsOnCurve. |
| `func (c *CurveFp) CountAffinePoints() *big.Int` | CountAffinePoints returns the number of affine points on the curve by direct enumeration: for each x in F_p it adds 1 + (rhs/p) using the Legendre… |
| `func (c *CurveFp) CurveOrderBSGS(rng *rand.Rand) (*big.Int, error)` | CurveOrderBSGS returns #E(F_p) using a Mestre-style baby-step/giant-step method: it takes random points, uses their point orders to constrain the… |
| `func (c *CurveFp) CurveOrderNaive() *big.Int` | CurveOrderNaive returns the exact order of the group E(F_p), namely the number of affine points plus one for the point at infinity. |
| `func (c *CurveFp) DiscreteLogBSGS(p, q PointFp, order *big.Int) (*big.Int, error)` | DiscreteLogBSGS returns x in [0, order) with x*p = q, using baby-step/giant-step, given the order of the base point p. |
| `func (c *CurveFp) Discriminant() *big.Int` | Discriminant returns the discriminant -16*(4*A^3 + 27*B^2) reduced modulo p. |
| `func (c *CurveFp) DivisionPolynomial2(x, y *big.Int) *big.Int` | DivisionPolynomial2 returns psi_2(x, y) = 2*y evaluated modulo p. |
| `func (c *CurveFp) DivisionPolynomial3(x, y *big.Int) *big.Int` | DivisionPolynomial3 returns psi_3(x, y) = 3*x^4 + 6*A*x^2 + 12*B*x - A^2 evaluated modulo p. |
| `func (c *CurveFp) DivisionPolynomial4(x, y *big.Int) *big.Int` | DivisionPolynomial4 returns psi_4(x, y) modulo p, equal to 4*y*(x^6 + 5*A*x^4 + 20*B*x^3 - 5*A^2*x^2 - 4*A*B*x - 8*B^2 - A^3). |
| `func (c *CurveFp) DivisionPolynomialSequence(n int, x, y *big.Int) ([]*big.Int, error)` | DivisionPolynomialSequence returns psi_0, psi_1, ..., psi_n evaluated at (x, y) modulo p. |
| `func (c *CurveFp) DivisionPolynomialValue(n int, x, y *big.Int) (*big.Int, error)` | DivisionPolynomialValue returns psi_n(x, y) evaluated at the point (x, y) modulo p, using the standard recursion. |
| `func (c *CurveFp) Double(pt PointFp) PointFp` | Double returns 2*pt using the tangent-line slope. |
| `func (c *CurveFp) EmbeddingDegree(n *big.Int) *big.Int` | EmbeddingDegree returns the embedding degree of the order-n subgroup, the least k >= 1 with n dividing p^k - 1. |
| `func (c *CurveFp) EnumeratePoints() []PointFp` | EnumeratePoints returns every point of E(F_p) including the point at infinity. |
| `func (c *CurveFp) Equal(other *CurveFp) bool` | Equal reports whether two curves have identical coefficients and modulus. |
| `func (c *CurveFp) FindPointOfOrder(n *big.Int, rng *rand.Rand) (PointFp, error)` | FindPointOfOrder returns a point of exact order n, or ErrOrderNotFound if no such point exists (n must divide the group order and divide the… |
| `func (c *CurveFp) FrobeniusCharacteristicPoly() []*big.Int` | FrobeniusCharacteristicPoly returns the coefficients (1, -a_p, p) of the characteristic polynomial T^2 - a_p*T + p of the Frobenius endomorphism, as… |
| `func (c *CurveFp) GroupExponent(rng *rand.Rand) (*big.Int, error)` | GroupExponent returns the exponent of E(F_p), the least common multiple of the orders of all points, computed from the group order and structure. |
| `func (c *CurveFp) GroupStructureFp(rng *rand.Rand) (*big.Int, *big.Int)` | GroupStructureFp returns the invariant factors (n1, n2) with n1 dividing n2 and n1*n2 = #E(F_p), so that E(F_p) is isomorphic to Z/n1 x Z/n2. |
| `func (c *CurveFp) HasFullTwoTorsion() bool` | HasFullTwoTorsion reports whether all three non-trivial 2-torsion points are rational over F_p, equivalently the cubic x^3 + A*x + B splits… |
| `func (c *CurveFp) HasseInterval() (*big.Int, *big.Int)` | HasseInterval returns the closed interval [p+1-2*sqrt(p), p+1+2*sqrt(p)] in which #E(F_p) must lie, as (lo, hi). |
| `func (c *CurveFp) IsAnomalous() bool` | IsAnomalous reports whether #E(F_p) = p, equivalently a_p = 1. |
| `func (c *CurveFp) IsCyclicFp(rng *rand.Rand) bool` | IsCyclicFp reports whether E(F_p) is cyclic, tested by comparing the group exponent against the group order. |
| `func (c *CurveFp) IsInfinity(pt PointFp) bool` | IsInfinity reports whether pt is the point at infinity. |
| `func (c *CurveFp) IsIsomorphicFp(other *CurveFp) bool` | IsIsomorphicFp reports whether two curves over the same field are isomorphic over F_p, i.e. |
| `func (c *CurveFp) IsNTorsionByDivisionPoly(n int, pt PointFp) (bool, error)` | IsNTorsionByDivisionPoly reports whether pt is an n-torsion point (n*pt = O) by testing whether psi_n vanishes at pt. |
| `func (c *CurveFp) IsNTorsionPoint(n *big.Int, pt PointFp) bool` | IsNTorsionPoint reports whether n*pt equals the point at infinity, computed directly with the group law. |
| `func (c *CurveFp) IsOnCurve(pt PointFp) bool` | IsOnCurve reports whether pt satisfies y^2 = x^3 + A*x + B mod p. |
| `func (c *CurveFp) IsOrdinary() bool` | IsOrdinary reports whether the curve is ordinary over F_p, the negation of IsSupersingular. |
| `func (c *CurveFp) IsQuadraticTwistFp(other *CurveFp) bool` | IsQuadraticTwistFp reports whether other is a quadratic twist of the receiver over F_p: the curves share a j-invariant but are not isomorphic over… |
| `func (c *CurveFp) IsSingular() bool` | IsSingular reports whether the curve discriminant vanishes modulo p. |
| `func (c *CurveFp) IsSupersingular() bool` | IsSupersingular reports whether the curve is supersingular over F_p, i.e. |
| `func (c *CurveFp) IsTwoTorsion(pt PointFp) bool` | IsTwoTorsion reports whether pt is a non-identity point of order two, i.e. |
| `func (c *CurveFp) IsomorphismScaleFp(other *CurveFp) (*big.Int, bool)` | IsomorphismScaleFp returns a field element u with A' = u^4*A and B' = u^6*B witnessing an isomorphism from the receiver to other, and true when such… |
| `func (c *CurveFp) JInvariant() (*big.Int, error)` | JInvariant returns the j-invariant 1728 * 4*A^3 / (4*A^3 + 27*B^2) reduced modulo p. |
| `func (c *CurveFp) LiftX(x *big.Int) (PointFp, PointFp, error)` | LiftX returns the up-to-two affine points whose x-coordinate equals x. |
| `func (c *CurveFp) MillerFunction(n *big.Int, p, z PointFp) (*big.Int, error)` | MillerFunction evaluates the Miller function f_{n,P} at the affine point Z, where div(f_{n,P}) = n*(P) - n*(O) when n*P = O. Z must avoid the support… |
| `func (c *CurveFp) Multiply(k *big.Int, pt PointFp) PointFp` | Multiply is an alias for ScalarMul. |
| `func (c *CurveFp) Neg(pt PointFp) PointFp` | Neg returns the additive inverse (x, -y) of pt. |
| `func (c *CurveFp) Negate(pt PointFp) PointFp` | Negate is an alias for Neg. |
| `func (c *CurveFp) NewPointFp(x, y *big.Int) (PointFp, error)` | NewPointFp constructs an affine point (x, y) on the curve, reducing the coordinates modulo p. |
| `func (c *CurveFp) NumberOfPoints() *big.Int` | NumberOfPoints returns #E(F_p) including the point at infinity; it is an alias for CurveOrderNaive. |
| `func (c *CurveFp) Order() *big.Int` | Order returns the order of the group E(F_p); it is an alias for CurveOrderNaive. |
| `func (c *CurveFp) PointEqual(p, q PointFp) bool` | PointEqual reports whether two points are equal, treating all points at infinity as equal. |
| `func (c *CurveFp) PointOrder(pt PointFp, groupOrder *big.Int) *big.Int` | PointOrder returns the order of pt computed from the group order: it factors #E(F_p) and returns the least divisor n with n*pt = O. This is exact and… |
| `func (c *CurveFp) PointOrderBSGS(pt PointFp) (*big.Int, error)` | PointOrderBSGS returns the order of pt, the least m > 0 with m*pt = O, using baby-step/giant-step over the Hasse interval followed by prime-factor… |
| `func (c *CurveFp) QuadraticTwistFp(d *big.Int) (*CurveFp, error)` | QuadraticTwistFp returns the quadratic twist of the curve by a non-residue d: the curve y^2 = x^3 + A*d^2*x + B*d^3. |
| `func (c *CurveFp) RandomPointFp(rng *rand.Rand) PointFp` | RandomPointFp returns a random affine point on the curve using the supplied source. |
| `func (c *CurveFp) RightHandSide(x *big.Int) *big.Int` | RightHandSide returns x^3 + A*x + B mod p, the right-hand side of the curve equation evaluated at x. |
| `func (c *CurveFp) SameJInvariantFp(other *CurveFp) bool` | SameJInvariantFp reports whether two curves over the same field share a j-invariant, a necessary condition for being isomorphic or twists. |
| `func (c *CurveFp) ScalarMul(k *big.Int, pt PointFp) PointFp` | ScalarMul returns k*pt using left-to-right double-and-add. |
| `func (c *CurveFp) ScaleCurveFp(u *big.Int) (*CurveFp, error)` | ScaleCurveFp returns the curve y^2 = x^3 + u^4*A*x + u^6*B over the same field, the image of the receiver under the scaling isomorphism with… |
| `func (c *CurveFp) String() string` | String renders the curve equation and field. |
| `func (c *CurveFp) Sub(p, q PointFp) PointFp` | Sub returns p - q. |
| `func (c *CurveFp) SubgroupPoints(gen PointFp) []PointFp` | SubgroupPoints returns the cyclic subgroup generated by gen as a slice of points, starting from the point at infinity and ending just before… |
| `func (c *CurveFp) TorsionBasisFp(n *big.Int, rng *rand.Rand) (PointFp, PointFp, bool)` | TorsionBasisFp returns two points P, Q generating the full n-torsion E`n` = Z/n x Z/n over F_p, when it is rational, together with true. |
| `func (c *CurveFp) TraceOfFrobenius() *big.Int` | TraceOfFrobenius returns a_p = p + 1 - #E(F_p) computed from the naive order. |
| `func (c *CurveFp) TransformPointFp(u *big.Int, pt PointFp) PointFp` | TransformPointFp maps a point under the isomorphism (x, y) -> (u^2*x, u^3*y) that carries the curve y^2 = x^3 + A*x + B to y^2 = x^3 + u^4*A*x +… |
| `func (c *CurveFp) TwistOrder() *big.Int` | TwistOrder returns the order of the quadratic twist of the curve, which is 2*(p+1) - #E(F_p) = p + 1 + a_p. |
| `func (c *CurveFp) TwoTorsionPoints() []PointFp` | TwoTorsionPoints returns the non-identity 2-torsion points, the affine points with y = 0, i.e. |
| `func (c *CurveFp) WeilPairing(n *big.Int, pPt, qPt PointFp, rng *rand.Rand) (*big.Int, error)` | WeilPairing returns the Weil pairing e_n(P, Q) of two n-torsion points as an element of F_p, using Miller's algorithm with random divisor shifts… |
| `func (c *CurveFp) WeilPairingRootOfUnity(v, n *big.Int) bool` | WeilPairingRootOfUnity reports whether v is an n-th root of unity in F_p, a sanity check on a Weil pairing value: v^n = 1. |

</details>

<details>
<summary><code>CurveQ</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewCurveQ(a, b *big.Rat) (*CurveQ, error)` | NewCurveQ constructs y^2 = x^3 + A*x + B over Q. It returns ErrSingularCurveQ when the discriminant vanishes. |
| `func NewCurveQInt(a, b *big.Int) (*CurveQ, error)` | NewCurveQInt constructs y^2 = x^3 + A*x + B over Q from integer coefficients. |
| `func (c *CurveQ) APAtPrime(p *big.Int) (*big.Int, error)` | APAtPrime returns the trace of Frobenius a_p = p + 1 - #E~(F_p) of the good reduction of an integral curve at the prime p. |
| `func (c *CurveQ) Add(p, q PointQ) PointQ` | Add returns p + q using the rational chord-and-tangent law. |
| `func (c *CurveQ) AreIndependent(points []PointQ, iterations int, tol float64) bool` | AreIndependent reports whether the given points are linearly independent in E(Q) modulo torsion, judged by the height-pairing regulator exceeding tol… |
| `func (c *CurveQ) BadPrimes() ([]*big.Int, error)` | BadPrimes returns the primes of bad reduction of an integral curve, namely the prime divisors of its discriminant. |
| `func (c *CurveQ) CanonicalHeightApprox(pt PointQ, iterations int) float64` | CanonicalHeightApprox returns an approximation to the Neron-Tate canonical height of pt, computed as the limit 4^{-n} * h(x(2^n * pt)) of the naive… |
| `func (c *CurveQ) Contains(pt PointQ) bool` | Contains reports whether pt lies on the rational curve; it is an alias for IsOnCurve. |
| `func (c *CurveQ) Discriminant() *big.Rat` | Discriminant returns the discriminant -16*(4*A^3 + 27*B^2) as a rational. |
| `func (c *CurveQ) DiscriminantInt() (*big.Int, bool)` | DiscriminantInt returns the discriminant of an integral curve as an integer, together with true. |
| `func (c *CurveQ) Double(pt PointQ) PointQ` | Double returns 2*pt using the tangent slope. |
| `func (c *CurveQ) HasGoodReduction(p *big.Int) bool` | HasGoodReduction reports whether an integral curve has good reduction at the prime p, i.e. |
| `func (c *CurveQ) HeightPairing(p, q PointQ, iterations int) float64` | HeightPairing returns the Neron-Tate height pairing <P, Q> = (h(P+Q) - h(P) - h(Q)) / 2 using canonical heights approximated with the given number of… |
| `func (c *CurveQ) HeightPairingMatrix(points []PointQ, iterations int) [][]float64` | HeightPairingMatrix returns the Gram matrix of the height pairing of the given points, an n-by-n symmetric matrix of canonical height pairings. |
| `func (c *CurveQ) IntCoefficients() (*big.Int, *big.Int, bool)` | IntCoefficients returns the integer coefficients (A, B) of an integral curve and true, or (nil, nil, false) when the coefficients are not integers. |
| `func (c *CurveQ) IsInfinity(pt PointQ) bool` | IsInfinity reports whether pt is the rational point at infinity. |
| `func (c *CurveQ) IsIntegralCurve() bool` | IsIntegralCurve reports whether A and B are integers, a precondition for the Nagell-Lutz theorem. |
| `func (c *CurveQ) IsIntegralPoint(pt PointQ) bool` | IsIntegralPoint reports whether both coordinates of pt are integers. |
| `func (c *CurveQ) IsIsomorphicQ(other *CurveQ) bool` | IsIsomorphicQ reports whether two rational curves are isomorphic over Q. |
| `func (c *CurveQ) IsOnCurve(pt PointQ) bool` | IsOnCurve reports whether pt satisfies the curve equation. |
| `func (c *CurveQ) IsQuadraticTwistQ(other *CurveQ) bool` | IsQuadraticTwistQ reports whether other is a quadratic twist of the receiver over Q, i.e. |
| `func (c *CurveQ) IsTorsionPoint(pt PointQ) bool` | IsTorsionPoint reports whether pt has finite order in E(Q). |
| `func (c *CurveQ) IsomorphismScaleQ(other *CurveQ) (*big.Rat, bool)` | IsomorphismScaleQ returns a rational u with A' = u^4*A and B' = u^6*B witnessing an isomorphism from the receiver to other over Q, and true when such… |
| `func (c *CurveQ) JInvariant() (*big.Rat, error)` | JInvariant returns the j-invariant 1728 * 4*A^3 / (4*A^3 + 27*B^2) as a rational. |
| `func (c *CurveQ) Multiply(k *big.Int, pt PointQ) PointQ` | Multiply is an alias for ScalarMul over the rationals. |
| `func (c *CurveQ) NagellLutzCandidates() ([]PointQ, error)` | NagellLutzCandidates returns the integral points (x, y) satisfying the Nagell-Lutz necessary conditions for a torsion point: x and y integers with… |
| `func (c *CurveQ) NaiveHeightQ(pt PointQ) float64` | NaiveHeightQ returns the naive logarithmic height of the point pt, defined as log(max(\|num\|, \|den\|)) of its x-coordinate in lowest terms. |
| `func (c *CurveQ) Neg(pt PointQ) PointQ` | Neg returns the additive inverse (x, -y) of pt. |
| `func (c *CurveQ) Negate(pt PointQ) PointQ` | Negate is an alias for Neg over the rationals. |
| `func (c *CurveQ) NewPointQ(x, y *big.Rat) (PointQ, error)` | NewPointQ constructs the affine point (x, y) and checks that it lies on the curve, returning ErrPointNotOnCurveQ otherwise. |
| `func (c *CurveQ) PointEqual(p, q PointQ) bool` | PointEqual reports whether two rational points are equal, treating all points at infinity as equal. |
| `func (c *CurveQ) PointOrderQ(pt PointQ) int` | PointOrderQ returns the order of pt in E(Q), the least k > 0 with k*pt = O, or 0 when pt has infinite order (no such k up to Mazur's bound of 12). |
| `func (c *CurveQ) QuadraticTwistQ(d *big.Rat) (*CurveQ, error)` | QuadraticTwistQ returns the quadratic twist of the curve by the rational d: the curve y^2 = x^3 + A*d^2*x + B*d^3. |
| `func (c *CurveQ) RankLowerBound(points []PointQ, iterations int, tol float64) int` | RankLowerBound returns a heuristic lower bound on the Mordell-Weil rank of the curve given a list of rational points: the numerical rank of the… |
| `func (c *CurveQ) ReduceModP(p *big.Int) (*CurveFp, error)` | ReduceModP returns the reduction of the rational curve modulo the prime p as a CurveFp. |
| `func (c *CurveQ) ReducePointModP(p *big.Int, pt PointQ) (PointFp, error)` | ReducePointModP reduces a rational point of the curve modulo the prime p to a point on the reduced curve, returning ErrBadReduction when a coordinate… |
| `func (c *CurveQ) Regulator(points []PointQ, iterations int) float64` | Regulator returns the determinant of the height-pairing Gram matrix of the given points, an approximation to the elliptic regulator. |
| `func (c *CurveQ) RightHandSide(x *big.Rat) *big.Rat` | RightHandSide returns x^3 + A*x + B evaluated at x as a rational. |
| `func (c *CurveQ) SameJInvariantQ(other *CurveQ) bool` | SameJInvariantQ reports whether two rational curves share a j-invariant. |
| `func (c *CurveQ) ScalarMul(k *big.Int, pt PointQ) PointQ` | ScalarMul returns k*pt using double-and-add. |
| `func (c *CurveQ) String() string` | String renders the rational curve equation. |
| `func (c *CurveQ) Sub(p, q PointQ) PointQ` | Sub returns p - q. |
| `func (c *CurveQ) TorsionExponent() (int, error)` | TorsionExponent returns the exponent of the torsion subgroup, the least common multiple of the orders of its elements. |
| `func (c *CurveQ) TorsionOrder() (int, error)` | TorsionOrder returns the order of the torsion subgroup of E(Q). |
| `func (c *CurveQ) TorsionSubgroup() ([]PointQ, error)` | TorsionSubgroup returns all torsion points of E(Q), including the point at infinity, obtained by filtering the Nagell-Lutz candidates for finite… |
| `func (c *CurveQ) TwoTorsionPointsQ() ([]PointQ, error)` | TwoTorsionPointsQ returns the rational 2-torsion points (the affine points with y = 0, i.e. |

</details>

<details>
<summary><code>PointFp</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PointAtInfinityFp() PointFp` | PointAtInfinityFp returns the identity element of the group, the point at infinity. |
| `func (pt PointFp) String() string` | String renders a point in affine coordinates, or "O" for infinity. |

</details>

<details>
<summary><code>PointQ</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PointAtInfinityQ() PointQ` | PointAtInfinityQ returns the identity element of E(Q). |
| `func (pt PointQ) String() string` | String renders a rational point, or "O" for infinity. |

</details>

### Variables

`ErrNotPrime`, `ErrNoSquareRoot`, `ErrNotInvertible`, `ErrZeroModulus`, `ErrAmbiguousOrder`, `ErrBadReduction`, `ErrCRTInconsistent`, `ErrDivisionPolyAtTwoTorsion`, `ErrNonIntegralCurve`, `ErrOrderNotFound`, `ErrPairingDegenerate`, `ErrPairingInput`, `ErrPointNotOnCurve`, `ErrPointNotOnCurveQ`, `ErrSingularCurve`, `ErrSingularCurveQ`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/ellipticcurves).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
