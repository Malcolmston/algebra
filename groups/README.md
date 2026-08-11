# groups — self-contained finite-group and modular-arithmetic toolkit built entirely on the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/groups.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/groups)

Package groups is a self-contained finite-group and modular-arithmetic toolkit
built entirely on the Go standard library.

It is a sibling subpackage of github.com/malcolmston/algebra but does not
depend on it: every routine works directly on machine integers, small
permutation slices, or explicit finite-group tables. All results are
deterministic and validated in the accompanying tests against closed-form
reference values.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/groups"
```

## Usage

This is the package's own `ExampleCompose`, so it compiles and its output is
asserted on every `go test ./groups/`.

```go
p := groups.PermFromCycles(3, [][]int{{0, 1}})
	q := groups.PermFromCycles(3, [][]int{{1, 2}})
	fmt.Println("p*q =", groups.Compose(p, q))
	fmt.Println("q*p =", groups.Compose(q, p))

	fmt.Println(groups.Compose(p, q).Equal(groups.Compose(q, p)))
```

```
p*q = (0 1 2)
q*p = (0 2 1)
false
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AlternatingGroupOrder(n int) int` | AlternatingGroupOrder returns \|A_n\|, the order of the alternating group on n points: n!/2 for n >= 2 and 1 for n in {0, 1}. |
| `func Coprime(a, b int) bool` | Coprime reports whether a and b are relatively prime, i.e. |
| `func CyclicOrder(n int) int` | CyclicOrder returns the order (number of elements) of the additive cyclic group Z/nZ, which is n. |
| `func DihedralElementOrder(n int, a DihedralElement) int` | DihedralElementOrder returns the order of a in D_n: 1 for the identity, 2 for every reflection, and n/gcd(k, n) for the rotation r^k. |
| `func DihedralOrder(n int) int` | DihedralOrder returns the order of the dihedral group D_n, which is 2n. |
| `func ElementOrderZn(a, n int) int` | ElementOrderZn returns the additive order of a in the cyclic group Z/nZ, that is the smallest k > 0 with k*a ≡ 0 (mod n). |
| `func EulerTotient(n int) int` | EulerTotient returns Euler's totient φ(n): the number of integers in [1, n] that are coprime to n, equivalently the order of the unit group (Z/nZ)*. |
| `func ExtendedGcd(a, b int) (g, x, y int)` | ExtendedGcd returns g = gcd(a, b) together with Bézout coefficients x and y satisfying a*x + b*y == g. |
| `func Factorial(n int) int` | Factorial returns n! as an int. |
| `func GFAdd(a, b, p int) int` | GFAdd returns a + b in the prime field GF(p). |
| `func GFDiv(a, b, p int) int` | GFDiv returns a / b in the prime field GF(p). |
| `func GFInv(a, p int) int` | GFInv returns the multiplicative inverse of a in the prime field GF(p). |
| `func GFMul(a, b, p int) int` | GFMul returns a * b in the prime field GF(p). |
| `func GFNeg(a, p int) int` | GFNeg returns the additive inverse -a in the prime field GF(p). |
| `func GFPow(a, e, p int) int` | GFPow returns a raised to the power e in the prime field GF(p). |
| `func GFSub(a, b, p int) int` | GFSub returns a - b in the prime field GF(p). |
| `func Gcd(a, b int) int` | Gcd returns the non-negative greatest common divisor of a and b using the Euclidean algorithm. |
| `func GcdMany(nums ...int) int` | GcdMany returns the greatest common divisor of all supplied integers. |
| `func IsPermutation(p Perm) bool` | IsPermutation reports whether p is a valid permutation: a bijection of {0, ..., len(p)-1} with every value appearing exactly once. |
| `func IsPrime(n int) bool` | IsPrime reports whether n is a prime number using trial division. |
| `func IsPrimitiveRoot(g, n int) bool` | IsPrimitiveRoot reports whether g is a primitive root modulo n, i.e. |
| `func IsUnitModN(a, n int) bool` | IsUnitModN reports whether a is a unit modulo n, i.e. |
| `func Lcm(a, b int) int` | Lcm returns the non-negative least common multiple of a and b. |
| `func LcmMany(nums ...int) int` | LcmMany returns the least common multiple of all supplied integers. |
| `func Mod(a, n int) int` | Mod returns the representative of a in Z/nZ lying in [0, n). |
| `func ModAdd(a, b, n int) int` | ModAdd returns (a + b) mod n in [0, n). |
| `func ModInverse(a, n int) (int, bool)` | ModInverse returns the multiplicative inverse of a modulo n in [0, n) and true when gcd(a, n) == 1. |
| `func ModMul(a, b, n int) int` | ModMul returns (a * b) mod n in [0, n) computed without intermediate overflow for operands whose reduced values fit in an int. |
| `func ModNeg(a, n int) int` | ModNeg returns (-a) mod n in [0, n). |
| `func ModPow(base, exp, n int) int` | ModPow returns base^exp mod n in [0, n) via binary exponentiation. |
| `func ModSub(a, b, n int) int` | ModSub returns (a - b) mod n in [0, n). |
| `func MultiplicativeOrder(a, n int) (int, bool)` | MultiplicativeOrder returns the multiplicative order of a modulo n, the smallest k > 0 with a^k ≡ 1 (mod n), together with true. |
| `func PolyEqual(a, b Poly) bool` | PolyEqual reports whether a and b are equal as polynomials, comparing coefficients within tolerance after trimming. |
| `func PolyEval(a Poly, x float64) float64` | PolyEval returns the value a(x) evaluated at the point x using Horner's method. |
| `func PolyModAdd(a, b []int, p int) []int` | PolyModAdd returns the sum a + b of two polynomials over GF(p), with coefficients reduced into [0, p). |
| `func PolyModDegree(a []int, p int) int` | PolyModDegree returns the degree of the polynomial a over GF(p): the index of its highest non-zero coefficient after reduction mod p, or -1 for the… |
| `func PolyModDivMod(a, b []int, p int) (q, r []int)` | PolyModDivMod returns the quotient q and remainder r of dividing a by b over GF(p), satisfying a = q·b + r with deg(r) < deg(b). |
| `func PolyModEval(a []int, x, p int) int` | PolyModEval returns the value a(x) mod p evaluated at x over GF(p) using Horner's method. |
| `func PolyModGCD(a, b []int, p int) []int` | PolyModGCD returns the monic greatest common divisor of a and b over GF(p) via the Euclidean algorithm. |
| `func PolyModMonic(a []int, p int) []int` | PolyModMonic returns a scaled so its leading coefficient is 1 over GF(p). |
| `func PolyModMul(a, b []int, p int) []int` | PolyModMul returns the product a·b of two polynomials over GF(p). |
| `func PolyModSub(a, b []int, p int) []int` | PolyModSub returns the difference a - b of two polynomials over GF(p). |
| `func PrimitiveRoots(n int) []int` | PrimitiveRoots returns every primitive root modulo n in [1, n) in increasing order. |
| `func SymmetricGroupOrder(n int) int` | SymmetricGroupOrder returns \|S_n\| = n!, the order of the symmetric group on n points. |
| `func UnitsModN(n int) []int` | UnitsModN returns the elements of the multiplicative group (Z/nZ)* in increasing order: every a in [1, n) with gcd(a, n) == 1. |

### Types

| Type | What it is |
| --- | --- |
| `DihedralElement` | DihedralElement represents an element of a dihedral group D_n, the symmetry group of a regular n-gon with 2n elements. |
| `FiniteGroup` | FiniteGroup is an abstract finite group whose elements are the integer indices 0, 1, ..., Order()-1 and whose binary operation is supplied as a… |
| `Perm` | Perm is a permutation of {0, 1, ..., n-1} in one-line notation: p`i` is the image of i under p. |
| `Poly` | Poly is a univariate polynomial over the real (rational) field, stored in coefficient slice form from lowest to highest degree: Poly{c0, c1, c2}… |

<details>
<summary><code>DihedralElement</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DihedralCompose(n int, a, b DihedralElement) DihedralElement` | DihedralCompose returns the product a·b in D_n using the defining relations r^n = s^2 = e and s·r = r^-1·s. |
| `func DihedralGroup(n int) []DihedralElement` | DihedralGroup returns all 2n elements of D_n: the n rotations r^0..r^(n-1) followed by the n reflections s·r^0..s·r^(n-1). |
| `func DihedralIdentity() DihedralElement` | DihedralIdentity returns the identity element of every dihedral group, the rotation by zero. |
| `func DihedralInverse(n int, a DihedralElement) DihedralElement` | DihedralInverse returns the inverse of a in D_n. |
| `func DihedralReflection(n, k int) DihedralElement` | DihedralReflection returns the reflection s·r^k in D_n with its exponent reduced modulo n. |
| `func DihedralRotation(n, k int) DihedralElement` | DihedralRotation returns the rotation r^k in D_n with its exponent reduced modulo n. |
| `func (a DihedralElement) String() string` | String renders a dihedral element as "r^k" for rotations and "s·r^k" for reflections, abbreviating the exponent-zero cases as "e", "r", "s". |

</details>

<details>
<summary><code>FiniteGroup</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CyclicGroupZn(n int) *FiniteGroup` | CyclicGroupZn returns the additive cyclic group Z/nZ as a FiniteGroup: the elements are 0..n-1, the operation is addition modulo n, and the identity… |
| `func NewFiniteGroup(n, identity int, op func(a, b int) int) *FiniteGroup` | NewFiniteGroup builds a FiniteGroup on the elements {0, ..., n-1} with the given operation and identity index. |
| `func SymmetricGroupTable(n int) *FiniteGroup` | SymmetricGroupTable returns the symmetric group S_n as a FiniteGroup whose element indices enumerate the permutations in the lexicographic order of… |
| `func (g *FiniteGroup) CayleyTable() [][]int` | CayleyTable returns the group's multiplication table as an n×n matrix where entry `a``b` holds the index of a·b. |
| `func (g *FiniteGroup) Center() []int` | Center returns the sorted indices of the center Z(G): the elements that commute with every element of the group. |
| `func (g *FiniteGroup) Cosets(subgroup []int) [][]int` | Cosets returns the partition of the group into distinct left cosets of the given subgroup, each coset sorted and the cosets ordered by their smallest… |
| `func (g *FiniteGroup) ElementOrder(a int) int` | ElementOrder returns the order of a: the smallest k > 0 with a^k equal to the identity. |
| `func (g *FiniteGroup) Elements() []int` | Elements returns the slice {0, 1, ..., Order()-1} of all element indices. |
| `func (g *FiniteGroup) GeneratedSubgroup(gens []int) []int` | GeneratedSubgroup returns the sorted set of elements of the subgroup generated by the given generators, computed by closing the generator set under… |
| `func (g *FiniteGroup) Identity() int` | Identity returns the index of the group's identity element. |
| `func (g *FiniteGroup) Index(subgroup []int) int` | Index returns the index [G:H] of the subgroup H in the group G, the number of distinct left cosets. |
| `func (g *FiniteGroup) Inverse(a int) int` | Inverse returns the index of the inverse of a. |
| `func (g *FiniteGroup) IsAbelian() bool` | IsAbelian reports whether the group is commutative: a·b == b·a for every pair of elements. |
| `func (g *FiniteGroup) IsCyclic() bool` | IsCyclic reports whether the group is cyclic, i.e. |
| `func (g *FiniteGroup) IsSubgroup(elems []int) bool` | IsSubgroup reports whether the given set of element indices forms a subgroup: it must be non-empty, contain the identity, and be closed under the… |
| `func (g *FiniteGroup) IsValid() bool` | IsValid reports whether the structure satisfies the group axioms over its n indices: closure, the declared identity acting as a two-sided identity,… |
| `func (g *FiniteGroup) LeftCoset(subgroup []int, x int) []int` | LeftCoset returns the sorted left coset g·H = {g·h : h in subgroup} for the element index x. |
| `func (g *FiniteGroup) Op(a, b int) int` | Op returns the product a·b of two element indices under the group operation. |
| `func (g *FiniteGroup) Order() int` | Order returns the number of elements in the group. |
| `func (g *FiniteGroup) RightCoset(subgroup []int, x int) []int` | RightCoset returns the sorted right coset H·g = {h·g : h in subgroup} for the element index x. |

</details>

<details>
<summary><code>Perm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AlternatingGroup(n int) []Perm` | AlternatingGroup returns all even permutations of degree n (the alternating group A_n) in lexicographic order. |
| `func Compose(p, q Perm) Perm` | Compose returns the composition p∘q, the permutation that applies q first and then p, so (Compose(p, q))`i` == p[q`i`]. |
| `func Identity(n int) Perm` | Identity returns the identity permutation of degree n, mapping every i to itself. |
| `func PermFromCycles(n int, cycles [][]int) Perm` | PermFromCycles builds a permutation of degree n from disjoint cycles given as slices of points. |
| `func SymmetricGroup(n int) []Perm` | SymmetricGroup returns all n! permutations of degree n in lexicographic order of their one-line notation. |
| `func Transposition(n, i, j int) Perm` | Transposition returns the permutation of degree n that swaps i and j and fixes every other point. |
| `func (p Perm) Apply(i int) int` | Apply returns the image p`i` of the point i under the permutation p. |
| `func (p Perm) CycleType() []int` | CycleType returns the multiset of cycle lengths of p sorted in non-increasing order (a partition of the degree n). |
| `func (p Perm) Cycles() [][]int` | Cycles returns the disjoint-cycle decomposition of p. |
| `func (p Perm) Degree() int` | Degree returns the size n of the set {0, ..., n-1} that p permutes, i.e. |
| `func (p Perm) Equal(q Perm) bool` | Equal reports whether p and q are the same permutation (equal degree and equal images). |
| `func (p Perm) Inverse() Perm` | Inverse returns the inverse permutation p^-1, satisfying Compose(p, p.Inverse()) == Identity(p.Degree()). |
| `func (p Perm) IsEven() bool` | IsEven reports whether p is an even permutation (sign +1). |
| `func (p Perm) IsIdentity() bool` | IsIdentity reports whether p fixes every point, i.e. |
| `func (p Perm) NumInversions() int` | NumInversions returns the number of inversions of p: pairs i < j with p`i` > p`j`. |
| `func (p Perm) Order() int` | Order returns the order of p in the symmetric group: the smallest k > 0 with p^k the identity. |
| `func (p Perm) Pow(k int) Perm` | Pow returns p raised to the integer power k. |
| `func (p Perm) Sign() int` | Sign returns the sign (parity) of p: +1 if p is an even permutation and -1 if it is odd. |
| `func (p Perm) String() string` | String renders p in disjoint-cycle notation, e.g. |

</details>

<details>
<summary><code>Poly</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PolyAdd(a, b Poly) Poly` | PolyAdd returns the sum a + b. |
| `func PolyDerivative(a Poly) Poly` | PolyDerivative returns the formal derivative a'(x) of the polynomial a. |
| `func PolyDivMod(a, b Poly) (q, r Poly)` | PolyDivMod returns the quotient q and remainder r of polynomial long division a = q·b + r with deg(r) < deg(b). |
| `func PolyGCD(a, b Poly) Poly` | PolyGCD returns the monic greatest common divisor of a and b via the Euclidean algorithm over polynomials. |
| `func PolyMod(a, b Poly) Poly` | PolyMod returns the remainder of a divided by b (the r from `PolyDivMod`). |
| `func PolyMonic(a Poly) Poly` | PolyMonic returns a scaled to be monic (leading coefficient 1). |
| `func PolyMul(a, b Poly) Poly` | PolyMul returns the product a·b via the convolution of coefficient slices. |
| `func PolyScale(a Poly, c float64) Poly` | PolyScale returns the polynomial a with every coefficient multiplied by the scalar c. |
| `func PolySub(a, b Poly) Poly` | PolySub returns the difference a - b. |
| `func PolyTrim(a Poly) Poly` | PolyTrim returns a copy of a with trailing (highest-degree) coefficients that are effectively zero removed, giving the canonical representation. |
| `func (a Poly) Degree() int` | Degree returns the degree of a: the index of its highest non-zero coefficient. |
| `func (a Poly) IsZero() bool` | IsZero reports whether a is the zero polynomial (every coefficient within tolerance of zero). |
| `func (a Poly) LeadingCoeff() float64` | LeadingCoeff returns the coefficient of the highest-degree term of a, or 0 for the zero polynomial. |
| `func (a Poly) String() string` | String renders a in conventional descending-degree notation, e.g. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/groups).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
