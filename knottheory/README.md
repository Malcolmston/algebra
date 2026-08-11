# knottheory — implements combinatorial knot and braid theory using only the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/knottheory.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/knottheory)

Package knottheory implements combinatorial knot and braid theory using only
the Go standard library.

- Laurent polynomials over the integers (the Laurent type), the coefficient
ring in which every polynomial knot invariant lives, with a full set of ring
operations, evaluation, substitution and normalisation helpers. - Permutations
(the Permutation type): composition, inversion, cycle decomposition, order and
sign, used to describe the underlying permutation and closure of a braid. -
Braid groups (the Braid type): Artin generators and relations, free reduction,
exponent sum (writhe), the induced permutation, mirror and reverse, torus
braids, the full twist and the Garside half twist, and the reduced Burau
representation over the Laurent ring. - Knot and link diagrams described by
signed Gauss codes (the GaussCode and Diagram types): validation, crossing
number, writhe, self-writhe and the Gauss linking number between components. -
Planar-diagram (PD) codes (the PDCode type): the Kauffman bracket via a
union-find state sum, the writhe-normalised bracket, and the Jones polynomial
for small diagrams, both in the standard variable t and in the variable
t^{1/2} for links. - The Alexander polynomial, computed two independent ways:
from the reduced Burau matrix of a braid, and directly from a diagram via the
Alexander (Fox-derivative) matrix, together with the knot determinant and
p-colourability counts. - Reidemeister-move recognisers (types I, II and III)
that operate combinatorially on Gauss codes. - Torus-knot invariants T(p, q):
crossing number, Seifert genus, bridge and unknotting numbers, determinant,
and the closed forms of the Alexander and Jones polynomials.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/knottheory"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BraidGroupRank(n int) int` | BraidGroupRank returns the number of Artin generators of B_n, namely n-1. |
| `func IsTorusLink(p, q int) bool` | IsTorusLink reports whether T(p, q) is a genuine link of more than one component (gcd(\|p\|,\|q\|) > 1). |
| `func IsTrivialTorusKnot(p, q int) bool` | IsTrivialTorusKnot reports whether T(p, q) is the unknot, which happens iff \|p\| <= 1 or \|q\| <= 1 (and the pair is coprime). |
| `func IsValidBraidWord(strands int, word []int) bool` | IsValidBraidWord reports whether word is a legal word for B_n, that is every entry is non-zero and has absolute value at most n-1. |
| `func TorusIsKnot(p, q int) bool` | TorusIsKnot reports whether T(p, q) is a knot, that is gcd(\|p\|,\|q\|) = 1. |
| `func TorusKnotArf(p, q int) int` | TorusKnotArf returns the Arf invariant of the torus knot T(p, q) with coprime p, q, derived from its determinant. |
| `func TorusKnotBraidWord(p, q int) []int` | TorusKnotBraidWord returns the generator word of the torus braid whose closure is T(p, q). |
| `func TorusKnotBridgeNumber(p, q int) int` | TorusKnotBridgeNumber returns the bridge number of T(p, q) with coprime p, q >= 2, which is min(\|p\|, \|q\|). |
| `func TorusKnotCrossingNumber(p, q int) int` | TorusKnotCrossingNumber returns the crossing number of the torus knot T(p, q) with coprime p, q >= 2, namely min(p(q-1), q(p-1)). |
| `func TorusKnotDeterminant(p, q int) int` | TorusKnotDeterminant returns the determinant of the torus knot T(p, q), the absolute value of its Alexander polynomial at t = -1. |
| `func TorusKnotGenus(p, q int) int` | TorusKnotGenus returns the Seifert genus of the torus knot T(p, q) with coprime p, q, which is (p-1)(q-1)/2. |
| `func TorusKnotSeifertGenus(p, q int) int` | TorusKnotSeifertGenus is an alias for TorusKnotGenus. |
| `func TorusKnotUnknottingNumber(p, q int) int` | TorusKnotUnknottingNumber returns the unknotting number of the torus knot T(p, q) with coprime p, q, equal to its genus (p-1)(q-1)/2 by the… |
| `func TorusLinkComponents(p, q int) int` | TorusLinkComponents returns the number of components of T(p, q), which is gcd(\|p\|, \|q\|). |

### Types

| Type | What it is |
| --- | --- |
| `Braid` | Braid is an element of the Artin braid group B_n on n strands, represented by a word in the standard generators. |
| `Diagram` | Diagram is a knot or link diagram given as one signed Gauss code per component. |
| `GaussCode` | GaussCode is the signed Gauss code of a single knot component: the cyclic sequence of crossing encounters met while traversing the knot once. |
| `GaussEntry` | GaussEntry is a single crossing encounter recorded while traversing a knot or link component. |
| `Laurent` | Laurent is a Laurent polynomial in a single variable with integer coefficients, that is a finite sum of terms c*X^e where e may be any integer… |
| `LaurentMatrix` | LaurentMatrix is a dense matrix whose entries are Laurent polynomials. |
| `LaurentTerm` | LaurentTerm is a single monomial c*X^e of a Laurent polynomial. |
| `PDCode` | PDCode is a knot or link diagram in planar-diagram notation together with a count of any crossingless unknotted circles. |
| `PDCrossing` | PDCrossing is one crossing of a planar-diagram (PD) code. |
| `Permutation` | Permutation is a bijection of the set {0, 1, ..., n-1} represented by the slice of images: p`i` is the image of i. |

<details>
<summary><code>Braid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FullTwist(n int) Braid` | FullTwist returns the full twist braid (sigma_1 ... |
| `func GarsideHalfTwist(n int) Braid` | GarsideHalfTwist returns the Garside half-twist braid Delta_n on n strands, the positive half twist whose square is the full twist. |
| `func Generator(n, i int) (Braid, error)` | Generator returns the braid consisting of the single generator sigma_i on n strands. |
| `func GeneratorInv(n, i int) (Braid, error)` | GeneratorInv returns the braid consisting of the single inverse generator sigma_i^{-1} on n strands. |
| `func IdentityBraid(n int) Braid` | IdentityBraid returns the trivial (empty) braid on n strands. |
| `func MustBraid(strands int, word ...int) Braid` | MustBraid is like NewBraid but panics on error. |
| `func NewBraid(strands int, word ...int) (Braid, error)` | NewBraid returns the braid on the given number of strands whose word is the supplied sequence of generators. |
| `func TorusBraid(p, q int) (Braid, error)` | TorusBraid returns the braid (sigma_1 sigma_2 ... |
| `func TorusKnotBraid(p, q int) (Braid, error)` | TorusKnotBraid returns a braid on \|p\| strands whose closure is T(p, q), namely (sigma_1 ... |
| `func (b Braid) AlexanderPolynomial() (Laurent, error)` | AlexanderPolynomial returns the Alexander polynomial of the braid closure, computed from the reduced Burau matrix via the Burau formula |
| `func (b Braid) ClosureCrossingNumber() int` | ClosureCrossingNumber returns the number of crossings in the standard closed braid diagram, which is the braid word length. |
| `func (b Braid) ClosureIsKnot() bool` | ClosureIsKnot reports whether the closure of the braid is a knot (a single component). |
| `func (b Braid) Concat(other Braid) (Braid, error)` | Concat returns the product b*other formed by concatenating the two braid words. |
| `func (b Braid) Conjugate(g Braid) (Braid, error)` | Conjugate returns the braid g * b * g^{-1}. |
| `func (b Braid) ExponentSum() int` | ExponentSum returns the algebraic sum of the generator signs, that is the number of positive generators minus the number of negative ones. |
| `func (b Braid) FreeReduce() Braid` | FreeReduce returns the braid obtained by repeatedly cancelling adjacent inverse pairs sigma_i sigma_i^{-1}. |
| `func (b Braid) Generators() []int` | Generators returns the sorted list of distinct generator indices that occur in the braid word (ignoring sign). |
| `func (b Braid) Inverse() Braid` | Inverse returns the inverse braid, whose word is the reversed word with every generator sign flipped. |
| `func (b Braid) IsFreelyReduced() bool` | IsFreelyReduced reports whether the braid word contains no adjacent inverse pair. |
| `func (b Braid) IsNegative() bool` | IsNegative reports whether every generator in the word is negative. |
| `func (b Braid) IsPositive() bool` | IsPositive reports whether every generator in the word is positive. |
| `func (b Braid) IsPureBraid() bool` | IsPureBraid reports whether the braid induces the identity permutation, so that every strand returns to its starting position. |
| `func (b Braid) IsTrivial() bool` | IsTrivial reports whether the braid word is empty. |
| `func (b Braid) Length() int` | Length returns the number of letters (generators) in the braid word, which equals the number of crossings in the associated braid diagram. |
| `func (b Braid) Mirror() Braid` | Mirror returns the braid obtained by flipping every crossing, that is by negating every generator. |
| `func (b Braid) NumComponents() int` | NumComponents returns the number of link components of the braid closure, which equals the number of cycles of the underlying permutation. |
| `func (b Braid) Permutation() Permutation` | Permutation returns the underlying permutation of the braid: strand starting at position i ends at position Permutation()`i`. |
| `func (b Braid) Power(k int) Braid` | Power returns the braid raised to the integer power k. |
| `func (b Braid) ReducedBurau() *LaurentMatrix` | ReducedBurau returns the reduced Burau matrix of the braid, an (n-1)x(n-1) matrix of Laurent polynomials in t obtained as the ordered product of the… |
| `func (b Braid) Reverse() Braid` | Reverse returns the braid whose word is read right to left (the reverse or "review" of the braid), keeping every generator sign. |
| `func (b Braid) Strands() int` | Strands returns the number of strands of the braid. |
| `func (b Braid) String() string` | String renders the braid word using sigma notation, for example "s1 s2^-1 s1". |
| `func (b Braid) Word() []int` | Word returns a copy of the generator word of the braid. |
| `func (b Braid) Writhe() int` | Writhe returns the writhe of the closure of the braid, which equals the exponent sum. |

</details>

<details>
<summary><code>Diagram</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewDiagram(components ...GaussCode) (Diagram, error)` | NewDiagram validates and returns a Diagram from its component Gauss codes. |
| `func (d Diagram) CrossingNumber() int` | CrossingNumber returns the number of crossings in the diagram. |
| `func (d Diagram) IsKnot() bool` | IsKnot reports whether the diagram has exactly one component. |
| `func (d Diagram) IsLink() bool` | IsLink reports whether the diagram has more than one component. |
| `func (d Diagram) LinkingNumber(i, j int) (int, error)` | LinkingNumber returns the Gauss linking number of components i and j, defined as half the sum of the signs of the crossings between the two… |
| `func (d Diagram) NumComponents() int` | NumComponents returns the number of components of the diagram. |
| `func (d Diagram) SelfWrithe(comp int) int` | SelfWrithe returns the writhe contributed only by self-crossings of the given component (crossings both of whose strands lie on that component). |
| `func (d Diagram) TotalLinkingNumber() int` | TotalLinkingNumber returns the sum of the pairwise linking numbers over all unordered pairs of distinct components. |
| `func (d Diagram) Writhe() int` | Writhe returns the writhe of the whole diagram, the sum of all crossing signs. |

</details>

<details>
<summary><code>GaussCode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewGaussCode(entries []GaussEntry) (GaussCode, error)` | NewGaussCode validates and returns a GaussCode. |
| `func (gc GaussCode) CrossingLabels() []int` | CrossingLabels returns the sorted list of distinct crossing labels. |
| `func (gc GaussCode) CrossingNumber() int` | CrossingNumber returns the number of crossings of the diagram. |
| `func (gc GaussCode) IsReduced() bool` | IsReduced reports whether the diagram admits neither a type I kink nor a reducible type II bigon, a necessary (though not sufficient) condition for… |
| `func (gc GaussCode) IsReidemeisterIIIApplicable() bool` | IsReidemeisterIIIApplicable reports whether the diagram contains a candidate type III triangle. |
| `func (gc GaussCode) IsReidemeisterIIReducible() bool` | IsReidemeisterIIReducible reports whether a reducible type II bigon exists. |
| `func (gc GaussCode) IsReidemeisterIReducible() bool` | IsReidemeisterIReducible reports whether the diagram contains at least one type I kink. |
| `func (gc GaussCode) Length() int` | Length returns the number of crossing encounters, which is twice the number of crossings. |
| `func (gc GaussCode) Mirror() GaussCode` | Mirror returns the Gauss code of the mirror image, obtained by exchanging over and under at every crossing and negating every sign. |
| `func (gc GaussCode) ReidemeisterIIITriangles() [][3]int` | ReidemeisterIIITriangles returns the triples of distinct crossing labels that are pairwise cyclically adjacent in the Gauss code. |
| `func (gc GaussCode) ReidemeisterIIReducible() [][2]int` | ReidemeisterIIReducible returns the pairs of crossing labels that form a reducible Reidemeister II bigon. |
| `func (gc GaussCode) ReidemeisterIKinks() []int` | ReidemeisterIKinks returns the crossing labels that form a Reidemeister I kink, that is a crossing whose two encounters are cyclically adjacent in… |
| `func (gc GaussCode) RemoveReidemeisterII(a, b int) (GaussCode, bool)` | RemoveReidemeisterII returns the Gauss code with the two crossings of a reducible bigon removed. |
| `func (gc GaussCode) RemoveReidemeisterIKink(crossing int) (GaussCode, bool)` | RemoveReidemeisterIKink returns the Gauss code with the kink at the given crossing removed. |
| `func (gc GaussCode) ReverseOrientation() GaussCode` | ReverseOrientation returns the Gauss code obtained by traversing the knot in the opposite direction. |
| `func (gc GaussCode) ToDiagram() Diagram` | ToDiagram wraps a single-component Gauss code as a one-component Diagram. |
| `func (gc GaussCode) Validate() error` | Validate checks the structural consistency of the Gauss code and returns a descriptive error if it is malformed. |
| `func (gc GaussCode) Writhe() int` | Writhe returns the writhe of the diagram, the sum of the crossing signs. |

</details>

<details>
<summary><code>Laurent</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LaurentConst(c int) Laurent` | LaurentConst returns the constant Laurent polynomial c (that is c*X^0). |
| `func LaurentProduct(ps ...Laurent) Laurent` | LaurentProduct returns the product of all the supplied Laurent polynomials. |
| `func LaurentSum(ps ...Laurent) Laurent` | LaurentSum returns the sum of all the supplied Laurent polynomials. |
| `func Monomial(coeff, exp int) Laurent` | Monomial returns the single-term Laurent polynomial coeff*X^exp. |
| `func NewLaurent(min int, coeffs []int) Laurent` | NewLaurent builds a Laurent polynomial whose coefficient of X^(min+i) is coeffs`i`. |
| `func OneLaurent() Laurent` | OneLaurent returns the constant Laurent polynomial 1. |
| `func SkeinDelta() Laurent` | SkeinDelta returns the loop value delta = -A^2 - A^{-2} of the Kauffman bracket skein relation, as a Laurent polynomial in A. |
| `func TorusKnotAlexander(p, q int) (Laurent, error)` | TorusKnotAlexander returns the Alexander polynomial of the torus knot T(p, q) with coprime p, q >= 2, given by the closed form |
| `func TorusKnotJones(p, q int) (Laurent, error)` | TorusKnotJones returns the Jones polynomial of the torus knot T(p, q) with coprime p, q >= 2, in the variable t, using the closed form |
| `func ZeroLaurent() Laurent` | ZeroLaurent returns the zero Laurent polynomial. |
| `func (L Laurent) Add(other Laurent) Laurent` | Add returns L+other. |
| `func (L Laurent) Clone() Laurent` | Clone returns an independent copy of L. |
| `func (L Laurent) Coeff(exp int) int` | Coeff returns the coefficient of X^exp in L. |
| `func (L Laurent) CoeffSlice() (min int, coeffs []int)` | CoeffSlice returns the minimum exponent and a copy of the dense coefficient slice of L, so that coeffs`i` is the coefficient of X^(min+i). |
| `func (L Laurent) ConstantTerm() int` | ConstantTerm returns the coefficient of X^0. |
| `func (L Laurent) Content() int` | Content returns the greatest common divisor of the coefficients of L (the integer content), and 0 for the zero polynomial. |
| `func (L Laurent) Derivative() Laurent` | Derivative returns the formal derivative dL/dX, itself a Laurent polynomial. |
| `func (L Laurent) DivExact(other Laurent) (q Laurent, ok bool)` | DivExact returns the quotient L/other when the division is exact, and reports ok=false otherwise. |
| `func (L Laurent) Equal(other Laurent) bool` | Equal reports whether L and other are the same polynomial. |
| `func (L Laurent) EqualUpToShift(other Laurent) bool` | EqualUpToShift reports whether L and other agree after shifting one of them by a power of X, that is whether they are equal as polynomials up to… |
| `func (L Laurent) Eval(x float64) float64` | Eval evaluates L at the real number x. |
| `func (L Laurent) EvalUnit(x int) int` | EvalUnit evaluates L at an integer root of unity argument, meaning it returns the integer value L(x) for x equal to +1 or -1, where all powers stay… |
| `func (L Laurent) IsConstant() bool` | IsConstant reports whether L is a constant (a single term of degree zero, or the zero polynomial). |
| `func (L Laurent) IsOne() bool` | IsOne reports whether L is the constant polynomial 1. |
| `func (L Laurent) IsPalindromic() bool` | IsPalindromic reports whether L(X) equals L(X^{-1}) up to multiplication by a power of X. The Alexander polynomial of a knot is always palindromic. |
| `func (L Laurent) IsZero() bool` | IsZero reports whether L is the zero polynomial. |
| `func (L Laurent) LeadingCoeff() int` | LeadingCoeff returns the coefficient of the highest-degree term, or 0 for the zero polynomial. |
| `func (L Laurent) MaxAbsCoeff() int` | MaxAbsCoeff returns the largest absolute value among the coefficients of L, and 0 for the zero polynomial. |
| `func (L Laurent) MaxDegree() int` | MaxDegree returns the largest exponent that occurs with a non-zero coefficient. |
| `func (L Laurent) MinDegree() int` | MinDegree returns the smallest exponent that occurs with a non-zero coefficient. |
| `func (L Laurent) Mul(other Laurent) Laurent` | Mul returns the product L*other. |
| `func (L Laurent) MulMonomial(c, e int) Laurent` | MulMonomial returns L multiplied by the single term c*X^e. |
| `func (L Laurent) Neg() Laurent` | Neg returns -L. |
| `func (L Laurent) NumTerms() int` | NumTerms returns the number of non-zero terms of L. |
| `func (L Laurent) Pow(n int) Laurent` | Pow returns L raised to the non-negative power n. |
| `func (L Laurent) Reverse() Laurent` | Reverse returns L(X^{-1}), the polynomial obtained by replacing X with its inverse. |
| `func (L Laurent) Scale(k int) Laurent` | Scale returns k*L. |
| `func (L Laurent) ShiftExp(d int) Laurent` | ShiftExp returns X^d * L, that is L with every exponent shifted by d. |
| `func (L Laurent) SpanWidth() int` | SpanWidth returns MaxDegree-MinDegree, the breadth of the polynomial. |
| `func (L Laurent) String() string` | String renders L in the variable X using conventional notation, for example "X^-2 - X + 3". |
| `func (L Laurent) StringVar(v string) string` | StringVar renders L using the supplied variable name. |
| `func (L Laurent) Sub(other Laurent) Laurent` | Sub returns L-other. |
| `func (L Laurent) SubstitutePow(k int) Laurent` | SubstitutePow returns L with the variable X replaced by X^k, mapping every term c*X^e to c*X^{e*k}. |
| `func (L Laurent) Terms() []LaurentTerm` | Terms returns the non-zero terms of L sorted by increasing exponent. |
| `func (L Laurent) TrailingCoeff() int` | TrailingCoeff returns the coefficient of the lowest-degree term, or 0 for the zero polynomial. |

</details>

<details>
<summary><code>LaurentMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IdentityLaurentMatrix(n int) *LaurentMatrix` | IdentityLaurentMatrix returns the n-by-n identity matrix over the Laurent ring. |
| `func NewLaurentMatrix(r, c int) *LaurentMatrix` | NewLaurentMatrix returns an r-by-c matrix with every entry the zero polynomial. |
| `func (m *LaurentMatrix) At(i, j int) Laurent` | At returns the entry in row i, column j. |
| `func (m *LaurentMatrix) Clone() *LaurentMatrix` | Clone returns an independent copy of the matrix. |
| `func (m *LaurentMatrix) Cols() int` | Cols returns the number of columns of the matrix. |
| `func (m *LaurentMatrix) Determinant() Laurent` | Determinant returns the determinant of a square matrix computed by Laplace (cofactor) expansion over the Laurent ring. |
| `func (m *LaurentMatrix) Mul(other *LaurentMatrix) *LaurentMatrix` | Mul returns the matrix product m*other. |
| `func (m *LaurentMatrix) Rows() int` | Rows returns the number of rows of the matrix. |
| `func (m *LaurentMatrix) Set(i, j int, v Laurent)` | Set stores v in row i, column j. |
| `func (m *LaurentMatrix) Sub(other *LaurentMatrix) *LaurentMatrix` | Sub returns the entrywise difference m-other. |

</details>

<details>
<summary><code>PDCode</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CinquefoilPD() PDCode` | CinquefoilPD returns the (2,5) torus knot 5_1 (Solomon's seal knot) with writhe +5. |
| `func FigureEightPD() PDCode` | FigureEightPD returns the figure-eight knot (4_1) diagram, an amphichiral knot with writhe 0. |
| `func HopfLinkPD(positive bool) PDCode` | HopfLinkPD returns a Hopf link diagram. |
| `func LeftTrefoilPD() PDCode` | LeftTrefoilPD returns the left-handed trefoil, the mirror image of TrefoilPD. |
| `func TrefoilPD() PDCode` | TrefoilPD returns a right-handed trefoil (3_1) diagram with writhe +3. |
| `func UnknotPD() PDCode` | UnknotPD returns the standard crossingless diagram of the unknot. |
| `func UnlinkPD(n int) PDCode` | UnlinkPD returns the diagram of the n-component unlink (n crossingless circles). |
| `func (pd PDCode) AlexanderMatrix() *LaurentMatrix` | AlexanderMatrix returns the Alexander matrix of the diagram over the Laurent ring in the variable t. |
| `func (pd PDCode) AlexanderPolynomial() (Laurent, error)` | AlexanderPolynomial returns the Alexander polynomial of the diagram computed from its Alexander matrix as the determinant of the minor obtained by… |
| `func (pd PDCode) ArcCount() int` | ArcCount returns the number of arcs of the diagram (maximal over-strands between under-crossings). |
| `func (pd PDCode) ArfInvariant() int` | ArfInvariant returns the Arf invariant (0 or 1) of the knot, computed from the determinant: it is 0 when the determinant is congruent to +/-1 modulo… |
| `func (pd PDCode) IsPColorable(p int) bool` | IsPColorable reports whether the diagram admits a non-trivial Fox p-colouring, which for a prime p is equivalent to p dividing the determinant. |
| `func (pd PDCode) IsThreeColorable() bool` | IsThreeColorable reports whether the diagram is Fox 3-colourable. |
| `func (pd PDCode) JonesPolynomial() (Laurent, error)` | JonesPolynomial returns the Jones polynomial in the variable t. |
| `func (pd PDCode) JonesPolynomialSqrt() Laurent` | JonesPolynomialSqrt returns the Jones polynomial expressed in the variable u = t^(1/2), so that a term c*u^k denotes c*t^(k/2). |
| `func (pd PDCode) JonesSpan() int` | JonesSpan returns the breadth of the Jones polynomial (in the variable t^{1/2}), a lower bound of twice the crossing number for a reduced alternating… |
| `func (pd PDCode) KauffmanBracket() Laurent` | KauffmanBracket returns the Kauffman bracket polynomial of the diagram as a Laurent polynomial in the variable A. It is computed by the state sum |
| `func (pd PDCode) KnotDeterminant() int` | KnotDeterminant returns the determinant of the diagram, the absolute value of the Alexander polynomial evaluated at t = -1. |
| `func (pd PDCode) Mirror() PDCode` | Mirror returns the PD code of the mirror image, obtained by swapping the over- and under-arcs (A,B,C,D -> B,C,D,A) and negating every sign. |
| `func (pd PDCode) NormalizedBracket() Laurent` | NormalizedBracket returns the writhe-normalised Kauffman bracket (-A^3)^(-w) <D>, an ambient-isotopy invariant (in the variable A) from which the… |
| `func (pd PDCode) NumCrossings() int` | NumCrossings returns the number of crossings of the diagram. |
| `func (pd PDCode) SignSequence() []int` | SignSequence returns the signs of the crossings of the PD diagram in order. |
| `func (pd PDCode) Validate() error` | Validate checks that every edge label occurs exactly twice and that every sign is +1 or -1. |
| `func (pd PDCode) Writhe() int` | Writhe returns the writhe of the diagram, the sum of the crossing signs. |

</details>

<details>
<summary><code>Permutation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IdentityPermutation(n int) Permutation` | IdentityPermutation returns the identity permutation on n points. |
| `func NewPermutation(images []int) (Permutation, error)` | NewPermutation returns the permutation with the given image slice after checking that it is a genuine bijection of {0, ..., n-1}. |
| `func TranspositionPermutation(n, i, j int) (Permutation, error)` | TranspositionPermutation returns the permutation on n points that swaps i and j and fixes everything else. |
| `func (p Permutation) Apply(i int) int` | Apply returns the image of i under the permutation. |
| `func (p Permutation) Compose(q Permutation) Permutation` | Compose returns the permutation p∘q defined by (p∘q)(i) = p(q(i)). |
| `func (p Permutation) Cycles() [][]int` | Cycles returns the cycle decomposition of p as a slice of cycles, each cycle being the list of points visited starting from its smallest element. |
| `func (p Permutation) Equal(q Permutation) bool` | Equal reports whether p and q are the same permutation. |
| `func (p Permutation) Inverse() Permutation` | Inverse returns the inverse permutation. |
| `func (p Permutation) IsIdentity() bool` | IsIdentity reports whether p fixes every point. |
| `func (p Permutation) NumCycles() int` | NumCycles returns the number of cycles in the decomposition of p, counting fixed points. |
| `func (p Permutation) Order() int` | Order returns the multiplicative order of p, that is the least positive k with p^k the identity. |
| `func (p Permutation) Power(k int) Permutation` | Power returns p raised to the integer power k, which may be negative. |
| `func (p Permutation) Sign() int` | Sign returns +1 if p is an even permutation and -1 if it is odd. |
| `func (p Permutation) Size() int` | Size returns the number of points the permutation acts on. |
| `func (p Permutation) String() string` | String renders p in one-line notation, for example "[1 2 0]". |
| `func (p Permutation) TranspositionCount() int` | TranspositionCount returns the minimum number of transpositions needed to express the permutation, equal to size minus the number of cycles. |
| `func (p Permutation) Transpositions() [][2]int` | Transpositions returns a decomposition of the permutation into transpositions. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/knottheory).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
