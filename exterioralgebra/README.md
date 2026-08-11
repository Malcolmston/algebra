# exterioralgebra — implements the exterior (Grassmann) algebra Λ(Rⁿ) together with the calculus

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/exterioralgebra.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/exterioralgebra)

Package exterioralgebra implements the exterior (Grassmann) algebra Λ(Rⁿ)
together with the calculus of differential forms on Rⁿ.

The central value type is `Form`, an element of the exterior algebra over an
n-dimensional real vector space. A Form is a (possibly mixed-grade)
multivector: a real linear combination of basis blades
e_{i₁}∧…∧e_{i_k}, with 0 ≤ i₁ < … < i_k < n. Each basis blade is
encoded compactly as a bitmask whose set bits are the participating indices,
so the grade of a blade is the population count of its mask. Coefficients are
stored sparsely, so only the nonzero terms occupy memory.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/exterioralgebra"
```

## Usage

This is the package's own `ExampleBasis1`, so it compiles and its output is
asserted on every `go test ./exterioralgebra/`.

```go
e0 := ea.Basis1(3, 0)
	e1 := ea.Basis1(3, 1)
	fmt.Println(e0.Wedge(e1))
	fmt.Println(e1.Wedge(e0))
	fmt.Println(e0.Wedge(e0))
```

```
1 e0∧e1
-1 e0∧e1
0
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AlgebraDimension(n int) int` | AlgebraDimension returns the dimension 2ⁿ of the full exterior algebra Λ(Rⁿ) as a real vector space. |
| `func Angle(f, g *Form) float64` | Angle returns the angle in radians between two nonzero homogeneous Forms of the same grade, defined via the Euclidean inner product. |
| `func Binomial(n, k int) int` | Binomial returns the binomial coefficient C(n,k), the number of k-element subsets of an n-element set. |
| `func CrossProduct(a, b []float64) ([]float64, error)` | CrossProduct returns the 3-dimensional cross product a×b, computed as the Hodge dual of the wedge a∧b in Λ(R³). |
| `func Determinant(rows ...[]float64) (float64, error)` | Determinant returns the determinant of the n×n matrix whose rows are the n supplied vectors of Rⁿ, computed as the top-grade coefficient of their… |
| `func Dot(f, g *Form) float64` | Dot is an alias for `InnerProduct`. |
| `func EuclideanSignature(n int) []int` | EuclideanSignature returns the all-+1 metric signature of length n. |
| `func FullMask(n int) uint` | FullMask returns the bitmask 0b1…1 with the low n bits set, encoding the top-grade basis blade e_0∧…∧e_{n-1} of Λ(Rⁿ). |
| `func GradeDimension(n, k int) int` | GradeDimension returns the dimension of the grade-k subspace Λ^k(Rⁿ), which equals C(n,k). |
| `func HodgeSignSquared(k, n int) int` | HodgeSignSquared returns the scalar ε with ★★ = ε·id on grade-k Forms of Λ(Rⁿ) under the Euclidean metric, namely (−1)^{k(n−k)}. |
| `func IndicesToMask(n int, idx ...int) (mask uint, sign int, ok bool)` | IndicesToMask converts an ordered list of basis indices into the bitmask of the corresponding blade together with the sign of the permutation that… |
| `func InnerProduct(f, g *Form) float64` | InnerProduct returns the Euclidean inner product ⟨f, g⟩ obtained by declaring the basis blades orthonormal: it is the sum over blades of the… |
| `func LorentzSignature(n int) []int` | LorentzSignature returns the Lorentzian signature (−1,+1,…,+1) of length n, with the time index 0 carrying the −1. |
| `func MaskToIndices(mask uint) []int` | MaskToIndices returns the sorted slice of basis indices whose bits are set in mask. |
| `func NumDirectionalDerivative(f ScalarField, x, v []float64, h float64) float64` | NumDirectionalDerivative returns the directional derivative of the scalar field f at x along the direction vector v, computed as a symmetric finite… |
| `func Popcount(mask uint) int` | Popcount returns the number of set bits in mask, i.e. |
| `func WedgeSign(p, q int) int` | WedgeSign returns the graded-commutation sign (−1)^{pq} relating α∧β to β∧α for a grade-p and a grade-q homogeneous Form. |

### Types

| Type | What it is |
| --- | --- |
| `Form` | Form is an element of the exterior algebra Λ(Rⁿ): a real linear combination of basis blades. |
| `FormField` | FormField is a differential form on Rⁿ given pointwise: it maps a position x to the constant `Form` of coefficients at that point. |
| `PForm` | PForm is a differential form on Rⁿ whose coefficients are polynomials: a finite sum Σ_I P_I(x) dx^I where each P_I is a `Poly` in the n… |
| `Poly` | Poly is a multivariate polynomial in n real variables x_0,…,x_{n-1}. |
| `ScalarField` | ScalarField is a smooth real-valued function on Rⁿ, used as the coefficient of a numerically differentiated differential form. |
| `Term` | Term is one basis blade of a `Form` together with its coefficient. |

<details>
<summary><code>Form</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AddAll(forms ...*Form) *Form` | AddAll returns the sum of any number of Forms, all of which must share an ambient dimension. |
| `func AllBasisBlades(n int) []*Form` | AllBasisBlades returns all 2ⁿ basis blades of Λ(Rⁿ) as unit Forms, ordered by grade and then by mask. |
| `func Basis(n, k int) []*Form` | Basis returns every grade-k basis blade of Λ(Rⁿ) as a unit-coefficient Form, listed in canonical (increasing mask) order. |
| `func Basis1(n, i int) *Form` | Basis1 returns the grade-1 basis Form e_i in Λ(Rⁿ). |
| `func BasisBlade(n int, idx ...int) (*Form, error)` | BasisBlade returns the basis blade e_{idx[0]}∧…∧e_{idx`k-1`} in Λ(Rⁿ) with unit coefficient, reduced to its canonical sorted representative… |
| `func BasisBladeCoeff(n int, c float64, idx ...int) (*Form, error)` | BasisBladeCoeff returns c·e_{idx[0]}∧…∧e_{idx`k-1`} in Λ(Rⁿ). |
| `func FromMasks(n int, terms map[uint]float64) (*Form, error)` | FromMasks builds a Form in Λ(Rⁿ) directly from a map of blade bitmasks to coefficients. |
| `func FromVector(v []float64) *Form` | FromVector returns the grade-1 Form whose e_i coefficient is v`i`. |
| `func HodgeStar(f *Form) *Form` | HodgeStar returns the Euclidean Hodge dual ★f, computed with the positively oriented orthonormal frame and volume element e_0∧…∧e_{n-1}. |
| `func HodgeStarMetric(f *Form, signature []int) (*Form, error)` | HodgeStarMetric returns the Hodge dual of f with respect to a pseudo-Euclidean metric whose diagonal signature entries are each +1 or −1. |
| `func InteriorProduct(v, omega *Form) *Form` | InteriorProduct returns the interior product (contraction) ι_v ω of the vector v into the Form omega, lowering grade by one. |
| `func InverseHodgeStar(f *Form) *Form` | InverseHodgeStar returns the inverse Euclidean Hodge dual, the unique operator with InverseHodgeStar(HodgeStar(f)) == f. |
| `func LinearCombination(coeffs []float64, forms []*Form) *Form` | LinearCombination returns Σ coeffs`i`·forms`i`. |
| `func New(n int) *Form` | New returns the zero Form of the exterior algebra Λ(Rⁿ). |
| `func NumExteriorDerivative(field FormField, x []float64, h float64) *Form` | NumExteriorDerivative returns the exterior derivative of the form field at the point x, evaluated with symmetric finite differences of step h. |
| `func NumGradient(f ScalarField, x []float64, h float64) *Form` | NumGradient returns the exterior derivative of the 0-form f at the point x, namely the grade-1 Form Σ_i (∂f/∂x_i) dx^i, using a symmetric finite… |
| `func One(n int) *Form` | One returns the multiplicative identity of Λ(Rⁿ), the scalar 1. |
| `func Scalar(n int, c float64) *Form` | Scalar returns the grade-0 Form with the given scalar coefficient in Λ(Rⁿ). |
| `func VolumeForm(n int) *Form` | VolumeForm returns the canonical top-grade unit volume element e_0∧…∧e_{n-1} of Λ(Rⁿ). |
| `func Wedge(a, b *Form) *Form` | Wedge returns the exterior product a∧b as a package-level convenience. |
| `func WedgeAll(forms ...*Form) *Form` | WedgeAll returns the exterior product forms[0]∧forms[1]∧… evaluated left to right. |
| `func WedgeVectors(vs ...[]float64) (*Form, error)` | WedgeVectors returns the exterior product of the given vectors, each a grade-1 Form built from a coefficient slice. |
| `func Zero(n int) *Form` | Zero is an alias for `New`: it returns the additive identity of Λ(Rⁿ). |
| `func (f *Form) Add(g *Form) *Form` | Add returns the sum f+g. |
| `func (f *Form) BivectorPart() *Form` | BivectorPart returns the grade-2 part of f. |
| `func (f *Form) Clone() *Form` | Clone returns an independent deep copy of f. |
| `func (f *Form) Coeff(idx ...int) float64` | Coeff returns the coefficient of the blade e_{idx[0]}∧…∧e_{idx`k-1`}, accounting for the sign introduced by sorting the indices. |
| `func (f *Form) CoeffMask(mask uint) float64` | CoeffMask returns the coefficient of the sorted basis blade encoded by mask. |
| `func (f *Form) Coeffs() map[uint]float64` | Coeffs returns a fresh copy of the blade→coefficient map underlying f. |
| `func (f *Form) Dim() int` | Dim returns the ambient dimension n of the exterior algebra Λ(Rⁿ) that f belongs to. |
| `func (f *Form) Equal(g *Form) bool` | Equal reports whether f and g are exactly equal: same ambient dimension and identical coefficients on every blade. |
| `func (f *Form) EqualTol(g *Form, tol float64) bool` | EqualTol reports whether f and g have the same ambient dimension and every blade coefficient agrees to within an absolute tolerance tol. |
| `func (f *Form) EvenPart() *Form` | EvenPart returns the sum of the even-grade blades of f (the even subalgebra component). |
| `func (f *Form) Grade() (k int, ok bool)` | Grade returns the common grade of a homogeneous Form. |
| `func (f *Form) GradeInvolution() *Form` | GradeInvolution returns the grade involution α(f), which negates every odd-grade blade: on a grade-k blade it multiplies by (−1)^k. |
| `func (f *Form) GradeProject(k int) *Form` | GradeProject returns the homogeneous grade-k part of f, a new Form containing only the blades of exactly k indices. |
| `func (f *Form) Grades() []int` | Grades returns the sorted list of distinct grades present in f. |
| `func (f *Form) Interior(v *Form) *Form` | Interior is the method form of `InteriorProduct`: omega.Interior(v) returns ι_v omega. |
| `func (f *Form) IsBlade() bool` | IsBlade reports whether f is a single basis blade times a scalar (at most one nonzero term), and hence trivially decomposable. |
| `func (f *Form) IsHomogeneous() bool` | IsHomogeneous reports whether every nonzero term of f has the same grade. |
| `func (f *Form) IsScalar() bool` | IsScalar reports whether f has only a grade-0 part (it may be the zero Form). |
| `func (f *Form) IsVector() bool` | IsVector reports whether every nonzero blade of f has grade 1. |
| `func (f *Form) IsZero() bool` | IsZero reports whether f is the zero Form (has no nonzero terms). |
| `func (f *Form) LeftContract(g *Form) *Form` | LeftContract returns the left contraction f ⌋ g under the Euclidean metric, the bilinear product characterised by ⟨f⌋g, h⟩ = ⟨g, f∧h⟩. |
| `func (f *Form) Masks() []uint` | Masks returns the bitmasks of the nonzero blades of f, ordered by grade and then by mask. |
| `func (f *Form) MaxGrade() int` | MaxGrade returns the largest grade present in f, or 0 for the zero Form. |
| `func (f *Form) MinGrade() int` | MinGrade returns the smallest grade present in f, or 0 for the zero Form. |
| `func (f *Form) Neg() *Form` | Neg returns the additive inverse -f. |
| `func (f *Form) Norm() float64` | Norm returns the Euclidean norm √⟨f, f⟩ of f. |
| `func (f *Form) NormSq() float64` | NormSq returns the squared Euclidean norm ⟨f, f⟩ of f. |
| `func (f *Form) Normalize() (*Form, float64)` | Normalize returns f scaled to unit Euclidean norm, together with the original norm. |
| `func (f *Form) NumTerms() int` | NumTerms returns the number of nonzero basis blades in f. |
| `func (f *Form) OddPart() *Form` | OddPart returns the sum of the odd-grade blades of f. |
| `func (f *Form) Reverse() *Form` | Reverse returns the reversion f̃ of f, which reverses the order of the vector factors in every blade. |
| `func (f *Form) RightContract(g *Form) *Form` | RightContract returns the right contraction f ⌊ g under the Euclidean metric, characterised by ⟨f⌊g, h⟩ = ⟨f, h∧g⟩. |
| `func (f *Form) ScalarPart() float64` | ScalarPart returns the grade-0 coefficient of f as a plain float64. |
| `func (f *Form) ScalarValue() float64` | ScalarValue is an alias for `Form.ScalarPart`: the grade-0 coefficient of f. |
| `func (f *Form) Scale(c float64) *Form` | Scale returns c·f, the scalar multiple of f by c. |
| `func (f *Form) Set(c float64, idx ...int) error` | Set sets the coefficient of the blade e_{idx[0]}∧…∧e_{idx`k-1`} to c, taking the sorting sign into account. |
| `func (f *Form) SetMask(mask uint, c float64)` | SetMask sets the coefficient of the sorted blade encoded by mask to c, removing the term when c == 0. |
| `func (f *Form) Star() *Form` | Star is the method form of `HodgeStar`: f.Star() returns ★f. |
| `func (f *Form) String() string` | String renders f in the usual e_{i…} basis notation, for example "2 + 3 e0∧e1". |
| `func (f *Form) Sub(g *Form) *Form` | Sub returns the difference f-g. |
| `func (f *Form) Terms() []Term` | Terms returns the nonzero blades of f as a slice of `Term`, ordered by grade and then by mask. |
| `func (f *Form) ToVector() []float64` | ToVector returns the grade-1 part of f as a coefficient slice of length n, so that the i-th entry is the coefficient of e_i. |
| `func (f *Form) TopPart() *Form` | TopPart returns the top-grade (grade-n) part of f as a new Form. |
| `func (f *Form) VectorPart() *Form` | VectorPart returns the grade-1 part of f as a new Form. |
| `func (f *Form) Wedge(g *Form) *Form` | Wedge returns the exterior product f∧g. |
| `func (f *Form) WedgePow(k int) *Form` | WedgePow returns the k-fold exterior power f∧…∧f (k factors). |

</details>

<details>
<summary><code>PForm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Dx(n, i int) *PForm` | Dx returns the constant basis 1-form dx^i on Rⁿ. |
| `func NewPForm(n int) *PForm` | NewPForm returns the zero differential form on Rⁿ. |
| `func PBasisBlade(n int, coeff *Poly, idx ...int) (*PForm, error)` | PBasisBlade returns the differential form coeff·dx^{idx[0]}∧… on Rⁿ, reducing the blade to canonical sorted order with the appropriate sign. |
| `func PConst(p *Poly) *PForm` | PConst returns the grade-0 differential form whose single coefficient is the polynomial p. |
| `func (w *PForm) Add(u *PForm) *PForm` | Add returns the sum w+u. |
| `func (w *PForm) Clone() *PForm` | Clone returns an independent deep copy of w. |
| `func (w *PForm) Codifferential() *PForm` | Codifferential returns the codifferential δw with respect to the Euclidean metric, defined grade-by-grade on a grade-k component by δ =… |
| `func (w *PForm) CoeffMask(mask uint) *Poly` | CoeffMask returns the polynomial coefficient of the sorted blade encoded by mask, or the zero polynomial if that blade is absent. |
| `func (w *PForm) D() *PForm` | D is a short alias for `PForm.ExteriorDerivative`. |
| `func (w *PForm) Dim() int` | Dim returns the ambient dimension n of the form. |
| `func (w *PForm) Equal(u *PForm) bool` | Equal reports whether w and u are the same polynomial form. |
| `func (w *PForm) EqualTol(u *PForm, tol float64) bool` | EqualTol reports whether w and u have the same ambient dimension and every coefficient polynomial agrees to within the absolute tolerance tol. |
| `func (w *PForm) Eval(x []float64) *Form` | Eval evaluates every coefficient polynomial of w at the point x and returns the resulting constant `Form`. |
| `func (w *PForm) ExteriorDerivative() *PForm` | ExteriorDerivative returns the exterior derivative dw. |
| `func (w *PForm) GradeProject(k int) *PForm` | GradeProject returns the homogeneous grade-k part of w. |
| `func (w *PForm) Grades() []int` | Grades returns the sorted list of distinct grades present in the polynomial form w. |
| `func (w *PForm) HodgeLaplacian() *PForm` | HodgeLaplacian returns the Hodge–de Rham Laplacian Δw = (dδ + δd)w with respect to the Euclidean metric. |
| `func (w *PForm) HodgeStar() *PForm` | HodgeStar returns the Euclidean Hodge dual ★w of the polynomial form, applying the blade-level Hodge map coefficient-wise. |
| `func (w *PForm) InteriorConst(v []float64) *PForm` | InteriorConst returns the interior product ι_v w of the polynomial form w with the constant vector v of length n, lowering grade by one. |
| `func (w *PForm) IsClosed() bool` | IsClosed reports whether w is a closed form, that is dw = 0. |
| `func (w *PForm) IsHomogeneous() bool` | IsHomogeneous reports whether every nonzero blade of w has the same grade. |
| `func (w *PForm) IsZero() bool` | IsZero reports whether w has no nonzero coefficient polynomials. |
| `func (w *PForm) Masks() []uint` | Masks returns the bitmasks of the nonzero blades of w, ordered by grade and then by mask. |
| `func (w *PForm) MaxGrade() int` | MaxGrade returns the largest grade present in w, or 0 for the zero form. |
| `func (w *PForm) Neg() *PForm` | Neg returns the additive inverse −w. |
| `func (w *PForm) NumTerms() int` | NumTerms returns the number of nonzero blades in w. |
| `func (w *PForm) Pullback(phi []*Poly) (*PForm, error)` | Pullback returns the pullback φ*w of w along the polynomial map φ whose n component polynomials are given in phi, each of the same arity m. |
| `func (w *PForm) Scale(c float64) *PForm` | Scale returns c·w, the real scalar multiple of w. |
| `func (w *PForm) ScalePoly(p *Poly) *PForm` | ScalePoly returns the product of w with the scalar polynomial p, multiplying every coefficient by p. |
| `func (w *PForm) String() string` | String renders w with each blade shown as its coefficient polynomial times the differential dx^I, for example "(2 x0) dx0 + (1) dx0∧dx1". |
| `func (w *PForm) Sub(u *PForm) *PForm` | Sub returns the difference w−u. |
| `func (w *PForm) Wedge(u *PForm) *PForm` | Wedge returns the exterior product w∧u of two polynomial forms. |

</details>

<details>
<summary><code>Poly</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ConstPoly(n int, c float64) *Poly` | ConstPoly returns the constant polynomial c in n variables. |
| `func Monomial(coeff float64, exp []int) *Poly` | Monomial returns the polynomial coeff·∏ x_i^{exp`i`} in n = len(exp) variables. |
| `func NewPoly(n int) *Poly` | NewPoly returns the zero polynomial in n variables. |
| `func PolyCurl(field []*Poly) ([]*Poly, error)` | PolyCurl returns the curl of a 3-dimensional vector field field = (P,Q,R), namely (∂R/∂x₁−∂Q/∂x₂, ∂P/∂x₂−∂R/∂x₀,… |
| `func PolyDivergence(field []*Poly) (*Poly, error)` | PolyDivergence returns the divergence Σ_i ∂F_i/∂x_i of the vector field whose components are field. |
| `func PolyGradient(f *Poly) []*Poly` | PolyGradient returns the gradient of the scalar field f as its n component polynomials, with entry i equal to ∂f/∂x_i. |
| `func PolyLaplacian(f *Poly) *Poly` | PolyLaplacian is a package-level alias for `Poly.Laplacian`. |
| `func Var(n, i int) *Poly` | Var returns the polynomial equal to the single variable x_i in n variables. |
| `func Vars(n int) []*Poly` | Vars returns the n coordinate polynomials x_0,…,x_{n-1} of arity n. |
| `func (p *Poly) Add(q *Poly) *Poly` | Add returns the sum p+q. |
| `func (p *Poly) Arity() int` | Arity returns the number of variables n of p. |
| `func (p *Poly) Clone() *Poly` | Clone returns an independent deep copy of p. |
| `func (p *Poly) Coefficient(exp []int) float64` | Coefficient returns the coefficient of the monomial ∏ x_i^{exp`i`} in p. |
| `func (p *Poly) Compose(subs []*Poly) (*Poly, error)` | Compose substitutes each variable x_i with the polynomial subs`i`, returning the composed polynomial in the common arity of the substitutes. |
| `func (p *Poly) ConstantTerm() float64` | ConstantTerm returns the coefficient of the constant monomial of p. |
| `func (p *Poly) Degree() int` | Degree returns the total degree of p, the largest sum of exponents over its monomials. |
| `func (p *Poly) Equal(q *Poly) bool` | Equal reports whether p and q are the same polynomial (same arity and coefficients). |
| `func (p *Poly) EqualTol(q *Poly, tol float64) bool` | EqualTol reports whether p and q have the same arity and all coefficients agree to within an absolute tolerance tol. |
| `func (p *Poly) Eval(x []float64) float64` | Eval evaluates p at the point x, which must have length n. |
| `func (p *Poly) Gradient() []*Poly` | Gradient returns the n partial derivatives of p as a slice, with entry i equal to ∂p/∂x_i. |
| `func (p *Poly) Hessian() [][]*Poly` | Hessian returns the n×n matrix of second partial derivatives of p, with entry `i``j` equal to ∂²p/∂x_i∂x_j. |
| `func (p *Poly) IsConstant() bool` | IsConstant reports whether p has no variable dependence (degree ≤ 0). |
| `func (p *Poly) IsZero() bool` | IsZero reports whether p is the zero polynomial. |
| `func (p *Poly) Laplacian() *Poly` | Laplacian returns the scalar Laplacian Σ_i ∂²p/∂x_i² of the polynomial p. |
| `func (p *Poly) Mul(q *Poly) *Poly` | Mul returns the product p·q. |
| `func (p *Poly) Neg() *Poly` | Neg returns −p. |
| `func (p *Poly) NumTerms() int` | NumTerms returns the number of nonzero monomials in p. |
| `func (p *Poly) Partial(i int) *Poly` | Partial returns the partial derivative ∂p/∂x_i. |
| `func (p *Poly) Pow(k int) *Poly` | Pow returns p raised to the non-negative integer power k. |
| `func (p *Poly) Scale(c float64) *Poly` | Scale returns c·p. |
| `func (p *Poly) String() string` | String renders p in conventional notation such as "2 x0^2 + 3 x1 - 1". |
| `func (p *Poly) Sub(q *Poly) *Poly` | Sub returns the difference p−q. |

</details>

### Variables

`ErrDim`, `ErrGrade`, `ErrIndex`, `ErrMap`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/exterioralgebra).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
