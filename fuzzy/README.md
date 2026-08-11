# fuzzy — implements fuzzy logic and fuzzy set theory using only the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/fuzzy.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/fuzzy)

Package fuzzy implements fuzzy logic and fuzzy set theory using only the Go
standard library.

The package provides membership functions (triangular, trapezoidal, Gaussian,
sigmoid, generalized bell, S/Z/Pi shapes and singletons), the standard family
of t-norms and t-conorms (minimum/maximum, algebraic product/probabilistic
sum, Lukasiewicz, drastic, Hamacher, Einstein and nilpotent variants), fuzzy
complements (standard, Sugeno, Yager), discrete fuzzy sets with set-theoretic
operations, alpha-cuts, linguistic hedges, scalar cardinality and convexity
tests, fuzzy relations with max-min and max-product (and general max-t)
composition, Mamdani and Sugeno (TSK) inference systems, and a full complement
of defuzzification methods (centroid, bisector, mean/smallest/largest of
maxima and weighted average).

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/fuzzy"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BellAt(x, a, b, c float64) float64` | BellAt evaluates the generalized bell membership function 1/(1+\|(x-c)/a\|^(2b)) at the point x. |
| `func CentroidMF(mf MF, a, b float64, n int) (float64, error)` | CentroidMF returns the centroid of the membership function mf sampled at n points over the closed interval [a, b] using the trapezoidal rule. |
| `func Clamp01(v float64) float64` | Clamp01 returns v restricted to the closed interval [0, 1]. |
| `func ComplementStandard(a float64) float64` | ComplementStandard is the standard (Zadeh) fuzzy complement 1 - a. |
| `func ComplementSugeno(lambda float64) func(float64) float64` | ComplementSugeno returns the Sugeno class complement (1-a)/(1+lambda*a) with parameter lambda > -1. |
| `func ComplementYager(w float64) func(float64) float64` | ComplementYager returns the Yager class complement (1 - a^w)^(1/w) with parameter w > 0. |
| `func DiffSigmoidAt(x, a1, c1, a2, c2 float64) float64` | DiffSigmoidAt evaluates the difference of two sigmoids, clamped to [0, 1], giving an asymmetric bump: sigmoid(a1, c1) minus sigmoid(a2, c2). |
| `func Gaussian2At(x, mean1, sigma1, mean2, sigma2 float64) float64` | Gaussian2At evaluates the two-sided Gaussian membership function that uses (mean1, sigma1) for the left flank and (mean2, sigma2) for the right… |
| `func GaussianAt(x, mean, sigma float64) float64` | GaussianAt evaluates the Gaussian membership function with center mean and width sigma at the point x, exp(-(x-mean)^2 / (2 sigma^2)). |
| `func Linspace(a, b float64, n int) []float64` | Linspace returns n evenly spaced points spanning the closed interval [a, b]. |
| `func PiShapeAt(x, a, b, c, d float64) float64` | PiShapeAt evaluates the Pi-shaped membership function that rises from a to b, stays at 1 on [b, c] and falls from c to d, at the point x. |
| `func ProdSigmoidAt(x, a1, c1, a2, c2 float64) float64` | ProdSigmoidAt evaluates the product of two sigmoids sigmoid(a1, c1) times sigmoid(a2, c2) at the point x, producing a smooth bump. |
| `func RampDownAt(x, a, b float64) float64` | RampDownAt evaluates the linear falling ramp that is 1 at or below a, 0 at or above b and interpolates linearly between, at the point x. |
| `func RampUpAt(x, a, b float64) float64` | RampUpAt evaluates the linear rising ramp that is 0 at or below a, 1 at or above b and interpolates linearly between, at the point x. |
| `func RectangularAt(x, a, b float64) float64` | RectangularAt evaluates the crisp rectangular (interval) membership function that is 1 on the closed interval [a, b] and 0 elsewhere, at the point x. |
| `func SShapeAt(x, a, b float64) float64` | SShapeAt evaluates the spline-based S-shaped membership function with feet a and b at the point x. |
| `func SampleMF(mf MF, xs []float64) []float64` | SampleMF evaluates mf at each point of xs and returns the grades in a new slice. |
| `func SigmoidAt(x, a, c float64) float64` | SigmoidAt evaluates the sigmoidal membership function 1/(1+exp(-a(x-c))) at the point x. |
| `func SingletonAt(x, c float64) float64` | SingletonAt returns 1 when x equals c and 0 otherwise, the membership of a fuzzy singleton (used to fuzzify a crisp input). |
| `func TConormDrastic(a, b float64) float64` | TConormDrastic is the drastic sum t-conorm, the largest t-conorm: it equals a when b is 0, b when a is 0 and 1 otherwise. |
| `func TConormEinstein(a, b float64) float64` | TConormEinstein is the Einstein sum t-conorm (a+b) / (1 + a*b). |
| `func TConormLukasiewicz(a, b float64) float64` | TConormLukasiewicz is the Lukasiewicz (bounded sum) t-conorm min(1, a+b). |
| `func TConormMax(a, b float64) float64` | TConormMax is the maximum (Godel) t-conorm max(a, b), the smallest t-conorm. |
| `func TConormN(s TConorm, vs ...float64) float64` | TConormN aggregates the values vs with the t-conorm s, folding left and starting from the identity 0. |
| `func TConormNilpotentMax(a, b float64) float64` | TConormNilpotentMax is the nilpotent maximum t-conorm: max(a, b) when a+b < 1 and 1 otherwise. |
| `func TConormProbabilistic(a, b float64) float64` | TConormProbabilistic is the probabilistic sum t-conorm a + b - a*b, dual to the algebraic product. |
| `func TNormDrastic(a, b float64) float64` | TNormDrastic is the drastic product t-norm, the smallest t-norm: it equals a when b is 1, b when a is 1 and 0 otherwise. |
| `func TNormEinstein(a, b float64) float64` | TNormEinstein is the Einstein product t-norm a*b / (2 - (a+b-a*b)). |
| `func TNormLukasiewicz(a, b float64) float64` | TNormLukasiewicz is the Lukasiewicz t-norm max(0, a+b-1). |
| `func TNormMin(a, b float64) float64` | TNormMin is the minimum (Godel) t-norm min(a, b), the largest t-norm. |
| `func TNormN(t TNorm, vs ...float64) float64` | TNormN aggregates the values vs with the t-norm t, folding left and starting from the identity 1. |
| `func TNormNilpotentMin(a, b float64) float64` | TNormNilpotentMin is the nilpotent minimum t-norm: min(a, b) when a+b > 1 and 0 otherwise. |
| `func TNormProduct(a, b float64) float64` | TNormProduct is the algebraic product t-norm a*b. |
| `func TrapezoidalAt(x, a, b, c, d float64) float64` | TrapezoidalAt evaluates the trapezoidal membership function with feet a and d and shoulders b and c at the point x. |
| `func TriangularAt(x, a, b, c float64) float64` | TriangularAt evaluates the triangular membership function with feet a and c and peak b at the point x. |
| `func ZShapeAt(x, a, b float64) float64` | ZShapeAt evaluates the spline-based Z-shaped membership function with shoulders a and b at the point x. |

### Types

| Type | What it is |
| --- | --- |
| `Clause` | Clause is one atomic antecedent proposition "Var is `not` Term". |
| `Defuzz` | Defuzz names a defuzzification strategy. |
| `Implication` | Implication selects how a rule's firing strength shapes its consequent set. |
| `MF` | MF is a membership function mapping a point of the universe of discourse to a membership grade. |
| `MamdaniRule` | MamdaniRule is one rule of a Mamdani system: an antecedent (clauses combined by Op), a consequent term OutTerm of the output variable and a rule… |
| `MamdaniSystem` | MamdaniSystem is a Mamdani (max-min) fuzzy inference system with several input variables and a single output variable. |
| `Relation` | Relation is a binary fuzzy relation on the product of universes X and Y. M`i` `j` is the grade to which X`i` is related to Y`j`. |
| `RuleOp` | RuleOp selects how the clauses of a rule antecedent are combined. |
| `Set` | Set is a discrete fuzzy set: a sorted universe of discourse X together with the membership grade Mu`i` of X`i`. |
| `SugenoConsequent` | SugenoConsequent is a first-order Takagi-Sugeno-Kang consequent z = Const + sum(Coeffs`var` * input`var`). |
| `SugenoRule` | SugenoRule is one rule of a Sugeno system: an antecedent (clauses combined by Op), a crisp/linear consequent and a rule weight. |
| `SugenoSystem` | SugenoSystem is a Takagi-Sugeno-Kang fuzzy inference system whose crisp output is the firing-strength weighted average of the rule consequents. |
| `TConorm` | TConorm is a triangular conorm (s-norm), a binary operator on [0, 1] used as the fuzzy logical disjunction (OR) and as fuzzy set union. |
| `TNorm` | TNorm is a triangular norm, a binary operator on [0, 1] used as the fuzzy logical conjunction (AND) and as fuzzy set intersection. |
| `Variable` | Variable is a linguistic variable: a named universe [Min, Max] together with a collection of named linguistic terms, each described by a membership… |

<details>
<summary><code>Clause</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func If(v, t string) Clause` | If is a convenience constructor for a positive antecedent clause "v is t". |
| `func IfNot(v, t string) Clause` | IfNot is a convenience constructor for a negated antecedent clause "v is not t". |

</details>

<details>
<summary><code>Defuzz</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (d Defuzz) String() string` | String returns the canonical name of the defuzzification strategy. |

</details>

<details>
<summary><code>MF</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Bell(a, b, c float64) MF` | Bell returns the generalized bell membership function with width a, slope b and center c. |
| `func ClipMF(mf MF, alpha float64) MF` | ClipMF returns the membership function mf capped at the level alpha, that is min(mf(x), alpha). |
| `func ComplementMF(mf MF) MF` | ComplementMF returns the standard complement 1 - mf(x) of a membership function. |
| `func ComplementMFWith(mf MF, c func(float64) float64) MF` | ComplementMFWith returns the complement of mf using an arbitrary fuzzy complement operator c (for example ComplementSugeno or ComplementYager). |
| `func ConcentrateMF(mf MF) MF` | ConcentrateMF returns the concentration hedge mf(x)^2, sharpening a fuzzy set (the linguistic hedge "very"). |
| `func Constant(v float64) MF` | Constant returns a membership function equal to the constant grade v (clamped to [0, 1]) everywhere. |
| `func DiffSigmoid(a1, c1, a2, c2 float64) MF` | DiffSigmoid returns the difference-of-sigmoids membership function. |
| `func DilateMF(mf MF) MF` | DilateMF returns the dilation hedge mf(x)^0.5, spreading a fuzzy set (the linguistic hedge "somewhat" or "more or less"). |
| `func ExtremelyMF(mf MF) MF` | ExtremelyMF returns the "extremely" hedge, mf(x)^3. |
| `func Gaussian(mean, sigma float64) MF` | Gaussian returns the Gaussian membership function with center mean and width sigma. |
| `func Gaussian2(mean1, sigma1, mean2, sigma2 float64) MF` | Gaussian2 returns the two-sided Gaussian membership function combining a left flank (mean1, sigma1) and a right flank (mean2, sigma2). |
| `func IntensifyMF(mf MF) MF` | IntensifyMF returns the contrast intensification hedge that maps grades below 0.5 to 2*mf^2 and grades at or above 0.5 to 1 - 2*(1-mf)^2, pushing… |
| `func IntersectionMF(a, b MF) MF` | IntersectionMF returns the intersection (min) of two membership functions. |
| `func IntersectionMFWith(a, b MF, t TNorm) MF` | IntersectionMFWith returns the intersection of two membership functions using the t-norm t as the conjunction operator. |
| `func LeftShoulder(a, b float64) MF` | LeftShoulder returns a membership function that is 1 up to a, falls linearly to 0 at b and remains 0 thereafter. |
| `func MoreOrLessMF(mf MF) MF` | MoreOrLessMF returns the "more or less" hedge, mf(x)^0.5, a synonym of SomewhatMF. |
| `func NotMF(mf MF) MF` | NotMF returns the "not" hedge, the standard complement 1 - mf(x). |
| `func PiShape(a, b, c, d float64) MF` | PiShape returns the Pi-shaped membership function with parameters a <= b <= c <= d built from an S-shaped rising edge and a Z-shaped falling edge. |
| `func PowMF(mf MF, p float64) MF` | PowMF returns the membership function mf raised pointwise to the power p, mf(x)^p. |
| `func ProdSigmoid(a1, c1, a2, c2 float64) MF` | ProdSigmoid returns the product-of-sigmoids membership function. |
| `func RampDown(a, b float64) MF` | RampDown returns the linear falling ramp membership function from a to b. |
| `func RampUp(a, b float64) MF` | RampUp returns the linear rising ramp membership function from a to b. |
| `func Rectangular(a, b float64) MF` | Rectangular returns the crisp rectangular membership function equal to 1 on [a, b] and 0 outside. |
| `func RightShoulder(a, b float64) MF` | RightShoulder returns a membership function that is 0 up to a, rises linearly to 1 at b and remains 1 thereafter. |
| `func SShape(a, b float64) MF` | SShape returns the spline-based S-shaped membership function with feet a and b. |
| `func ScaleMF(mf MF, k float64) MF` | ScaleMF returns the membership function mf multiplied pointwise by the factor k and clamped to [0, 1]. |
| `func Sigmoid(a, c float64) MF` | Sigmoid returns the sigmoidal membership function with slope a and inflection point c. |
| `func Singleton(c float64) MF` | Singleton returns the fuzzy singleton membership function concentrated at c. |
| `func SlightlyMF(mf MF) MF` | SlightlyMF returns the "slightly" hedge, mf(x)^(1/3). |
| `func SomewhatMF(mf MF) MF` | SomewhatMF returns the "somewhat" hedge, mf(x)^0.5. |
| `func Trapezoidal(a, b, c, d float64) MF` | Trapezoidal returns the trapezoidal membership function with feet a and d and shoulders b and c. |
| `func Triangular(a, b, c float64) MF` | Triangular returns the triangular membership function with feet a and c and peak b. |
| `func UnionMF(a, b MF) MF` | UnionMF returns the union (max) of two membership functions. |
| `func UnionMFWith(a, b MF, s TConorm) MF` | UnionMFWith returns the union of two membership functions using the t-conorm s as the disjunction operator. |
| `func VeryMF(mf MF) MF` | VeryMF returns the "very" hedge, mf(x)^2. |
| `func ZShape(a, b float64) MF` | ZShape returns the spline-based Z-shaped membership function with shoulders a and b. |

</details>

<details>
<summary><code>MamdaniSystem</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMamdani(output *Variable) *MamdaniSystem` | NewMamdani creates a Mamdani system with the given output variable and the standard operator set (min t-norm, max t-conorm, standard complement, min… |
| `func (m *MamdaniSystem) AddInput(v *Variable) *MamdaniSystem` | AddInput registers an input variable and returns the system for chaining. |
| `func (m *MamdaniSystem) AddRule(op RuleOp, outTerm string, weight float64, clauses ...Clause) *MamdaniSystem` | AddRule appends a Mamdani rule and returns the system for chaining. |
| `func (m *MamdaniSystem) AggregateOutput(inputs map[string]float64) (Set, error)` | AggregateOutput returns the aggregated output fuzzy set produced by firing all rules against the crisp inputs, sampled on the output universe at… |
| `func (m *MamdaniSystem) FiringStrength(rule MamdaniRule, inputs map[string]float64) (float64, error)` | FiringStrength returns the degree to which the rule's antecedent is satisfied by the crisp inputs, combining its clauses with the system's And or Or… |
| `func (m *MamdaniSystem) Infer(inputs map[string]float64) (float64, error)` | Infer runs the full Mamdani pipeline (fuzzify, apply rules, aggregate, defuzzify) and returns the crisp output. |

</details>

<details>
<summary><code>Relation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CartesianProduct(a, b Set, tn TNorm) Relation` | CartesianProduct returns the fuzzy relation a x b formed by the t-norm tn, with grade tn(Mu_a(X`i`), Mu_b(Y`j`)). |
| `func CylindricalExtension(a Set, y []float64) Relation` | CylindricalExtension returns the cylindrical extension of the fuzzy set a (defined on X) onto the product X x y, with grade r`i``j` = Mu_a(X`i`) for… |
| `func NewRelation(x, y []float64, m [][]float64) (Relation, error)` | NewRelation builds a fuzzy relation from universes x and y and a grade matrix m of shape len(x) by len(y). |
| `func RelationFromFunc(x, y []float64, f func(a, b float64) float64) Relation` | RelationFromFunc builds a fuzzy relation whose grade at (X`i`, Y`j`) is f(X`i`, Y`j`), clamped to [0, 1]. |
| `func (r Relation) At(i, j int) float64` | At returns the grade of the (i, j) cell. |
| `func (r Relation) Clone() Relation` | Clone returns a deep copy of the relation. |
| `func (r Relation) Cols() int` | Cols reports the number of columns (size of universe Y). |
| `func (r Relation) Complement() Relation` | Complement returns the standard complement 1 - M of the relation. |
| `func (r Relation) ComposeSet(a Set) (Set, error)` | ComposeSet returns the fuzzy set on Y obtained by the max-min composition of the fuzzy set a (on X) with the relation r, grade_k = max_i… |
| `func (r Relation) ComposeSetT(a Set, tn TNorm) (Set, error)` | ComposeSetT returns the sup-t composition of the fuzzy set a with the relation r using the t-norm tn, grade_k = max_i tn(Mu_a(X`i`), r`i``k`). |
| `func (r Relation) Intersection(o Relation) (Relation, error)` | Intersection returns the elementwise intersection (min) of two relations of equal shape. |
| `func (r Relation) IsFuzzyEquivalence(tol float64) bool` | IsFuzzyEquivalence reports whether the relation is a fuzzy equivalence (similarity) relation within tol: reflexive, symmetric and max-min transitive. |
| `func (r Relation) IsMaxMinTransitive(tol float64) bool` | IsMaxMinTransitive reports whether the relation is max-min transitive within tol, that is r`i``k` >= max_j min(r`i``j`, r`j``k`) for all i, k. |
| `func (r Relation) IsReflexive(tol float64) bool` | IsReflexive reports whether the relation is reflexive within tol, that is whether it is square (X equals Y) and every diagonal grade is 1. |
| `func (r Relation) IsSymmetric(tol float64) bool` | IsSymmetric reports whether the relation is symmetric within tol. |
| `func (r Relation) MaxMinComposition(s Relation) (Relation, error)` | MaxMinComposition returns the max-min composition r o s, whose (i, k) grade is max_j min(r`i``j`, s`j``k`). |
| `func (r Relation) MaxProductComposition(s Relation) (Relation, error)` | MaxProductComposition returns the max-product composition r o s, whose (i, k) grade is max_j (r`i``j` * s`j``k`). |
| `func (r Relation) MaxTComposition(s Relation, tn TNorm) (Relation, error)` | MaxTComposition returns the sup-t composition r o s using the t-norm tn, whose (i, k) grade is max_j tn(r`i``j`, s`j``k`). |
| `func (r Relation) ProjectFirst() Set` | ProjectFirst returns the projection of the relation onto its X universe, the fuzzy set whose grade at X`i` is max_j r`i``j`. |
| `func (r Relation) ProjectSecond() Set` | ProjectSecond returns the projection of the relation onto its Y universe, the fuzzy set whose grade at Y`j` is max_i r`i``j`. |
| `func (r Relation) Rows() int` | Rows reports the number of rows (size of universe X). |
| `func (r Relation) TransitiveClosure() (Relation, error)` | TransitiveClosure returns the max-min transitive closure of a square relation, the union of r, r o r, r o r o r, ... |
| `func (r Relation) Transpose() Relation` | Transpose returns the inverse relation with X and Y swapped. |
| `func (r Relation) Union(o Relation) (Relation, error)` | Union returns the elementwise union (max) of two relations of equal shape. |

</details>

<details>
<summary><code>Set</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FromMF(mf MF, xs []float64) Set` | FromMF samples the membership function mf over the universe xs and returns the resulting discrete fuzzy set. |
| `func NewSet(x, mu []float64) (Set, error)` | NewSet builds a discrete fuzzy set from parallel universe and membership slices. |
| `func (s Set) AlphaCut(alpha float64) []float64` | AlphaCut returns the universe points whose grade is at least alpha, the (weak) alpha-cut of the set. |
| `func (s Set) Bisector() (float64, error)` | Bisector returns the bisector of area of the discrete fuzzy set, the universe point that splits the total membership sum into two (as equal as… |
| `func (s Set) Cardinality() float64` | Cardinality returns the scalar (sigma) cardinality of the set, the sum of all membership grades. |
| `func (s Set) Centroid() (float64, error)` | Centroid returns the center of gravity of the discrete fuzzy set, sum(X`i`*Mu`i`) / sum(Mu`i`). |
| `func (s Set) CentroidTrapz() (float64, error)` | CentroidTrapz returns the continuous centroid of the fuzzy set treating the membership as piecewise linear between universe points and integrating… |
| `func (s Set) Clone() Set` | Clone returns a deep copy of the set. |
| `func (s Set) Complement() Set` | Complement returns the standard complement 1 - Mu of the set on the same universe. |
| `func (s Set) ComplementWith(c func(float64) float64) Set` | ComplementWith returns the complement of the set using the fuzzy complement operator c. |
| `func (s Set) Concentrate() Set` | Concentrate returns the concentration hedge Mu^2 ("very"). |
| `func (s Set) Core() []float64` | Core returns the universe points whose grade equals 1, the crisp core of the set. |
| `func (s Set) Defuzzify(method Defuzz) (float64, error)` | Defuzzify applies the named defuzzification strategy to the set and returns the crisp result. |
| `func (s Set) DegreeOfSubsethood(t Set) (float64, error)` | DegreeOfSubsethood returns the Kosko degree to which s is a subset of t, cardinality(s AND t) / cardinality(s). |
| `func (s Set) Difference(t Set) (Set, error)` | Difference returns the bounded fuzzy difference s AND NOT t (min of s and the complement of t) on a shared universe. |
| `func (s Set) Dilate() Set` | Dilate returns the dilation hedge Mu^0.5 ("somewhat"). |
| `func (s Set) Equal(t Set, tol float64) bool` | Equal reports whether two sets share the same universe and have membership grades that agree within tol. |
| `func (s Set) Extremely() Set` | Extremely returns the "extremely" hedge, Mu^3. |
| `func (s Set) Height() float64` | Height returns the largest membership grade in the set, 0 for an empty set. |
| `func (s Set) Intensify() Set` | Intensify returns the contrast intensification hedge applied to every grade. |
| `func (s Set) Intersection(t Set) (Set, error)` | Intersection returns the intersection (min) of two sets defined on the same universe. |
| `func (s Set) IntersectionWith(t Set, tn TNorm) (Set, error)` | IntersectionWith returns the intersection of two sets using the t-norm tn. |
| `func (s Set) IsConvex() bool` | IsConvex reports whether the set is fuzzy-convex, that is whether its membership sequence over the sorted universe is unimodal (never rises after it… |
| `func (s Set) IsEmpty() bool` | IsEmpty reports whether every membership grade is 0. |
| `func (s Set) IsNormal() bool` | IsNormal reports whether the set attains grade 1 at some point. |
| `func (s Set) IsSubset(t Set, tol float64) bool` | IsSubset reports whether s is a fuzzy subset of t, that is Mu_s(x) <= Mu_t(x) within tol at every shared universe point. |
| `func (s Set) LargestOfMaxima() (float64, error)` | LargestOfMaxima returns the largest of maxima (LOM), the greatest universe point at which the membership attains its maximum. |
| `func (s Set) Len() int` | Len reports the number of points in the set's universe. |
| `func (s Set) MeanOfMaxima() (float64, error)` | MeanOfMaxima returns the mean of maxima (MOM), the average of the universe points at which the membership attains its maximum. |
| `func (s Set) Membership(x float64) float64` | Membership returns the grade of x. |
| `func (s Set) MoreOrLess() Set` | MoreOrLess returns the "more or less" hedge, Mu^0.5. |
| `func (s Set) Normalize() Set` | Normalize returns the set rescaled so its height becomes 1, dividing every grade by the current height. |
| `func (s Set) Not() Set` | Not returns the standard complement, a synonym for Complement suited to linguistic negation. |
| `func (s Set) Pow(p float64) Set` | Pow returns the set with every grade raised to the power p, the basis of the concentration (p>1) and dilation (0<p<1) hedges. |
| `func (s Set) RelativeCardinality() float64` | RelativeCardinality returns the scalar cardinality divided by the number of universe points, the mean membership grade. |
| `func (s Set) Slightly() Set` | Slightly returns the "slightly" hedge, Mu^(1/3). |
| `func (s Set) SmallestOfMaxima() (float64, error)` | SmallestOfMaxima returns the smallest of maxima (SOM), the least universe point at which the membership attains its maximum. |
| `func (s Set) Somewhat() Set` | Somewhat returns the "somewhat" hedge, Mu^0.5. |
| `func (s Set) String() string` | String renders the set in Zadeh's notation, mu1/x1 + mu2/x2 + ..., listing only points with positive grade. |
| `func (s Set) StrongAlphaCut(alpha float64) []float64` | StrongAlphaCut returns the universe points whose grade is strictly greater than alpha, the strong alpha-cut of the set. |
| `func (s Set) Support() []float64` | Support returns the universe points whose grade is strictly positive, the crisp support of the set. |
| `func (s Set) Union(t Set) (Set, error)` | Union returns the union (max) of two sets defined on the same universe. |
| `func (s Set) UnionWith(t Set, sn TConorm) (Set, error)` | UnionWith returns the union of two sets using the t-conorm sn. |
| `func (s Set) Very() Set` | Very returns the "very" hedge, Mu^2. |
| `func (s Set) WeightedAverage() (float64, error)` | WeightedAverage returns the weighted average defuzzification sum(X`i`*Mu`i`) / sum(Mu`i`). |

</details>

<details>
<summary><code>SugenoConsequent</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (c SugenoConsequent) Eval(inputs map[string]float64) (float64, error)` | Eval computes the crisp consequent value for the given crisp inputs. |

</details>

<details>
<summary><code>SugenoSystem</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewSugeno() *SugenoSystem` | NewSugeno creates a Sugeno system with the standard operator set (min t-norm, max t-conorm and standard complement). |
| `func (s *SugenoSystem) AddConstRule(op RuleOp, z, weight float64, clauses ...Clause) *SugenoSystem` | AddConstRule appends a zeroth-order Sugeno rule whose consequent is the constant z and returns the system for chaining. |
| `func (s *SugenoSystem) AddInput(v *Variable) *SugenoSystem` | AddInput registers an input variable and returns the system for chaining. |
| `func (s *SugenoSystem) AddRule(op RuleOp, cons SugenoConsequent, weight float64, clauses ...Clause) *SugenoSystem` | AddRule appends a first-order Sugeno rule with an explicit consequent and returns the system for chaining. |
| `func (s *SugenoSystem) FiringStrength(rule SugenoRule, inputs map[string]float64) (float64, error)` | FiringStrength returns the degree to which the rule's antecedent is satisfied by the crisp inputs, combining its clauses with the system's And or Or… |
| `func (s *SugenoSystem) Infer(inputs map[string]float64) (float64, error)` | Infer runs the Sugeno pipeline and returns the crisp output, the firing-strength weighted average of the rule consequents, sum(w_i * z_i) / sum(w_i). |

</details>

<details>
<summary><code>TConorm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DualTConorm(t TNorm) TConorm` | DualTConorm returns the t-conorm dual to the t-norm t under the standard complement, s(a, b) = 1 - t(1-a, 1-b). |
| `func TConormHamacher(gamma float64) TConorm` | TConormHamacher returns the Hamacher sum t-conorm with parameter gamma >= 0: (a + b + (gamma-1)*a*b) / (1 + gamma*a*b). |

</details>

<details>
<summary><code>TNorm</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DualTNorm(s TConorm) TNorm` | DualTNorm returns the t-norm dual to the t-conorm s under the standard complement, t(a, b) = 1 - s(1-a, 1-b). |
| `func TNormHamacher(gamma float64) TNorm` | TNormHamacher returns the Hamacher product t-norm with parameter gamma >= 0: a*b / (gamma + (1-gamma)(a+b-a*b)). |

</details>

<details>
<summary><code>Variable</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewVariable(name string, min, max float64) *Variable` | NewVariable creates a linguistic variable named name over the universe [min, max] with no terms. |
| `func (v *Variable) AddTerm(term string, mf MF) *Variable` | AddTerm registers the linguistic term named term with membership function mf and returns the variable for chaining. |
| `func (v *Variable) Fuzzify(term string, x float64) (float64, error)` | Fuzzify returns the membership grade of the crisp value x in the term named term. |

</details>

### Variables

`ErrComposeShape`, `ErrDomainMismatch`, `ErrEmptySet`, `ErrLengthMismatch`, `ErrNoArea`, `ErrNoInput`, `ErrNoRules`, `ErrRelationShape`, `ErrUnknownTerm`, `ErrUnknownVariable`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/fuzzy).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
