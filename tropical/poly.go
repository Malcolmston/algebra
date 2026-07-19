package tropical

import (
	"math"
	"strconv"
	"strings"
)

// Poly is a univariate tropical polynomial. Coefficient i multiplies x raised
// to the tropical power i, so evaluation is the tropical sum over i of
// coeff[i] (*) x^i, that is min_i(coeff[i] + i*x) for min-plus and the
// corresponding max for max-plus. Trailing tropical-zero coefficients are
// trimmed on construction.
type Poly struct {
	coeffs []float64
	sr     Semiring
}

// NewPoly returns the tropical polynomial over sr with the given coefficients,
// coeffs[i] belonging to x^i. It stores a trimmed copy.
func NewPoly(sr Semiring, coeffs []float64) Poly {
	c := make([]float64, len(coeffs))
	copy(c, coeffs)
	return Poly{coeffs: trimPoly(c, sr), sr: sr}
}

// MinPlusPoly returns a min-plus tropical polynomial with the given
// coefficients.
func MinPlusPoly(coeffs []float64) Poly { return NewPoly(MinPlusSemiring(), coeffs) }

// MaxPlusPoly returns a max-plus tropical polynomial with the given
// coefficients.
func MaxPlusPoly(coeffs []float64) Poly { return NewPoly(MaxPlusSemiring(), coeffs) }

func trimPoly(c []float64, sr Semiring) []float64 {
	z := sr.Zero()
	end := len(c)
	for end > 1 && c[end-1] == z {
		end--
	}
	return c[:end]
}

// Degree returns the degree of the polynomial: the highest power with a
// non-zero coefficient. The tropical-zero polynomial has degree 0.
func (p Poly) Degree() int { return len(p.coeffs) - 1 }

// Semiring returns the semiring under which the polynomial is interpreted.
func (p Poly) Semiring() Semiring { return p.sr }

// Coeff returns the coefficient of x^i, or the tropical zero when i is out of
// range.
func (p Poly) Coeff(i int) float64 {
	if i < 0 || i >= len(p.coeffs) {
		return p.sr.Zero()
	}
	return p.coeffs[i]
}

// Coeffs returns a fresh copy of the coefficient slice, index i holding the
// coefficient of x^i.
func (p Poly) Coeffs() []float64 {
	c := make([]float64, len(p.coeffs))
	copy(c, p.coeffs)
	return c
}

// Clone returns an independent copy of the polynomial.
func (p Poly) Clone() Poly { return NewPoly(p.sr, p.coeffs) }

// Equal reports whether p and q are over the same semiring and have identical
// coefficients.
func (p Poly) Equal(q Poly) bool {
	if p.sr.kind != q.sr.kind || len(p.coeffs) != len(q.coeffs) {
		return false
	}
	for i := range p.coeffs {
		if p.coeffs[i] != q.coeffs[i] {
			return false
		}
	}
	return true
}

// Eval returns the value of the tropical polynomial at x.
func (p Poly) Eval(x float64) float64 {
	r := p.sr.Zero()
	for i, a := range p.coeffs {
		r = p.sr.Add(r, p.sr.Mul(a, p.sr.Pow(x, i)))
	}
	return r
}

// EvalArg returns the value of the polynomial at x together with the index of a
// term attaining it (the active monomial). The index is -1 for the tropical
// zero polynomial.
func (p Poly) EvalArg(x float64) (float64, int) {
	best := p.sr.Zero()
	idx := -1
	for i, a := range p.coeffs {
		t := p.sr.Mul(a, p.sr.Pow(x, i))
		if idx == -1 || (t != best && p.sr.AtLeastAsGood(t, best)) {
			best = t
			idx = i
		}
	}
	return best, idx
}

// EvalAll evaluates the polynomial at every x in xs and returns the values.
func (p Poly) EvalAll(xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = p.Eval(x)
	}
	return out
}

// Add returns the tropical sum of p and q, computed coefficientwise.
func (p Poly) Add(q Poly) Poly {
	n := len(p.coeffs)
	if len(q.coeffs) > n {
		n = len(q.coeffs)
	}
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = p.sr.Add(p.Coeff(i), q.Coeff(i))
	}
	return NewPoly(p.sr, out)
}

// Mul returns the tropical product of p and q, computed as a tropical
// convolution of their coefficients.
func (p Poly) Mul(q Poly) Poly {
	if len(p.coeffs) == 0 || len(q.coeffs) == 0 {
		return NewPoly(p.sr, []float64{p.sr.Zero()})
	}
	n := len(p.coeffs) + len(q.coeffs) - 1
	out := make([]float64, n)
	for i := range out {
		out[i] = p.sr.Zero()
	}
	for i, a := range p.coeffs {
		for j, b := range q.coeffs {
			out[i+j] = p.sr.Add(out[i+j], p.sr.Mul(a, b))
		}
	}
	return NewPoly(p.sr, out)
}

// Pow returns the tropical power p^n. The zeroth power is the tropical one
// polynomial. It panics for negative n.
func (p Poly) Pow(n int) Poly {
	if n < 0 {
		panic("tropical: Poly.Pow requires n >= 0")
	}
	result := NewPoly(p.sr, []float64{p.sr.One()})
	base := p.Clone()
	for n > 0 {
		if n&1 == 1 {
			result = result.Mul(base)
		}
		n >>= 1
		if n > 0 {
			base = base.Mul(base)
		}
	}
	return result
}

// ScalarMul returns the polynomial with every coefficient tropically multiplied
// by c.
func (p Poly) ScalarMul(c float64) Poly {
	out := make([]float64, len(p.coeffs))
	for i, a := range p.coeffs {
		out[i] = p.sr.Mul(a, c)
	}
	return NewPoly(p.sr, out)
}

// TropicallyVanishes reports whether the tropical polynomial vanishes at x,
// i.e. whether the optimal monomial value is attained by at least two distinct
// terms to within tol. These points are exactly the tropical roots.
func (p Poly) TropicallyVanishes(x, tol float64) bool {
	best := p.sr.Zero()
	for i, a := range p.coeffs {
		t := p.sr.Mul(a, p.sr.Pow(x, i))
		best = p.sr.Add(best, t)
	}
	if p.sr.IsZero(best) {
		return false
	}
	count := 0
	for i, a := range p.coeffs {
		t := p.sr.Mul(a, p.sr.Pow(x, i))
		if math.Abs(t-best) <= tol {
			count++
		}
	}
	return count >= 2
}

// IsRoot reports whether x is a tropical root of the polynomial, i.e. whether
// the polynomial tropically vanishes at x to within tol.
func (p Poly) IsRoot(x, tol float64) bool { return p.TropicallyVanishes(x, tol) }

// PolygonVertex is a vertex of a Newton polygon: a coefficient index paired
// with its (finite) coefficient value.
type PolygonVertex struct {
	Index int
	Coeff float64
}

// NewtonPolygon returns the vertices of the Newton polygon of the polynomial:
// the lower convex hull of the finite coefficient points (i, coeff[i]) for
// min-plus and the upper convex hull for max-plus, ordered by increasing index.
// The slopes of its edges are the negatives of the tropical roots and the edge
// widths are the root multiplicities.
func (p Poly) NewtonPolygon() []PolygonVertex {
	pts := make([]PolygonVertex, 0, len(p.coeffs))
	z := p.sr.Zero()
	for i, a := range p.coeffs {
		if a != z {
			pts = append(pts, PolygonVertex{Index: i, Coeff: a})
		}
	}
	if len(pts) <= 1 {
		return pts
	}
	lower := p.sr.IsMinPlus()
	var h []PolygonVertex
	for _, q := range pts {
		for len(h) >= 2 {
			o := h[len(h)-2]
			b := h[len(h)-1]
			cr := float64(b.Index-o.Index)*(q.Coeff-o.Coeff) - (b.Coeff-o.Coeff)*float64(q.Index-o.Index)
			if (lower && cr <= 0) || (!lower && cr >= 0) {
				h = h[:len(h)-1]
			} else {
				break
			}
		}
		h = append(h, q)
	}
	return h
}

// Root is a tropical root together with its multiplicity (the horizontal width
// of the Newton-polygon edge that produced it).
type Root struct {
	Value        float64
	Multiplicity int
}

// Roots returns the tropical roots of the polynomial with their multiplicities,
// read off the Newton polygon and ordered by increasing root value. A
// polynomial of degree d whose constant and leading coefficients are both
// finite has multiplicities summing to d. Constant and monomial polynomials
// have no roots.
func (p Poly) Roots() []Root {
	hull := p.NewtonPolygon()
	if len(hull) < 2 {
		return nil
	}
	roots := make([]Root, 0, len(hull)-1)
	for k := 0; k+1 < len(hull); k++ {
		a := hull[k]
		b := hull[k+1]
		width := b.Index - a.Index
		val := (a.Coeff - b.Coeff) / float64(width)
		roots = append(roots, Root{Value: val, Multiplicity: width})
	}
	sortRootsByValue(roots)
	return roots
}

// RootValues returns the distinct tropical roots repeated according to their
// multiplicity, ordered by increasing value.
func (p Poly) RootValues() []float64 {
	roots := p.Roots()
	var out []float64
	for _, r := range roots {
		for k := 0; k < r.Multiplicity; k++ {
			out = append(out, r.Value)
		}
	}
	return out
}

// String renders the polynomial in the form "a0 (+) a1 x (+) a2 x^2 (+) ...",
// omitting terms with a tropical-zero coefficient.
func (p Poly) String() string {
	z := p.sr.Zero()
	var parts []string
	for i, a := range p.coeffs {
		if a == z {
			continue
		}
		term := p.sr.FormatScalar(a)
		switch i {
		case 0:
			// constant term, keep as is
		case 1:
			term += " x"
		default:
			term += " x^" + strconv.Itoa(i)
		}
		parts = append(parts, term)
	}
	if len(parts) == 0 {
		return p.sr.FormatScalar(z)
	}
	return strings.Join(parts, " (+) ")
}

// FromRoots returns the monic tropical polynomial over sr whose tropical roots
// are exactly the given values, that is the tropical product of the linear
// factors (x (+) r) over all r. The leading coefficient is the tropical one.
func FromRoots(sr Semiring, roots []float64) Poly {
	result := NewPoly(sr, []float64{sr.One()})
	for _, r := range roots {
		factor := NewPoly(sr, []float64{r, sr.One()})
		result = result.Mul(factor)
	}
	return result
}

// MinPlusFromRoots returns the monic min-plus polynomial with the given roots.
func MinPlusFromRoots(roots []float64) Poly { return FromRoots(MinPlusSemiring(), roots) }

// MaxPlusFromRoots returns the monic max-plus polynomial with the given roots.
func MaxPlusFromRoots(roots []float64) Poly { return FromRoots(MaxPlusSemiring(), roots) }

func sortRootsByValue(r []Root) {
	for i := 1; i < len(r); i++ {
		for j := i; j > 0 && r[j-1].Value > r[j].Value; j-- {
			r[j-1], r[j] = r[j], r[j-1]
		}
	}
}
