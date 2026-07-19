package exterioralgebra

import (
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"
)

// Poly is a multivariate polynomial in n real variables x_0,…,x_{n-1}. It is
// stored sparsely as a map from exponent vectors to nonzero coefficients, so it
// can represent polynomials of any degree without a dense coefficient array.
// Poly is the coefficient ring used by [PForm] for differential forms with
// polynomial coefficients.
//
// The zero value is not usable; build polynomials with [NewPoly], [ConstPoly],
// [Var] and [Monomial].
type Poly struct {
	n     int
	terms map[string]polyTerm
}

type polyTerm struct {
	exp   []int
	coeff float64
}

// encodeExp turns an exponent vector into a stable string map key.
func encodeExp(exp []int) string {
	var b strings.Builder
	for i, e := range exp {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(strconv.Itoa(e))
	}
	return b.String()
}

// NewPoly returns the zero polynomial in n variables. It panics if n < 0.
func NewPoly(n int) *Poly {
	if n < 0 {
		panic(ErrDim)
	}
	return &Poly{n: n, terms: make(map[string]polyTerm)}
}

// ConstPoly returns the constant polynomial c in n variables.
func ConstPoly(n int, c float64) *Poly {
	p := NewPoly(n)
	if c != 0 {
		exp := make([]int, n)
		p.terms[encodeExp(exp)] = polyTerm{exp: exp, coeff: c}
	}
	return p
}

// Var returns the polynomial equal to the single variable x_i in n variables.
// It panics if i is out of range.
func Var(n, i int) *Poly {
	if i < 0 || i >= n {
		panic(ErrIndex)
	}
	p := NewPoly(n)
	exp := make([]int, n)
	exp[i] = 1
	p.terms[encodeExp(exp)] = polyTerm{exp: exp, coeff: 1}
	return p
}

// Monomial returns the polynomial coeff·∏ x_i^{exp[i]} in n = len(exp)
// variables.
func Monomial(coeff float64, exp []int) *Poly {
	p := NewPoly(len(exp))
	if coeff != 0 {
		e := append([]int(nil), exp...)
		p.terms[encodeExp(e)] = polyTerm{exp: e, coeff: coeff}
	}
	return p
}

// Arity returns the number of variables n of p.
func (p *Poly) Arity() int { return p.n }

// IsZero reports whether p is the zero polynomial.
func (p *Poly) IsZero() bool { return len(p.terms) == 0 }

// Clone returns an independent deep copy of p.
func (p *Poly) Clone() *Poly {
	q := NewPoly(p.n)
	for k, t := range p.terms {
		q.terms[k] = polyTerm{exp: append([]int(nil), t.exp...), coeff: t.coeff}
	}
	return q
}

// add accumulates coeff onto the monomial with the given exponent vector,
// pruning entries whose running coefficient reaches zero.
func (p *Poly) add(exp []int, coeff float64) {
	if coeff == 0 {
		return
	}
	k := encodeExp(exp)
	if t, ok := p.terms[k]; ok {
		v := t.coeff + coeff
		if v == 0 {
			delete(p.terms, k)
		} else {
			p.terms[k] = polyTerm{exp: t.exp, coeff: v}
		}
		return
	}
	p.terms[k] = polyTerm{exp: append([]int(nil), exp...), coeff: coeff}
}

// Add returns the sum p+q. It panics if the arities differ.
func (p *Poly) Add(q *Poly) *Poly {
	if p.n != q.n {
		panic(ErrDim)
	}
	r := p.Clone()
	for _, t := range q.terms {
		r.add(t.exp, t.coeff)
	}
	return r
}

// Sub returns the difference p−q. It panics if the arities differ.
func (p *Poly) Sub(q *Poly) *Poly {
	if p.n != q.n {
		panic(ErrDim)
	}
	r := p.Clone()
	for _, t := range q.terms {
		r.add(t.exp, -t.coeff)
	}
	return r
}

// Neg returns −p.
func (p *Poly) Neg() *Poly {
	r := NewPoly(p.n)
	for k, t := range p.terms {
		r.terms[k] = polyTerm{exp: t.exp, coeff: -t.coeff}
	}
	return r
}

// Scale returns c·p.
func (p *Poly) Scale(c float64) *Poly {
	r := NewPoly(p.n)
	if c == 0 {
		return r
	}
	for k, t := range p.terms {
		r.terms[k] = polyTerm{exp: t.exp, coeff: c * t.coeff}
	}
	return r
}

// Mul returns the product p·q. It panics if the arities differ.
func (p *Poly) Mul(q *Poly) *Poly {
	if p.n != q.n {
		panic(ErrDim)
	}
	r := NewPoly(p.n)
	for _, a := range p.terms {
		for _, b := range q.terms {
			exp := make([]int, p.n)
			for i := range exp {
				exp[i] = a.exp[i] + b.exp[i]
			}
			r.add(exp, a.coeff*b.coeff)
		}
	}
	return r
}

// Pow returns p raised to the non-negative integer power k. Pow(0) is the
// constant 1. It panics if k < 0.
func (p *Poly) Pow(k int) *Poly {
	if k < 0 {
		panic(ErrGrade)
	}
	res := ConstPoly(p.n, 1)
	base := p.Clone()
	for k > 0 {
		if k&1 == 1 {
			res = res.Mul(base)
		}
		k >>= 1
		if k > 0 {
			base = base.Mul(base)
		}
	}
	return res
}

// Partial returns the partial derivative ∂p/∂x_i. It panics if i is out of
// range.
func (p *Poly) Partial(i int) *Poly {
	if i < 0 || i >= p.n {
		panic(ErrIndex)
	}
	r := NewPoly(p.n)
	for _, t := range p.terms {
		if t.exp[i] == 0 {
			continue
		}
		exp := append([]int(nil), t.exp...)
		c := t.coeff * float64(exp[i])
		exp[i]--
		r.add(exp, c)
	}
	return r
}

// Gradient returns the n partial derivatives of p as a slice, with entry i
// equal to ∂p/∂x_i.
func (p *Poly) Gradient() []*Poly {
	g := make([]*Poly, p.n)
	for i := range g {
		g[i] = p.Partial(i)
	}
	return g
}

// Eval evaluates p at the point x, which must have length n. It panics on a
// length mismatch.
func (p *Poly) Eval(x []float64) float64 {
	if len(x) != p.n {
		panic(ErrDim)
	}
	var s float64
	for _, t := range p.terms {
		term := t.coeff
		for i, e := range t.exp {
			if e != 0 {
				term *= math.Pow(x[i], float64(e))
			}
		}
		s += term
	}
	return s
}

// Compose substitutes each variable x_i with the polynomial subs[i], returning
// the composed polynomial in the common arity of the substitutes. It requires
// len(subs) == n and that every substitute has the same arity; it returns
// [ErrMap] otherwise.
func (p *Poly) Compose(subs []*Poly) (*Poly, error) {
	if len(subs) != p.n {
		return nil, ErrMap
	}
	m := -1
	for _, s := range subs {
		if m == -1 {
			m = s.n
		} else if s.n != m {
			return nil, ErrMap
		}
	}
	if m == -1 {
		// No variables: p is a constant; report it in arity 0.
		return ConstPoly(0, p.ConstantTerm()), nil
	}
	res := NewPoly(m)
	for _, t := range p.terms {
		term := ConstPoly(m, t.coeff)
		for i, e := range t.exp {
			if e != 0 {
				term = term.Mul(subs[i].Pow(e))
			}
		}
		res = res.Add(term)
	}
	return res, nil
}

// ConstantTerm returns the coefficient of the constant monomial of p.
func (p *Poly) ConstantTerm() float64 {
	if p.n == 0 {
		for _, t := range p.terms {
			return t.coeff
		}
		return 0
	}
	return p.terms[encodeExp(make([]int, p.n))].coeff
}

// Degree returns the total degree of p, the largest sum of exponents over its
// monomials. The zero polynomial has degree −1 by convention.
func (p *Poly) Degree() int {
	deg := -1
	for _, t := range p.terms {
		s := 0
		for _, e := range t.exp {
			s += e
		}
		if s > deg {
			deg = s
		}
	}
	return deg
}

// NumTerms returns the number of nonzero monomials in p.
func (p *Poly) NumTerms() int { return len(p.terms) }

// Equal reports whether p and q are the same polynomial (same arity and
// coefficients).
func (p *Poly) Equal(q *Poly) bool {
	if p.n != q.n || len(p.terms) != len(q.terms) {
		return false
	}
	for k, t := range p.terms {
		u, ok := q.terms[k]
		if !ok || u.coeff != t.coeff {
			return false
		}
	}
	return true
}

// EqualTol reports whether p and q have the same arity and all coefficients
// agree to within an absolute tolerance tol.
func (p *Poly) EqualTol(q *Poly, tol float64) bool {
	if p.n != q.n {
		return false
	}
	for k, t := range p.terms {
		u := q.terms[k]
		if math.Abs(t.coeff-u.coeff) > tol {
			return false
		}
	}
	for k, t := range q.terms {
		if _, ok := p.terms[k]; !ok && math.Abs(t.coeff) > tol {
			return false
		}
	}
	return true
}

// sortedTerms returns the monomials of p ordered by descending total degree and
// then lexicographically by exponent vector.
func (p *Poly) sortedTerms() []polyTerm {
	out := make([]polyTerm, 0, len(p.terms))
	for _, t := range p.terms {
		out = append(out, t)
	}
	sort.Slice(out, func(i, j int) bool {
		di, dj := 0, 0
		for _, e := range out[i].exp {
			di += e
		}
		for _, e := range out[j].exp {
			dj += e
		}
		if di != dj {
			return di > dj
		}
		for k := range out[i].exp {
			if out[i].exp[k] != out[j].exp[k] {
				return out[i].exp[k] > out[j].exp[k]
			}
		}
		return false
	})
	return out
}

// String renders p in conventional notation such as "2 x0^2 + 3 x1 - 1".
func (p *Poly) String() string {
	if p.IsZero() {
		return "0"
	}
	var b strings.Builder
	for i, t := range p.sortedTerms() {
		c := t.coeff
		if i > 0 {
			if c < 0 {
				b.WriteString(" - ")
				c = -c
			} else {
				b.WriteString(" + ")
			}
		} else if c < 0 {
			b.WriteString("-")
			c = -c
		}
		var mono strings.Builder
		for j, e := range t.exp {
			if e == 0 {
				continue
			}
			if mono.Len() > 0 {
				mono.WriteByte(' ')
			}
			if e == 1 {
				fmt.Fprintf(&mono, "x%d", j)
			} else {
				fmt.Fprintf(&mono, "x%d^%d", j, e)
			}
		}
		if mono.Len() == 0 {
			fmt.Fprintf(&b, "%g", c)
		} else if c == 1 {
			b.WriteString(mono.String())
		} else {
			fmt.Fprintf(&b, "%g %s", c, mono.String())
		}
	}
	return b.String()
}
