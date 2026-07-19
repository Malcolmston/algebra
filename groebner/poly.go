package groebner

import (
	"math/big"
	"sort"
	"strings"
)

// Poly is a sparse multivariate polynomial with rational coefficients over a
// fixed number of variables. Internally the terms are kept in a canonical form:
// combined like terms, no zero coefficients, and sorted in descending
// lexicographic order of their monomials. This canonical form makes structural
// equality a simple term-by-term comparison; the leading term with respect to a
// particular monomial order is computed on demand.
type Poly struct {
	nvars int
	terms []Term
}

// NewPoly builds a polynomial in n variables from the given terms, reducing it
// to canonical form. Terms with mismatched arity are padded or truncated to n
// variables.
func NewPoly(n int, terms ...Term) Poly {
	p := Poly{nvars: n}
	for _, t := range terms {
		p.terms = append(p.terms, Term{Coeff: CloneRat(t.Coeff), Mono: resize(t.Mono, n)})
	}
	p.canonicalize()
	return p
}

func resize(m Monomial, n int) Monomial {
	r := make(Monomial, n)
	for i := 0; i < n && i < len(m); i++ {
		r[i] = m[i]
	}
	return r
}

// Zero returns the zero polynomial in n variables.
func Zero(n int) Poly { return Poly{nvars: n} }

// One returns the constant polynomial 1 in n variables.
func One(n int) Poly { return Constant(n, bigOne) }

// Constant returns the constant polynomial equal to c in n variables.
func Constant(n int, c *big.Rat) Poly {
	if c.Sign() == 0 {
		return Zero(n)
	}
	return Poly{nvars: n, terms: []Term{{Coeff: CloneRat(c), Mono: ZeroMonomial(n)}}}
}

// ConstantInt returns the constant polynomial equal to the integer c.
func ConstantInt(n int, c int64) Poly { return Constant(n, RatFromInt(c)) }

// Var returns the polynomial equal to the i-th variable x_i in n variables.
func Var(n, i int) Poly {
	return Poly{nvars: n, terms: []Term{{Coeff: CloneRat(bigOne), Mono: VarMonomial(n, i)}}}
}

// Vars returns the polynomials x_1, ..., x_n in a ring with n variables.
func Vars(n int) []Poly {
	v := make([]Poly, n)
	for i := range v {
		v[i] = Var(n, i)
	}
	return v
}

// Monomial returns the polynomial consisting of a single monomial m with
// coefficient 1.
func MonomialPoly(n int, m Monomial) Poly {
	return Poly{nvars: n, terms: []Term{{Coeff: CloneRat(bigOne), Mono: resize(m, n)}}}
}

// FromTerm returns the one-term polynomial equal to t.
func FromTerm(n int, t Term) Poly { return NewPoly(n, t) }

// canonicalize sorts, combines like terms and drops zeros.
func (p *Poly) canonicalize() {
	if len(p.terms) == 0 {
		return
	}
	sort.SliceStable(p.terms, func(i, j int) bool {
		return CompareLex(p.terms[i].Mono, p.terms[j].Mono) > 0
	})
	out := p.terms[:0]
	var cur Term
	have := false
	for _, t := range p.terms {
		if have && cur.Mono.Equal(t.Mono) {
			cur.Coeff = ratAdd(cur.Coeff, t.Coeff)
			continue
		}
		if have && cur.Coeff.Sign() != 0 {
			out = append(out, cur)
		}
		cur = Term{Coeff: CloneRat(t.Coeff), Mono: t.Mono.Clone()}
		have = true
	}
	if have && cur.Coeff.Sign() != 0 {
		out = append(out, cur)
	}
	p.terms = out
}

// Nvars returns the number of variables of the ambient ring.
func (p Poly) Nvars() int { return p.nvars }

// Len returns the number of (nonzero) terms in the polynomial.
func (p Poly) Len() int { return len(p.terms) }

// IsZero reports whether the polynomial is identically zero.
func (p Poly) IsZero() bool { return len(p.terms) == 0 }

// IsConstant reports whether the polynomial is a constant (zero or a single
// term of degree 0).
func (p Poly) IsConstant() bool {
	return len(p.terms) == 0 || (len(p.terms) == 1 && p.terms[0].Mono.IsConstant())
}

// IsOne reports whether the polynomial equals the constant 1.
func (p Poly) IsOne() bool {
	return len(p.terms) == 1 && p.terms[0].Mono.IsConstant() && RatIsOne(p.terms[0].Coeff)
}

// Clone returns an independent deep copy of the polynomial.
func (p Poly) Clone() Poly {
	c := Poly{nvars: p.nvars, terms: make([]Term, len(p.terms))}
	for i, t := range p.terms {
		c.terms[i] = t.Clone()
	}
	return c
}

// Terms returns a copy of the polynomial's terms in canonical (descending
// lexicographic) order.
func (p Poly) Terms() []Term {
	out := make([]Term, len(p.terms))
	for i, t := range p.terms {
		out[i] = t.Clone()
	}
	return out
}

// TotalDegree returns the maximum total degree among the polynomial's terms, or
// -1 for the zero polynomial.
func (p Poly) TotalDegree() int {
	d := -1
	for _, t := range p.terms {
		if td := t.Mono.Degree(); td > d {
			d = td
		}
	}
	return d
}

// Coeff returns the coefficient of the monomial m in the polynomial, or 0 if it
// does not appear.
func (p Poly) Coeff(m Monomial) *big.Rat {
	for _, t := range p.terms {
		if t.Mono.Equal(m) {
			return CloneRat(t.Coeff)
		}
	}
	return big.NewRat(0, 1)
}

// ConstantTerm returns the constant coefficient of the polynomial.
func (p Poly) ConstantTerm() *big.Rat {
	return p.Coeff(ZeroMonomial(p.nvars))
}

// Monomials returns the monomials of the polynomial in canonical order.
func (p Poly) Monomials() []Monomial {
	out := make([]Monomial, len(p.terms))
	for i, t := range p.terms {
		out[i] = t.Mono.Clone()
	}
	return out
}

// Coefficients returns the coefficients of the polynomial in canonical order.
func (p Poly) Coefficients() []*big.Rat {
	out := make([]*big.Rat, len(p.terms))
	for i, t := range p.terms {
		out[i] = CloneRat(t.Coeff)
	}
	return out
}

// LeadingTerm returns the term of the polynomial whose monomial is largest with
// respect to the order o. It returns the zero term for the zero polynomial.
func (p Poly) LeadingTerm(o Order) Term {
	if len(p.terms) == 0 {
		return Term{Coeff: big.NewRat(0, 1), Mono: ZeroMonomial(p.nvars)}
	}
	best := 0
	for i := 1; i < len(p.terms); i++ {
		if o(p.terms[i].Mono, p.terms[best].Mono) > 0 {
			best = i
		}
	}
	return p.terms[best].Clone()
}

// LeadingMonomial returns the monomial of the leading term with respect to o.
func (p Poly) LeadingMonomial(o Order) Monomial { return p.LeadingTerm(o).Mono }

// LeadingCoeff returns the coefficient of the leading term with respect to o.
func (p Poly) LeadingCoeff(o Order) *big.Rat { return p.LeadingTerm(o).Coeff }

// Multidegree returns the exponent vector of the leading monomial with respect
// to o.
func (p Poly) Multidegree(o Order) Monomial { return p.LeadingMonomial(o) }

// Equal reports whether two polynomials are equal as elements of the ring.
func (p Poly) Equal(o Poly) bool {
	if p.nvars != o.nvars || len(p.terms) != len(o.terms) {
		return false
	}
	for i := range p.terms {
		if !p.terms[i].Equal(o.terms[i]) {
			return false
		}
	}
	return true
}

// String renders the polynomial with default variable names x1, x2, ....
func (p Poly) String() string { return p.Format(defaultVarNames(p.nvars)) }

// Format renders the polynomial with the supplied variable names, terms in
// descending lexicographic order. The zero polynomial renders as "0".
func (p Poly) Format(vars []string) string {
	if len(p.terms) == 0 {
		return "0"
	}
	var b strings.Builder
	for i, t := range p.terms {
		s := t.Format(vars)
		if i == 0 {
			b.WriteString(s)
			continue
		}
		if strings.HasPrefix(s, "-") {
			b.WriteString(" - ")
			b.WriteString(s[1:])
		} else {
			b.WriteString(" + ")
			b.WriteString(s)
		}
	}
	return b.String()
}
