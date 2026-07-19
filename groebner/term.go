package groebner

import (
	"math/big"
	"strings"
)

// Term is a single monomial with a rational coefficient, the building block of
// a polynomial.
type Term struct {
	Coeff *big.Rat
	Mono  Monomial
}

// NewTerm returns a term with the given coefficient and monomial. The
// coefficient is copied; the monomial is cloned.
func NewTerm(c *big.Rat, m Monomial) Term {
	return Term{Coeff: CloneRat(c), Mono: m.Clone()}
}

// Clone returns an independent copy of the term.
func (t Term) Clone() Term {
	return Term{Coeff: CloneRat(t.Coeff), Mono: t.Mono.Clone()}
}

// IsZero reports whether the term's coefficient is zero.
func (t Term) IsZero() bool { return t.Coeff.Sign() == 0 }

// Degree returns the total degree of the term's monomial.
func (t Term) Degree() int { return t.Mono.Degree() }

// Neg returns the term with negated coefficient.
func (t Term) Neg() Term {
	return Term{Coeff: ratNeg(t.Coeff), Mono: t.Mono.Clone()}
}

// Mul returns the product of two terms: coefficients multiply and monomials
// combine by adding exponents.
func (t Term) Mul(o Term) Term {
	return Term{Coeff: ratMul(t.Coeff, o.Coeff), Mono: t.Mono.Mul(o.Mono)}
}

// Equal reports whether two terms have equal coefficients and equal monomials.
func (t Term) Equal(o Term) bool {
	return t.Coeff.Cmp(o.Coeff) == 0 && t.Mono.Equal(o.Mono)
}

// Divides reports whether term t divides term o, i.e. its monomial divides the
// monomial of o (coefficients are units in a field so are ignored).
func (t Term) Divides(o Term) bool { return t.Mono.Divides(o.Mono) }

// Div returns the quotient o/t as a term when the monomial of t divides that of
// o. The boolean result is false when the monomial division is not exact.
func (t Term) Div(o Term) (Term, bool) {
	m, ok := t.Mono.Div(o.Mono)
	if !ok {
		return Term{}, false
	}
	return Term{Coeff: ratDiv(o.Coeff, t.Coeff), Mono: m}, true
}

// String renders the term with default variable names.
func (t Term) String() string { return t.Format(defaultVarNames(len(t.Mono))) }

// Format renders the term with the supplied variable names.
func (t Term) Format(vars []string) string {
	if t.Mono.IsConstant() {
		return t.Coeff.RatString()
	}
	var b strings.Builder
	switch {
	case t.Coeff.Cmp(bigOne) == 0:
		// omit coefficient 1
	case t.Coeff.Cmp(big.NewRat(-1, 1)) == 0:
		b.WriteString("-")
	default:
		b.WriteString(t.Coeff.RatString())
		b.WriteString("*")
	}
	b.WriteString(t.Mono.Format(vars))
	return b.String()
}
