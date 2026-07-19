package groebner

import (
	"math/big"
	"math/cmplx"
)

// Add returns the sum of two polynomials.
func (p Poly) Add(o Poly) Poly {
	r := Poly{nvars: p.nvars}
	r.terms = make([]Term, 0, len(p.terms)+len(o.terms))
	for _, t := range p.terms {
		r.terms = append(r.terms, t.Clone())
	}
	for _, t := range o.terms {
		r.terms = append(r.terms, t.Clone())
	}
	r.canonicalize()
	return r
}

// Sub returns the difference p - o.
func (p Poly) Sub(o Poly) Poly { return p.Add(o.Neg()) }

// Neg returns the polynomial with every coefficient negated.
func (p Poly) Neg() Poly {
	r := Poly{nvars: p.nvars, terms: make([]Term, len(p.terms))}
	for i, t := range p.terms {
		r.terms[i] = t.Neg()
	}
	return r
}

// ScalarMul returns the polynomial multiplied by the rational scalar c.
func (p Poly) ScalarMul(c *big.Rat) Poly {
	if c.Sign() == 0 {
		return Zero(p.nvars)
	}
	r := Poly{nvars: p.nvars, terms: make([]Term, len(p.terms))}
	for i, t := range p.terms {
		r.terms[i] = Term{Coeff: ratMul(t.Coeff, c), Mono: t.Mono.Clone()}
	}
	return r
}

// MulTerm returns the polynomial multiplied by the single term t.
func (p Poly) MulTerm(t Term) Poly {
	if t.IsZero() {
		return Zero(p.nvars)
	}
	r := Poly{nvars: p.nvars, terms: make([]Term, len(p.terms))}
	for i, s := range p.terms {
		r.terms[i] = s.Mul(t)
	}
	r.canonicalize()
	return r
}

// MulMonomial returns the polynomial multiplied by the monomial m (coefficient
// 1).
func (p Poly) MulMonomial(m Monomial) Poly {
	return p.MulTerm(Term{Coeff: CloneRat(bigOne), Mono: m})
}

// Mul returns the product of two polynomials.
func (p Poly) Mul(o Poly) Poly {
	if p.IsZero() || o.IsZero() {
		return Zero(p.nvars)
	}
	r := Poly{nvars: p.nvars}
	r.terms = make([]Term, 0, len(p.terms)*len(o.terms))
	for _, a := range p.terms {
		for _, b := range o.terms {
			r.terms = append(r.terms, a.Mul(b))
		}
	}
	r.canonicalize()
	return r
}

// Pow returns the polynomial raised to the non-negative integer power k, using
// binary exponentiation. Pow(0) returns the constant 1.
func (p Poly) Pow(k int) Poly {
	if k < 0 {
		return Zero(p.nvars)
	}
	result := One(p.nvars)
	base := p.Clone()
	for k > 0 {
		if k&1 == 1 {
			result = result.Mul(base)
		}
		k >>= 1
		if k > 0 {
			base = base.Mul(base)
		}
	}
	return result
}

// Monic returns the polynomial scaled so that its leading coefficient with
// respect to o is 1. The zero polynomial is returned unchanged.
func (p Poly) Monic(o Order) Poly {
	if p.IsZero() {
		return p.Clone()
	}
	lc := p.LeadingCoeff(o)
	return p.ScalarMul(ratDiv(bigOne, lc))
}

// Eval evaluates the polynomial at the rational point (one value per variable)
// and returns the exact rational result.
func (p Poly) Eval(point []*big.Rat) *big.Rat {
	acc := new(big.Rat)
	for _, t := range p.terms {
		v := CloneRat(t.Coeff)
		for i, e := range t.Mono {
			for k := 0; k < e; k++ {
				v = ratMul(v, point[i])
			}
		}
		acc.Add(acc, v)
	}
	return acc
}

// EvalComplex evaluates the polynomial at a complex point and returns the
// complex result. It is used by the numerical variety solver.
func (p Poly) EvalComplex(point []complex128) complex128 {
	var acc complex128
	for _, t := range p.terms {
		c, _ := t.Coeff.Float64()
		v := complex(c, 0)
		for i, e := range t.Mono {
			for k := 0; k < e; k++ {
				v *= point[i]
			}
		}
		acc += v
	}
	return acc
}

// Derivative returns the formal partial derivative of the polynomial with
// respect to the i-th variable.
func (p Poly) Derivative(i int) Poly {
	r := Poly{nvars: p.nvars}
	for _, t := range p.terms {
		e := t.Mono[i]
		if e == 0 {
			continue
		}
		m := t.Mono.Clone()
		m[i] = e - 1
		r.terms = append(r.terms, Term{Coeff: ratMul(t.Coeff, RatFromInt(int64(e))), Mono: m})
	}
	r.canonicalize()
	return r
}

// Subst substitutes the constant rational value val for the i-th variable and
// returns the resulting polynomial (still in n variables, with that variable
// eliminated from every term).
func (p Poly) Subst(i int, val *big.Rat) Poly {
	r := Poly{nvars: p.nvars}
	for _, t := range p.terms {
		e := t.Mono[i]
		c := CloneRat(t.Coeff)
		for k := 0; k < e; k++ {
			c = ratMul(c, val)
		}
		m := t.Mono.Clone()
		m[i] = 0
		r.terms = append(r.terms, Term{Coeff: c, Mono: m})
	}
	r.canonicalize()
	return r
}

// DependsOn reports whether the i-th variable appears with a positive exponent
// in some term of the polynomial.
func (p Poly) DependsOn(i int) bool {
	for _, t := range p.terms {
		if t.Mono[i] > 0 {
			return true
		}
	}
	return false
}

// UsedVars returns the sorted indices of variables that appear in the
// polynomial.
func (p Poly) UsedVars() []int {
	seen := make([]bool, p.nvars)
	for _, t := range p.terms {
		for i, e := range t.Mono {
			if e > 0 {
				seen[i] = true
			}
		}
	}
	var out []int
	for i, b := range seen {
		if b {
			out = append(out, i)
		}
	}
	return out
}

// DegreeIn returns the highest exponent of the i-th variable appearing in the
// polynomial, or 0 for the zero polynomial.
func (p Poly) DegreeIn(i int) int {
	d := 0
	for _, t := range p.terms {
		if t.Mono[i] > d {
			d = t.Mono[i]
		}
	}
	return d
}

// AddSlice returns the sum of a slice of polynomials in n variables.
func AddSlice(n int, ps []Poly) Poly {
	acc := Zero(n)
	for _, p := range ps {
		acc = acc.Add(p)
	}
	return acc
}

// EvalComplexAbs returns the magnitude |p(point)| of the complex evaluation,
// used to test whether a candidate solution satisfies the polynomial.
func (p Poly) EvalComplexAbs(point []complex128) float64 {
	return cmplx.Abs(p.EvalComplex(point))
}
