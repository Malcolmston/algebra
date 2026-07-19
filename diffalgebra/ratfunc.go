package diffalgebra

import (
	"math/big"
	"strings"
)

// RatFunc is a rational function p/q over Q, kept in reduced form with a monic
// denominator. The zero rational function has numerator zero and denominator
// one.
type RatFunc struct {
	num Poly
	den Poly
}

// NewRatFunc builds the reduced rational function num/den. It returns
// ErrDivByZero when den is the zero polynomial.
func NewRatFunc(num, den Poly) (RatFunc, error) {
	if den.IsZero() {
		return RatFunc{}, ErrDivByZero
	}
	return reduceRatFunc(num, den), nil
}

// mustRat is an internal constructor used when the denominator is known to be
// non-zero.
func mustRat(num, den Poly) RatFunc {
	return reduceRatFunc(num, den)
}

// reduceRatFunc cancels the gcd and normalises the denominator to be monic.
func reduceRatFunc(num, den Poly) RatFunc {
	if num.IsZero() {
		return RatFunc{num: ZeroPoly(), den: OnePoly()}
	}
	g := num.GCD(den)
	n, _, _ := num.DivMod(g)
	d, _, _ := den.DivMod(g)
	// make denominator monic
	lc := d.LeadingCoeff()
	if lc.Cmp(ratInt(1)) != 0 {
		inv := ratInv(lc)
		n = n.ScalarMul(inv)
		d = d.ScalarMul(inv)
	}
	return RatFunc{num: n, den: d}
}

// RatFuncFromPoly returns the rational function p/1.
func RatFuncFromPoly(p Poly) RatFunc { return RatFunc{num: p.Clone(), den: OnePoly()} }

// ConstRatFunc returns the constant rational function equal to r.
func ConstRatFunc(r *big.Rat) RatFunc { return RatFuncFromPoly(ConstPoly(r)) }

// ZeroRatFunc returns the zero rational function.
func ZeroRatFunc() RatFunc { return RatFunc{num: ZeroPoly(), den: OnePoly()} }

// OneRatFunc returns the constant rational function 1.
func OneRatFunc() RatFunc { return RatFuncFromPoly(OnePoly()) }

// XRatFunc returns the rational function x.
func XRatFunc() RatFunc { return RatFuncFromPoly(XPoly()) }

// Num returns a copy of the numerator.
func (f RatFunc) Num() Poly { return f.num.Clone() }

// Den returns a copy of the denominator.
func (f RatFunc) Den() Poly { return f.den.Clone() }

// IsZero reports whether f is the zero rational function.
func (f RatFunc) IsZero() bool { return f.num.IsZero() }

// IsPolynomial reports whether the denominator is a constant.
func (f RatFunc) IsPolynomial() bool { return f.den.IsConstant() }

// IsProper reports whether deg(num) < deg(den).
func (f RatFunc) IsProper() bool { return f.num.Degree() < f.den.Degree() }

// Equal reports whether f and g are equal as rational functions.
func (f RatFunc) Equal(g RatFunc) bool {
	return f.num.Mul(g.den).Equal(g.num.Mul(f.den))
}

// Neg returns -f.
func (f RatFunc) Neg() RatFunc { return RatFunc{num: f.num.Neg(), den: f.den.Clone()} }

// Add returns f+g.
func (f RatFunc) Add(g RatFunc) RatFunc {
	num := f.num.Mul(g.den).Add(g.num.Mul(f.den))
	den := f.den.Mul(g.den)
	return reduceRatFunc(num, den)
}

// Sub returns f-g.
func (f RatFunc) Sub(g RatFunc) RatFunc { return f.Add(g.Neg()) }

// Mul returns f*g.
func (f RatFunc) Mul(g RatFunc) RatFunc {
	return reduceRatFunc(f.num.Mul(g.num), f.den.Mul(g.den))
}

// ScalarMul returns r*f for a rational scalar r.
func (f RatFunc) ScalarMul(r *big.Rat) RatFunc {
	return reduceRatFunc(f.num.ScalarMul(r), f.den.Clone())
}

// Inv returns 1/f. It returns ErrDivByZero when f is zero.
func (f RatFunc) Inv() (RatFunc, error) {
	if f.IsZero() {
		return RatFunc{}, ErrDivByZero
	}
	return reduceRatFunc(f.den, f.num), nil
}

// Div returns f/g. It returns ErrDivByZero when g is zero.
func (f RatFunc) Div(g RatFunc) (RatFunc, error) {
	if g.IsZero() {
		return RatFunc{}, ErrDivByZero
	}
	return reduceRatFunc(f.num.Mul(g.den), f.den.Mul(g.num)), nil
}

// Pow raises f to the integer power n (n may be negative). A negative power of
// the zero function returns ErrDivByZero.
func (f RatFunc) Pow(n int) (RatFunc, error) {
	if n < 0 {
		inv, err := f.Inv()
		if err != nil {
			return RatFunc{}, err
		}
		return inv.Pow(-n)
	}
	result := OneRatFunc()
	base := f
	for n > 0 {
		if n&1 == 1 {
			result = result.Mul(base)
		}
		base = base.Mul(base)
		n >>= 1
	}
	return result, nil
}

// Derivative returns the derivative df/dx via the quotient rule.
func (f RatFunc) Derivative() RatFunc {
	// (n/d)' = (n' d - n d') / d^2
	num := f.num.Derivative().Mul(f.den).Sub(f.num.Mul(f.den.Derivative()))
	den := f.den.Mul(f.den)
	return reduceRatFunc(num, den)
}

// LogDerivative returns the logarithmic derivative f'/f. It returns
// ErrDivByZero when f is zero.
func (f RatFunc) LogDerivative() (RatFunc, error) {
	if f.IsZero() {
		return RatFunc{}, ErrDivByZero
	}
	return f.Derivative().Div(f)
}

// EvalRat evaluates f at the rational point x. It returns ErrDivByZero when the
// denominator vanishes there.
func (f RatFunc) EvalRat(x *big.Rat) (*big.Rat, error) {
	d := f.den.EvalRat(x)
	if ratZero(d) {
		return nil, ErrDivByZero
	}
	return ratDiv(f.num.EvalRat(x), d), nil
}

// EvalFloat evaluates f at the floating-point point x.
func (f RatFunc) EvalFloat(x float64) float64 {
	return f.num.EvalFloat(x) / f.den.EvalFloat(x)
}

// EvalComplex evaluates f at the complex point z.
func (f RatFunc) EvalComplex(z complex128) complex128 {
	return f.num.EvalComplex(z) / f.den.EvalComplex(z)
}

// PolynomialPart returns the polynomial quotient and proper remainder so that
// f = quotient + remainder with the remainder proper.
func (f RatFunc) PolynomialPart() (quotient Poly, remainder RatFunc) {
	q, r, _ := f.num.DivMod(f.den)
	return q, mustRat(r, f.den)
}

// String renders f as "(num)/(den)" or just the numerator when polynomial.
func (f RatFunc) String() string {
	if f.IsPolynomial() {
		return f.num.ScalarMul(ratInv(f.den.LeadingCoeff())).String()
	}
	var b strings.Builder
	b.WriteString("(")
	b.WriteString(f.num.String())
	b.WriteString(")/(")
	b.WriteString(f.den.String())
	b.WriteString(")")
	return b.String()
}

// PartialFractionTerm is a summand A/(F^k) of a partial-fraction decomposition,
// where F is a monic square-free polynomial appearing to power Power and
// deg(Numerator) < deg(F).
type PartialFractionTerm struct {
	Numerator Poly
	Factor    Poly
	Power     int
}

// PartialFractions returns the polynomial part together with the
// partial-fraction terms of f over Q, grouping the denominator by its
// square-free factorisation and reducing each power. The decomposition is exact
// over Q (residues that would be irrational are left folded inside the
// numerator over the corresponding square-free factor).
func (f RatFunc) PartialFractions() (Poly, []PartialFractionTerm) {
	poly, proper := f.PolynomialPart()
	if proper.IsZero() {
		return poly, nil
	}
	num := proper.num
	sqf := proper.den.SquareFreeFactorization()
	if len(sqf) == 0 {
		return poly, nil
	}
	// Build the coprime moduli M_i = F_i^i and CRT split the numerator.
	type modpart struct {
		factor Poly
		power  int
		mod    Poly // factor^power
		numer  Poly
	}
	var parts []modpart
	for _, sf := range sqf {
		parts = append(parts, modpart{factor: sf.Factor, power: sf.Mult, mod: sf.Factor.Pow(sf.Mult)})
	}
	D := proper.den
	for i := range parts {
		Ni, _, _ := D.DivMod(parts[i].mod) // Ni = D / mod, coprime to mod
		inv := invMod(Ni, parts[i].mod)    // Ni^{-1} mod (mod)
		ai := num.Mul(inv)
		_, r, _ := ai.DivMod(parts[i].mod)
		parts[i].numer = r
	}
	// Now expand each A_i/F_i^power into descending powers of F_i.
	var terms []PartialFractionTerm
	for _, mp := range parts {
		terms = append(terms, expandPrimePower(mp.numer, mp.factor, mp.power)...)
	}
	return poly, terms
}

// invMod returns the inverse of a modulo m (assuming gcd(a,m)=1).
func invMod(a, m Poly) Poly {
	g, s, _ := a.ExtendedGCD(m)
	// g should be constant; s*a + t*m = g, so a^{-1} = s/g mod m
	inv := s.ScalarMul(ratInv(g.LeadingCoeff()))
	_, r, _ := inv.DivMod(m)
	return r
}

// expandPrimePower writes a/f^power as a sum of terms b_j/f^j with deg b_j <
// deg f, by repeated division by f.
func expandPrimePower(a, f Poly, power int) []PartialFractionTerm {
	var terms []PartialFractionTerm
	cur := a
	for j := power; j >= 1; j-- {
		q, r, _ := cur.DivMod(f)
		if !r.IsZero() {
			terms = append(terms, PartialFractionTerm{Numerator: r, Factor: f.Clone(), Power: j})
		}
		cur = q
		if cur.IsZero() {
			break
		}
	}
	return terms
}
