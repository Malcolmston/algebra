package galois

import (
	"errors"
	"math/big"
	"strings"
)

// Poly is a dense univariate polynomial over the prime field GF(p).
// Coefficients are stored little-endian: Coeff[i] is the coefficient of x^i,
// reduced into [0, p). A normalised polynomial has no trailing zero
// coefficients, so the zero polynomial has an empty Coeff slice.
type Poly struct {
	P     *big.Int
	Coeff []*big.Int
}

// normalize reduces every coefficient modulo P and trims trailing zeros.
func (a *Poly) normalize() *Poly {
	for i := range a.Coeff {
		a.Coeff[i] = new(big.Int).Mod(a.Coeff[i], a.P)
	}
	n := len(a.Coeff)
	for n > 0 && a.Coeff[n-1].Sign() == 0 {
		n--
	}
	a.Coeff = a.Coeff[:n]
	return a
}

// NewPoly builds a polynomial over GF(p) from int64 coefficients given
// little-endian (constant term first).
func NewPoly(p *big.Int, coeffs ...int64) *Poly {
	c := make([]*big.Int, len(coeffs))
	for i, v := range coeffs {
		c[i] = big.NewInt(v)
	}
	return (&Poly{P: clone(p), Coeff: c}).normalize()
}

// NewPolyBig builds a polynomial over GF(p) from big.Int coefficients given
// little-endian. The input slice is copied.
func NewPolyBig(p *big.Int, coeffs []*big.Int) *Poly {
	c := make([]*big.Int, len(coeffs))
	for i, v := range coeffs {
		c[i] = clone(v)
	}
	return (&Poly{P: clone(p), Coeff: c}).normalize()
}

// PolyZero returns the zero polynomial over GF(p).
func PolyZero(p *big.Int) *Poly {
	return &Poly{P: clone(p), Coeff: nil}
}

// PolyOne returns the constant polynomial 1 over GF(p).
func PolyOne(p *big.Int) *Poly {
	return NewPoly(p, 1)
}

// PolyConst returns the constant polynomial c over GF(p).
func PolyConst(p, c *big.Int) *Poly {
	return NewPolyBig(p, []*big.Int{c})
}

// PolyX returns the monomial x over GF(p).
func PolyX(p *big.Int) *Poly {
	return NewPoly(p, 0, 1)
}

// PolyMonomial returns the monomial coeff·x^deg over GF(p).
func PolyMonomial(p, coeff *big.Int, deg int) *Poly {
	if deg < 0 {
		deg = 0
	}
	c := make([]*big.Int, deg+1)
	for i := 0; i < deg; i++ {
		c[i] = big.NewInt(0)
	}
	c[deg] = clone(coeff)
	return (&Poly{P: clone(p), Coeff: c}).normalize()
}

// Degree returns the degree of the polynomial. The zero polynomial has
// degree -1 by convention.
func (a *Poly) Degree() int {
	return len(a.Coeff) - 1
}

// IsZero reports whether the polynomial is the zero polynomial.
func (a *Poly) IsZero() bool {
	return len(a.Coeff) == 0
}

// IsOne reports whether the polynomial is the constant 1.
func (a *Poly) IsOne() bool {
	return len(a.Coeff) == 1 && a.Coeff[0].Cmp(big1) == 0
}

// IsConstant reports whether the polynomial has degree 0 or is zero.
func (a *Poly) IsConstant() bool {
	return len(a.Coeff) <= 1
}

// LeadingCoeff returns the coefficient of the highest-degree term, or 0 for the
// zero polynomial.
func (a *Poly) LeadingCoeff() *big.Int {
	if a.IsZero() {
		return big.NewInt(0)
	}
	return clone(a.Coeff[len(a.Coeff)-1])
}

// Coefficient returns the coefficient of x^i, or 0 when i is out of range.
func (a *Poly) Coefficient(i int) *big.Int {
	if i < 0 || i >= len(a.Coeff) {
		return big.NewInt(0)
	}
	return clone(a.Coeff[i])
}

// Copy returns an independent deep copy of the polynomial.
func (a *Poly) Copy() *Poly {
	return NewPolyBig(a.P, a.Coeff)
}

// Equal reports whether a and b are equal as polynomials over the same field.
func (a *Poly) Equal(b *Poly) bool {
	if a.P.Cmp(b.P) != 0 || len(a.Coeff) != len(b.Coeff) {
		return false
	}
	for i := range a.Coeff {
		if a.Coeff[i].Cmp(b.Coeff[i]) != 0 {
			return false
		}
	}
	return true
}

// PolyEqual reports whether two polynomials are equal.
func PolyEqual(a, b *Poly) bool { return a.Equal(b) }

// Add returns the sum a + b over GF(p).
func (a *Poly) Add(b *Poly) *Poly {
	n := len(a.Coeff)
	if len(b.Coeff) > n {
		n = len(b.Coeff)
	}
	c := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		c[i] = AddMod(a.Coefficient(i), b.Coefficient(i), a.P)
	}
	return (&Poly{P: clone(a.P), Coeff: c}).normalize()
}

// Sub returns the difference a - b over GF(p).
func (a *Poly) Sub(b *Poly) *Poly {
	n := len(a.Coeff)
	if len(b.Coeff) > n {
		n = len(b.Coeff)
	}
	c := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		c[i] = SubMod(a.Coefficient(i), b.Coefficient(i), a.P)
	}
	return (&Poly{P: clone(a.P), Coeff: c}).normalize()
}

// Neg returns the additive inverse -a over GF(p).
func (a *Poly) Neg() *Poly {
	c := make([]*big.Int, len(a.Coeff))
	for i := range a.Coeff {
		c[i] = NegMod(a.Coeff[i], a.P)
	}
	return (&Poly{P: clone(a.P), Coeff: c}).normalize()
}

// ScalarMul returns the polynomial scaled by the field element s.
func (a *Poly) ScalarMul(s *big.Int) *Poly {
	c := make([]*big.Int, len(a.Coeff))
	for i := range a.Coeff {
		c[i] = MulMod(a.Coeff[i], s, a.P)
	}
	return (&Poly{P: clone(a.P), Coeff: c}).normalize()
}

// Mul returns the product a * b over GF(p).
func (a *Poly) Mul(b *Poly) *Poly {
	if a.IsZero() || b.IsZero() {
		return PolyZero(a.P)
	}
	c := make([]*big.Int, len(a.Coeff)+len(b.Coeff)-1)
	for i := range c {
		c[i] = big.NewInt(0)
	}
	for i, ai := range a.Coeff {
		if ai.Sign() == 0 {
			continue
		}
		for j, bj := range b.Coeff {
			t := new(big.Int).Mul(ai, bj)
			c[i+j].Add(c[i+j], t)
		}
	}
	return (&Poly{P: clone(a.P), Coeff: c}).normalize()
}

// DivMod returns the quotient q and remainder r of a divided by b, satisfying
// a = q*b + r with deg(r) < deg(b). It returns an error when b is zero.
func (a *Poly) DivMod(b *Poly) (q, r *Poly, err error) {
	if b.IsZero() {
		return nil, nil, errors.New("galois: polynomial division by zero")
	}
	binv, err := InvMod(b.LeadingCoeff(), a.P)
	if err != nil {
		return nil, nil, err
	}
	r = a.Copy()
	q = PolyZero(a.P)
	db := b.Degree()
	for !r.IsZero() && r.Degree() >= db {
		shift := r.Degree() - db
		coef := MulMod(r.LeadingCoeff(), binv, a.P)
		term := PolyMonomial(a.P, coef, shift)
		q = q.Add(term)
		r = r.Sub(term.Mul(b))
	}
	return q, r, nil
}

// Quo returns the quotient of a divided by b.
func (a *Poly) Quo(b *Poly) (*Poly, error) {
	q, _, err := a.DivMod(b)
	return q, err
}

// Rem returns the remainder of a divided by b (a mod b).
func (a *Poly) Rem(b *Poly) (*Poly, error) {
	_, r, err := a.DivMod(b)
	return r, err
}

// Pow returns a raised to the non-negative integer power e.
func (a *Poly) Pow(e int) *Poly {
	result := PolyOne(a.P)
	base := a.Copy()
	for e > 0 {
		if e&1 == 1 {
			result = result.Mul(base)
		}
		base = base.Mul(base)
		e >>= 1
	}
	return result
}

// PowMod returns a^e reduced modulo the polynomial m, using square-and-multiply.
// The exponent e must be non-negative.
func (a *Poly) PowMod(e *big.Int, m *Poly) (*Poly, error) {
	if e.Sign() < 0 {
		return nil, errors.New("galois: PowMod requires a non-negative exponent")
	}
	result := PolyOne(a.P)
	if _, result, _ = result.DivMod(m); result == nil {
		result = PolyOne(a.P)
	}
	base := a.Copy()
	if _, base, _ = base.DivMod(m); base == nil {
		base = PolyZero(a.P)
	}
	ee := clone(e)
	for ee.Sign() > 0 {
		if ee.Bit(0) == 1 {
			result = result.Mul(base)
			if _, rr, err := result.DivMod(m); err == nil {
				result = rr
			}
		}
		ee.Rsh(ee, 1)
		if ee.Sign() > 0 {
			base = base.Mul(base)
			if _, bb, err := base.DivMod(m); err == nil {
				base = bb
			}
		}
	}
	return result, nil
}

// Gcd returns the monic greatest common divisor of a and b over GF(p). The gcd
// of two zero polynomials is zero.
func (a *Poly) Gcd(b *Poly) *Poly {
	x := a.Copy()
	y := b.Copy()
	for !y.IsZero() {
		_, r, err := x.DivMod(y)
		if err != nil {
			break
		}
		x, y = y, r
	}
	if x.IsZero() {
		return x
	}
	return x.Monic()
}

// PolyGcd returns the monic greatest common divisor of a and b.
func PolyGcd(a, b *Poly) *Poly { return a.Gcd(b) }

// PolyLcm returns the monic least common multiple of a and b.
func PolyLcm(a, b *Poly) *Poly {
	if a.IsZero() || b.IsZero() {
		return PolyZero(a.P)
	}
	g := a.Gcd(b)
	q, _ := a.Mul(b).Quo(g)
	return q.Monic()
}

// ExtendedGcd returns g = gcd(a, b) together with cofactors s and t satisfying
// a*s + b*t = g. The returned g is monic (unless it is zero).
func (a *Poly) ExtendedGcd(b *Poly) (g, s, t *Poly) {
	oldR, r := a.Copy(), b.Copy()
	oldS, sCur := PolyOne(a.P), PolyZero(a.P)
	oldT, tCur := PolyZero(a.P), PolyOne(a.P)
	for !r.IsZero() {
		q, rem, err := oldR.DivMod(r)
		if err != nil {
			break
		}
		oldR, r = r, rem
		oldS, sCur = sCur, oldS.Sub(q.Mul(sCur))
		oldT, tCur = tCur, oldT.Sub(q.Mul(tCur))
	}
	if oldR.IsZero() {
		return oldR, oldS, oldT
	}
	inv, _ := InvMod(oldR.LeadingCoeff(), a.P)
	return oldR.ScalarMul(inv), oldS.ScalarMul(inv), oldT.ScalarMul(inv)
}

// Derivative returns the formal derivative of the polynomial over GF(p).
func (a *Poly) Derivative() *Poly {
	if a.Degree() < 1 {
		return PolyZero(a.P)
	}
	c := make([]*big.Int, len(a.Coeff)-1)
	for i := 1; i < len(a.Coeff); i++ {
		c[i-1] = MulMod(a.Coeff[i], big.NewInt(int64(i)), a.P)
	}
	return (&Poly{P: clone(a.P), Coeff: c}).normalize()
}

// Eval evaluates the polynomial at the field element x using Horner's rule,
// returning a value in [0, p).
func (a *Poly) Eval(x *big.Int) *big.Int {
	acc := big.NewInt(0)
	xm := new(big.Int).Mod(x, a.P)
	for i := len(a.Coeff) - 1; i >= 0; i-- {
		acc.Mul(acc, xm)
		acc.Add(acc, a.Coeff[i])
		acc.Mod(acc, a.P)
	}
	return acc
}

// Compose returns the composition a(b(x)) over GF(p).
func (a *Poly) Compose(b *Poly) *Poly {
	result := PolyZero(a.P)
	for i := len(a.Coeff) - 1; i >= 0; i-- {
		result = result.Mul(b).Add(PolyConst(a.P, a.Coeff[i]))
	}
	return result
}

// Monic returns a scaled so that its leading coefficient is 1. The zero
// polynomial is returned unchanged.
func (a *Poly) Monic() *Poly {
	if a.IsZero() {
		return a.Copy()
	}
	inv, err := InvMod(a.LeadingCoeff(), a.P)
	if err != nil {
		return a.Copy()
	}
	return a.ScalarMul(inv)
}

// IsMonic reports whether the leading coefficient is 1.
func (a *Poly) IsMonic() bool {
	return !a.IsZero() && a.LeadingCoeff().Cmp(big1) == 0
}

// PolyMulMod returns (a * b) mod m over GF(p).
func PolyMulMod(a, b, m *Poly) *Poly {
	_, r, err := a.Mul(b).DivMod(m)
	if err != nil {
		return PolyZero(a.P)
	}
	return r
}

// PolyPowMod returns a^e mod m over GF(p) for a non-negative exponent e.
func PolyPowMod(a *Poly, e *big.Int, m *Poly) (*Poly, error) {
	return a.PowMod(e, m)
}

// String renders the polynomial in descending-degree infix form such as
// "x^2 + 3x + 1", or "0" for the zero polynomial.
func (a *Poly) String() string {
	if a.IsZero() {
		return "0"
	}
	var parts []string
	for i := len(a.Coeff) - 1; i >= 0; i-- {
		c := a.Coeff[i]
		if c.Sign() == 0 {
			continue
		}
		var term string
		switch i {
		case 0:
			term = c.String()
		case 1:
			if c.Cmp(big1) == 0 {
				term = "x"
			} else {
				term = c.String() + "x"
			}
		default:
			if c.Cmp(big1) == 0 {
				term = "x^" + itoa(i)
			} else {
				term = c.String() + "x^" + itoa(i)
			}
		}
		parts = append(parts, term)
	}
	return strings.Join(parts, " + ")
}

func itoa(n int) string {
	return big.NewInt(int64(n)).String()
}
