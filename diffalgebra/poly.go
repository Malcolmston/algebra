package diffalgebra

import (
	"math/big"
	"math/cmplx"
	"strings"
)

// Poly is a dense univariate polynomial over the rational numbers Q. The zero
// polynomial is represented by an empty coefficient slice. Coefficient i is the
// coefficient of x^i; the slice is always kept in normalised form with no
// trailing (highest-degree) zero coefficients.
type Poly struct {
	c []*big.Rat
}

// normalizePoly trims trailing zero coefficients and returns p.
func normalizePoly(c []*big.Rat) Poly {
	n := len(c)
	for n > 0 && ratZero(c[n-1]) {
		n--
	}
	return Poly{c: c[:n]}
}

// NewPoly builds a polynomial from coefficients given in ascending degree order
// (coeffs[0] is the constant term). The inputs are copied.
func NewPoly(coeffs ...*big.Rat) Poly {
	c := make([]*big.Rat, len(coeffs))
	for i, v := range coeffs {
		c[i] = cloneRat(v)
	}
	return normalizePoly(c)
}

// PolyFromInts builds a polynomial from integer coefficients in ascending
// degree order.
func PolyFromInts(coeffs ...int64) Poly {
	c := make([]*big.Rat, len(coeffs))
	for i, v := range coeffs {
		c[i] = ratInt(v)
	}
	return normalizePoly(c)
}

// ZeroPoly returns the zero polynomial.
func ZeroPoly() Poly { return Poly{} }

// ConstPoly returns the constant polynomial equal to r.
func ConstPoly(r *big.Rat) Poly { return normalizePoly([]*big.Rat{cloneRat(r)}) }

// ConstPolyInt returns the constant polynomial equal to the integer n.
func ConstPolyInt(n int64) Poly { return normalizePoly([]*big.Rat{ratInt(n)}) }

// OnePoly returns the constant polynomial 1.
func OnePoly() Poly { return ConstPolyInt(1) }

// XPoly returns the polynomial x.
func XPoly() Poly { return PolyFromInts(0, 1) }

// Monomial returns the polynomial coeff * x^deg.
func Monomial(coeff *big.Rat, deg int) Poly {
	if deg < 0 || ratZero(coeff) {
		return ZeroPoly()
	}
	c := make([]*big.Rat, deg+1)
	for i := 0; i < deg; i++ {
		c[i] = ratInt(0)
	}
	c[deg] = cloneRat(coeff)
	return normalizePoly(c)
}

// Degree returns the degree of p. The zero polynomial has degree -1.
func (p Poly) Degree() int { return len(p.c) - 1 }

// IsZero reports whether p is the zero polynomial.
func (p Poly) IsZero() bool { return len(p.c) == 0 }

// IsConstant reports whether p has degree at most zero.
func (p Poly) IsConstant() bool { return len(p.c) <= 1 }

// Coeff returns the coefficient of x^i, or zero when i is out of range.
func (p Poly) Coeff(i int) *big.Rat {
	if i < 0 || i >= len(p.c) {
		return ratInt(0)
	}
	return cloneRat(p.c[i])
}

// Coeffs returns a fresh slice of the coefficients in ascending degree order.
func (p Poly) Coeffs() []*big.Rat {
	out := make([]*big.Rat, len(p.c))
	for i, v := range p.c {
		out[i] = cloneRat(v)
	}
	return out
}

// LeadingCoeff returns the coefficient of the highest-degree term, or zero for
// the zero polynomial.
func (p Poly) LeadingCoeff() *big.Rat {
	if p.IsZero() {
		return ratInt(0)
	}
	return cloneRat(p.c[len(p.c)-1])
}

// ConstantTerm returns the coefficient of x^0.
func (p Poly) ConstantTerm() *big.Rat { return p.Coeff(0) }

// Clone returns a deep copy of p.
func (p Poly) Clone() Poly { return NewPoly(p.c...) }

// Equal reports whether p and q are equal as polynomials.
func (p Poly) Equal(q Poly) bool {
	if len(p.c) != len(q.c) {
		return false
	}
	for i := range p.c {
		if p.c[i].Cmp(q.c[i]) != 0 {
			return false
		}
	}
	return true
}

// Neg returns -p.
func (p Poly) Neg() Poly {
	c := make([]*big.Rat, len(p.c))
	for i, v := range p.c {
		c[i] = ratNeg(v)
	}
	return normalizePoly(c)
}

// Add returns p+q.
func (p Poly) Add(q Poly) Poly {
	n := len(p.c)
	if len(q.c) > n {
		n = len(q.c)
	}
	c := make([]*big.Rat, n)
	for i := 0; i < n; i++ {
		c[i] = ratAdd(p.Coeff(i), q.Coeff(i))
	}
	return normalizePoly(c)
}

// Sub returns p-q.
func (p Poly) Sub(q Poly) Poly { return p.Add(q.Neg()) }

// ScalarMul returns r*p for a rational scalar r.
func (p Poly) ScalarMul(r *big.Rat) Poly {
	if ratZero(r) {
		return ZeroPoly()
	}
	c := make([]*big.Rat, len(p.c))
	for i, v := range p.c {
		c[i] = ratMul(v, r)
	}
	return normalizePoly(c)
}

// Mul returns p*q.
func (p Poly) Mul(q Poly) Poly {
	if p.IsZero() || q.IsZero() {
		return ZeroPoly()
	}
	c := make([]*big.Rat, len(p.c)+len(q.c)-1)
	for i := range c {
		c[i] = ratInt(0)
	}
	for i, a := range p.c {
		if ratZero(a) {
			continue
		}
		for j, b := range q.c {
			c[i+j] = ratAdd(c[i+j], ratMul(a, b))
		}
	}
	return normalizePoly(c)
}

// Pow returns p raised to the non-negative integer power n.
func (p Poly) Pow(n int) Poly {
	result := OnePoly()
	base := p.Clone()
	for n > 0 {
		if n&1 == 1 {
			result = result.Mul(base)
		}
		base = base.Mul(base)
		n >>= 1
	}
	return result
}

// DivMod returns the quotient and remainder of Euclidean division p = q*d + r
// with deg r < deg d. It returns ErrDivByZero when d is zero.
func (p Poly) DivMod(d Poly) (q, r Poly, err error) {
	if d.IsZero() {
		return ZeroPoly(), ZeroPoly(), ErrDivByZero
	}
	rem := p.Clone()
	dDeg := d.Degree()
	dLead := d.LeadingCoeff()
	quo := make([]*big.Rat, 0)
	for !rem.IsZero() && rem.Degree() >= dDeg {
		shift := rem.Degree() - dDeg
		factor := ratDiv(rem.LeadingCoeff(), dLead)
		term := Monomial(factor, shift)
		rem = rem.Sub(term.Mul(d))
		// accumulate quotient coefficient
		for len(quo) <= shift {
			quo = append(quo, ratInt(0))
		}
		quo[shift] = ratAdd(quo[shift], factor)
	}
	return normalizePoly(quo), rem, nil
}

// Quo returns the quotient of p divided by d (see DivMod).
func (p Poly) Quo(d Poly) (Poly, error) {
	q, _, err := p.DivMod(d)
	return q, err
}

// Rem returns the remainder of p divided by d (see DivMod).
func (p Poly) Rem(d Poly) (Poly, error) {
	_, r, err := p.DivMod(d)
	return r, err
}

// Monic returns p scaled to be monic (leading coefficient 1). The zero
// polynomial is returned unchanged.
func (p Poly) Monic() Poly {
	if p.IsZero() {
		return p
	}
	return p.ScalarMul(ratInv(p.LeadingCoeff()))
}

// IsMonic reports whether p is monic.
func (p Poly) IsMonic() bool {
	return !p.IsZero() && p.LeadingCoeff().Cmp(ratInt(1)) == 0
}

// GCD returns the monic greatest common divisor of p and q. The GCD of two zero
// polynomials is zero.
func (p Poly) GCD(q Poly) Poly {
	a, b := p.Clone(), q.Clone()
	for !b.IsZero() {
		_, r, _ := a.DivMod(b)
		a, b = b, r
	}
	return a.Monic()
}

// ExtendedGCD returns the monic gcd g of p and q together with cofactors s and
// t satisfying s*p + t*q = g.
func (p Poly) ExtendedGCD(q Poly) (g, s, t Poly) {
	oldR, r := p.Clone(), q.Clone()
	oldS, sc := OnePoly(), ZeroPoly()
	oldT, tc := ZeroPoly(), OnePoly()
	for !r.IsZero() {
		quo, rem, _ := oldR.DivMod(r)
		oldR, r = r, rem
		oldS, sc = sc, oldS.Sub(quo.Mul(sc))
		oldT, tc = tc, oldT.Sub(quo.Mul(tc))
	}
	if oldR.IsZero() {
		return oldR, oldS, oldT
	}
	inv := ratInv(oldR.LeadingCoeff())
	return oldR.ScalarMul(inv), oldS.ScalarMul(inv), oldT.ScalarMul(inv)
}

// Derivative returns the formal derivative dp/dx.
func (p Poly) Derivative() Poly {
	if len(p.c) <= 1 {
		return ZeroPoly()
	}
	c := make([]*big.Rat, len(p.c)-1)
	for i := 1; i < len(p.c); i++ {
		c[i-1] = ratMul(p.c[i], ratInt(int64(i)))
	}
	return normalizePoly(c)
}

// Integral returns the antiderivative of p with zero constant of integration.
func (p Poly) Integral() Poly {
	if p.IsZero() {
		return ZeroPoly()
	}
	c := make([]*big.Rat, len(p.c)+1)
	c[0] = ratInt(0)
	for i := 0; i < len(p.c); i++ {
		c[i+1] = ratDiv(p.c[i], ratInt(int64(i+1)))
	}
	return normalizePoly(c)
}

// EvalRat evaluates p at the rational point x using Horner's rule.
func (p Poly) EvalRat(x *big.Rat) *big.Rat {
	acc := ratInt(0)
	for i := len(p.c) - 1; i >= 0; i-- {
		acc = ratAdd(ratMul(acc, x), p.c[i])
	}
	return acc
}

// EvalFloat evaluates p at the floating-point point x.
func (p Poly) EvalFloat(x float64) float64 {
	acc := 0.0
	for i := len(p.c) - 1; i >= 0; i-- {
		acc = acc*x + RatToFloat(p.c[i])
	}
	return acc
}

// EvalComplex evaluates p at the complex point z.
func (p Poly) EvalComplex(z complex128) complex128 {
	acc := complex(0, 0)
	for i := len(p.c) - 1; i >= 0; i-- {
		acc = acc*z + complex(RatToFloat(p.c[i]), 0)
	}
	return acc
}

// Compose returns p(q(x)).
func (p Poly) Compose(q Poly) Poly {
	acc := ZeroPoly()
	for i := len(p.c) - 1; i >= 0; i-- {
		acc = acc.Mul(q).Add(ConstPoly(p.c[i]))
	}
	return acc
}

// Shift returns p(x + a) for a rational shift a.
func (p Poly) Shift(a *big.Rat) Poly {
	return p.Compose(NewPoly(a, ratInt(1)))
}

// Reverse returns the reversal x^deg * p(1/x), i.e. p with its coefficient list
// reversed.
func (p Poly) Reverse() Poly {
	if p.IsZero() {
		return p
	}
	n := len(p.c)
	c := make([]*big.Rat, n)
	for i := 0; i < n; i++ {
		c[i] = cloneRat(p.c[n-1-i])
	}
	return normalizePoly(c)
}

// Content returns the (positive) rational content of p: the gcd of the integer
// numerators divided by the lcm of the denominators. It returns zero for the
// zero polynomial.
func (p Poly) Content() *big.Rat {
	if p.IsZero() {
		return ratInt(0)
	}
	// lcm of denominators
	lcm := big.NewInt(1)
	for _, v := range p.c {
		d := v.Denom()
		g := new(big.Int).GCD(nil, nil, lcm, d)
		lcm.Div(lcm, g)
		lcm.Mul(lcm, d)
	}
	// gcd of numerators after scaling by lcm
	num := big.NewInt(0)
	for _, v := range p.c {
		scaled := new(big.Int).Mul(v.Num(), new(big.Int).Div(lcm, v.Denom()))
		scaled.Abs(scaled)
		num.GCD(nil, nil, num, scaled)
	}
	return new(big.Rat).SetFrac(num, lcm)
}

// PrimitivePart returns p divided by its content, so the result has integer
// coprime coefficients and positive leading sign preserved.
func (p Poly) PrimitivePart() Poly {
	if p.IsZero() {
		return p
	}
	return p.ScalarMul(ratInv(p.Content()))
}

// SquareFreePart returns the square-free part of p, namely p / gcd(p, p').
func (p Poly) SquareFreePart() Poly {
	if p.IsZero() {
		return p
	}
	g := p.GCD(p.Derivative())
	q, _, _ := p.DivMod(g)
	return q.Monic()
}

// SquareFreeFactor is one factor of a square-free factorisation: the monic
// square-free polynomial Factor appears to the power Mult in the original.
type SquareFreeFactor struct {
	Factor Poly
	Mult   int
}

// SquareFreeFactorization returns the square-free factorisation of p using
// Yun's algorithm. The returned factors are monic, pairwise coprime and
// square-free, and their product raised to the listed multiplicities equals the
// monic part of p. Constant and zero inputs yield an empty slice.
func (p Poly) SquareFreeFactorization() []SquareFreeFactor {
	if p.Degree() < 1 {
		return nil
	}
	f := p.Monic()
	fp := f.Derivative()
	a0 := f.GCD(fp)
	b, _, _ := f.DivMod(a0)
	c, _, _ := fp.DivMod(a0)
	d := c.Sub(b.Derivative())
	var out []SquareFreeFactor
	i := 1
	for b.Degree() >= 1 {
		a := b.GCD(d)
		if a.Degree() >= 1 {
			out = append(out, SquareFreeFactor{Factor: a.Monic(), Mult: i})
		}
		b, _, _ = b.DivMod(a)
		c, _, _ = d.DivMod(a)
		d = c.Sub(b.Derivative())
		i++
	}
	return out
}

// Resultant returns the resultant Res(p, q) computed with the Euclidean
// recurrence. It is zero exactly when p and q share a non-constant factor.
func (p Poly) Resultant(q Poly) *big.Rat {
	return resultant(p.Clone(), q.Clone())
}

func resultant(a, b Poly) *big.Rat {
	if a.IsZero() || b.IsZero() {
		return ratInt(0)
	}
	if a.Degree() == 0 {
		return ratPow(a.LeadingCoeff(), b.Degree())
	}
	if b.Degree() == 0 {
		return ratPow(b.LeadingCoeff(), a.Degree())
	}
	_, r, _ := a.DivMod(b)
	if r.IsZero() {
		return ratInt(0)
	}
	da, db := a.Degree(), b.Degree()
	sign := ratInt(1)
	if (da*db)%2 == 1 {
		sign = ratInt(-1)
	}
	coeff := ratPow(b.LeadingCoeff(), da-r.Degree())
	return ratMul(ratMul(sign, coeff), resultant(b, r))
}

// Discriminant returns the discriminant of p, defined via the resultant of p
// and its derivative.
func (p Poly) Discriminant() *big.Rat {
	n := p.Degree()
	if n < 1 {
		return ratInt(0)
	}
	res := p.Resultant(p.Derivative())
	sign := ratInt(1)
	if (n*(n-1)/2)%2 == 1 {
		sign = ratInt(-1)
	}
	return ratDiv(ratMul(sign, res), p.LeadingCoeff())
}

// RationalRoots returns all distinct rational roots of p (without multiplicity)
// using the rational-root theorem. The zero polynomial yields nil.
func (p Poly) RationalRoots() []*big.Rat {
	if p.IsZero() || p.Degree() == 0 {
		return nil
	}
	q := p.PrimitivePart()
	// strip x factors (root 0)
	var roots []*big.Rat
	k := 0
	for k < len(q.c) && ratZero(q.c[k]) {
		k++
	}
	if k > 0 {
		roots = append(roots, ratInt(0))
		q = normalizePoly(q.c[k:])
	}
	if q.Degree() == 0 {
		return roots
	}
	a0 := new(big.Int).Abs(q.c[0].Num())
	an := new(big.Int).Abs(q.c[len(q.c)-1].Num())
	pDivs := intDivisors(a0)
	qDivs := intDivisors(an)
	seen := map[string]bool{}
	for _, pd := range pDivs {
		for _, qd := range qDivs {
			for _, sign := range []int64{1, -1} {
				cand := new(big.Rat).SetFrac(new(big.Int).Mul(pd, big.NewInt(sign)), qd)
				key := cand.RatString()
				if seen[key] {
					continue
				}
				seen[key] = true
				if ratZero(q.EvalRat(cand)) {
					roots = append(roots, cand)
				}
			}
		}
	}
	return roots
}

// intDivisors returns the positive divisors of |n| (with 1 for n == 0).
func intDivisors(n *big.Int) []*big.Int {
	if n.Sign() == 0 {
		return []*big.Int{big.NewInt(1)}
	}
	abs := new(big.Int).Abs(n)
	var divs []*big.Int
	one := big.NewInt(1)
	i := big.NewInt(1)
	for i.Cmp(abs) <= 0 {
		rem := new(big.Int).Mod(abs, i)
		if rem.Sign() == 0 {
			divs = append(divs, new(big.Int).Set(i))
		}
		i = new(big.Int).Add(i, one)
	}
	return divs
}

// String renders p in descending-degree human-readable form, e.g. "x^2 - 2".
func (p Poly) String() string {
	if p.IsZero() {
		return "0"
	}
	var b strings.Builder
	first := true
	for i := len(p.c) - 1; i >= 0; i-- {
		v := p.c[i]
		if ratZero(v) {
			continue
		}
		neg := v.Sign() < 0
		abs := v
		if neg {
			abs = ratNeg(v)
		}
		if first {
			if neg {
				b.WriteString("-")
			}
			first = false
		} else {
			if neg {
				b.WriteString(" - ")
			} else {
				b.WriteString(" + ")
			}
		}
		showCoeff := abs.Cmp(ratInt(1)) != 0 || i == 0
		if showCoeff {
			b.WriteString(abs.RatString())
		}
		if i > 0 {
			if showCoeff {
				b.WriteString("*")
			}
			b.WriteString("x")
			if i > 1 {
				b.WriteString("^")
				b.WriteString(itoa(i))
			}
		}
	}
	return b.String()
}

func itoa(n int) string {
	return new(big.Int).SetInt64(int64(n)).String()
}

// ComplexRootsFloat returns the roots of p as complex numbers using the
// Durand-Kerner iteration seeded by the given value. The zero polynomial and
// constants return nil.
func (p Poly) ComplexRootsFloat(seed int64) []complex128 {
	n := p.Degree()
	if n < 1 {
		return nil
	}
	coeffs := make([]complex128, n+1)
	for i := 0; i <= n; i++ {
		coeffs[i] = complex(RatToFloat(p.c[i]), 0)
	}
	return durandKerner(coeffs, seed)
}

// polyMaxAbs returns the maximum magnitude of the coefficients, used for
// scaling convergence tests.
func cabs(z complex128) float64 { return cmplx.Abs(z) }
