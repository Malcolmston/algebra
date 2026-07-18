package algebra

import (
	"errors"
	"math/big"
	"sort"
	"strconv"
	"strings"
)

// This file adds a first-class dense univariate polynomial type over exact
// rational coefficients, together with polynomial arithmetic, GCD/LCM,
// square-free and rational factorization, resultant/discriminant invariants and
// partial-fraction decomposition. All results are exact where the coefficients
// are rational; non-rational coefficient inputs are rejected rather than being
// silently evaluated numerically. There is no global state; every routine is
// pure and deterministic.
//
// Internally a polynomial is stored as a slice of *big.Rat with index i holding
// the coefficient of Var^i, trimmed so the highest stored entry is non-zero. The
// zero polynomial is the empty slice. Unexported helpers are prefixed with "pa"
// (poly-algebra) to avoid collisions with the rest of the package.

// Poly is a dense univariate polynomial in a single named variable over the
// exact rationals. The coefficient of Var^i is stored at index i, using
// arbitrary-precision math/big rationals internally, so all arithmetic is exact.
// The zero polynomial has degree -1. Poly values are immutable once returned:
// every method builds and returns a fresh Poly rather than mutating the receiver.
type Poly struct {
	v      string
	coeffs []*big.Rat // index i = coefficient of Var^i, trimmed (highest entry non-zero); empty means zero
}

// PolyFactor is one factor of a polynomial together with its multiplicity, as
// returned by [Poly.SquareFree] and [Poly.Factor]. The factored polynomial
// equals the product of Base raised to Mult over all returned factors.
type PolyFactor struct {
	// Base is an irreducible or square-free factor (usually monic).
	Base *Poly
	// Mult is the multiplicity with which Base divides the original polynomial.
	Mult int
}

// PFTerm is a single term Coeff/(Denom^Power) of a partial-fraction
// decomposition, as returned by [PartialFractions]. Coeff is an expression whose
// degree in the variable is strictly less than the degree of Denom, and Denom is
// an irreducible (linear or quadratic) monic polynomial.
type PFTerm struct {
	// Coeff is the numerator expression of this term.
	Coeff Expr
	// Denom is the irreducible denominator polynomial.
	Denom *Poly
	// Power is the exponent applied to Denom in this term.
	Power int
}

// --- construction and conversion -------------------------------------------

// NewPoly returns the polynomial in varName whose coefficient of x^i is
// coeffs[i]. Each coefficient must be a numeric integer or rational expression;
// NewPoly panics if any coefficient is not exactly rational. Trailing zero
// coefficients are trimmed automatically, so NewPoly("x", Int(0)) is the zero
// polynomial.
func NewPoly(varName string, coeffs ...Expr) *Poly {
	rc := make([]*big.Rat, len(coeffs))
	for i, c := range coeffs {
		r, ok := paToRat(c)
		if !ok {
			panic("algebra: NewPoly requires integer or rational coefficients")
		}
		rc[i] = r
	}
	return &Poly{v: varName, coeffs: paTrim(rc)}
}

// PolyFrom extracts the polynomial in the symbol v from an arbitrary expression
// e. It returns an error if e is not polynomial in v (for example if v appears
// inside a function or with a negative or non-integer power) or if any
// coefficient is not exactly rational. v must be a [Symbol].
func PolyFrom(e Expr, v Expr) (*Poly, error) {
	s, ok := v.(*Symbol)
	if !ok {
		return nil, errors.New("algebra: PolyFrom requires a symbol")
	}
	cs, ok := polyCoeffs(e, s.Name)
	if !ok {
		return nil, errors.New("algebra: expression is not polynomial in " + s.Name)
	}
	rc, ok := ratCoeffs(cs)
	if !ok {
		return nil, errors.New("algebra: polynomial has non-rational coefficients")
	}
	return &Poly{v: s.Name, coeffs: paTrim(rc)}, nil
}

// Expr rebuilds an [Expr] equal to the polynomial, as a canonical sum of
// coefficient-times-power terms assembled through [Add], [Mul] and [Pow]. The
// zero polynomial yields Int(0).
func (p *Poly) Expr() Expr {
	if p.IsZero() {
		return Int(0)
	}
	terms := make([]Expr, 0, len(p.coeffs))
	x := Sym(p.v)
	for i, c := range p.coeffs {
		if c.Sign() == 0 {
			continue
		}
		terms = append(terms, Mul(newRational(new(big.Rat).Set(c)), Pow(x, Int(int64(i)))))
	}
	return Add(terms...)
}

// String renders the polynomial as a readable expression such as
// "x^2 - 2*x + 1", with terms ordered from the highest power down and the zero
// polynomial rendered as "0".
func (p *Poly) String() string {
	if p.IsZero() {
		return "0"
	}
	var b strings.Builder
	first := true
	for i := len(p.coeffs) - 1; i >= 0; i-- {
		c := p.coeffs[i]
		if c.Sign() == 0 {
			continue
		}
		neg := c.Sign() < 0
		mag := new(big.Rat).Abs(c)
		if first {
			if neg {
				b.WriteString("-")
			}
			first = false
		} else if neg {
			b.WriteString(" - ")
		} else {
			b.WriteString(" + ")
		}
		magStr := mag.RatString()
		isOne := mag.Cmp(big.NewRat(1, 1)) == 0
		switch {
		case i == 0:
			b.WriteString(magStr)
		case i == 1:
			if isOne {
				b.WriteString(p.v)
			} else {
				b.WriteString(magStr + "*" + p.v)
			}
		default:
			vp := p.v + "^" + strconv.Itoa(i)
			if isOne {
				b.WriteString(vp)
			} else {
				b.WriteString(magStr + "*" + vp)
			}
		}
	}
	return b.String()
}

// Var returns the name of the polynomial's variable.
func (p *Poly) Var() string { return p.v }

// Degree returns the degree of the polynomial, or -1 for the zero polynomial.
func (p *Poly) Degree() int { return len(p.coeffs) - 1 }

// Coeff returns the coefficient of x^i as an exact integer or rational
// expression. Indices outside the stored range (including negative indices)
// yield Int(0).
func (p *Poly) Coeff(i int) Expr {
	if i < 0 || i >= len(p.coeffs) {
		return Int(0)
	}
	return newRational(new(big.Rat).Set(p.coeffs[i]))
}

// LeadingCoeff returns the coefficient of the highest power as an exact integer
// or rational expression, or Int(0) for the zero polynomial.
func (p *Poly) LeadingCoeff() Expr {
	if p.IsZero() {
		return Int(0)
	}
	return p.Coeff(p.Degree())
}

// IsZero reports whether the polynomial is identically zero.
func (p *Poly) IsZero() bool { return len(p.coeffs) == 0 }

// Eval returns the value of the polynomial at x, computed by Horner's method and
// returned as a canonical [Expr]. For a numeric x the result folds to a number.
func (p *Poly) Eval(x Expr) Expr {
	if p.IsZero() {
		return Int(0)
	}
	acc := Expr(Int(0))
	for i := len(p.coeffs) - 1; i >= 0; i-- {
		acc = Add(Mul(acc, x), newRational(new(big.Rat).Set(p.coeffs[i])))
	}
	return Simplify(acc)
}

// --- arithmetic ------------------------------------------------------------

// Add returns the sum p + q. The result variable is taken from whichever operand
// is non-constant; Add panics if both operands are non-constant in different
// variables.
func (p *Poly) Add(q *Poly) *Poly {
	v := paBinVar(p, q)
	return &Poly{v: v, coeffs: paTrim(paAddC(p.coeffs, q.coeffs))}
}

// Sub returns the difference p - q. The result variable is chosen as in
// [Poly.Add].
func (p *Poly) Sub(q *Poly) *Poly {
	v := paBinVar(p, q)
	return &Poly{v: v, coeffs: paTrim(paSubC(p.coeffs, q.coeffs))}
}

// MulP returns the product p * q. The result variable is chosen as in
// [Poly.Add].
func (p *Poly) MulP(q *Poly) *Poly {
	v := paBinVar(p, q)
	return &Poly{v: v, coeffs: paTrim(paMulC(p.coeffs, q.coeffs))}
}

// Scale returns the polynomial with every coefficient multiplied by the scalar
// c, which must be an exact integer or rational; Scale panics otherwise.
func (p *Poly) Scale(c Expr) *Poly {
	r, ok := paToRat(c)
	if !ok {
		panic("algebra: Poly.Scale requires a rational scalar")
	}
	out := make([]*big.Rat, len(p.coeffs))
	for i, x := range p.coeffs {
		out[i] = new(big.Rat).Mul(x, r)
	}
	return &Poly{v: p.v, coeffs: paTrim(out)}
}

// DivMod returns the quotient and remainder of Euclidean division of p by q over
// the rationals, so that p == quo*q + rem with deg(rem) < deg(q). It returns an
// error if q is the zero polynomial.
func (p *Poly) DivMod(q *Poly) (quo, rem *Poly, err error) {
	if q.IsZero() {
		return nil, nil, errors.New("algebra: Poly.DivMod by zero polynomial")
	}
	v := paBinVar(p, q)
	qc, rc := paDivRat(p.coeffs, q.coeffs)
	return &Poly{v: v, coeffs: paTrim(qc)}, &Poly{v: v, coeffs: paTrim(rc)}, nil
}

// Derivative returns the formal derivative of the polynomial.
func (p *Poly) Derivative() *Poly {
	return &Poly{v: p.v, coeffs: paTrim(paDerivC(p.coeffs))}
}

// Monic returns the polynomial divided by its leading coefficient, yielding a
// monic polynomial with the same roots. The zero polynomial is returned
// unchanged.
func (p *Poly) Monic() *Poly {
	if p.IsZero() {
		return &Poly{v: p.v}
	}
	return &Poly{v: p.v, coeffs: paMonicC(p.coeffs)}
}

// --- gcd / lcm -------------------------------------------------------------

// PolyGCD returns the monic greatest common divisor of a and b, computed by the
// Euclidean algorithm over the rationals. The GCD of two zero polynomials is the
// zero polynomial.
func PolyGCD(a, b *Poly) *Poly {
	v := paBinVar(a, b)
	return &Poly{v: v, coeffs: paGCDC(a.coeffs, b.coeffs)}
}

// PolyLCM returns the monic least common multiple of a and b, computed as
// (a*b)/gcd(a,b). If either argument is the zero polynomial the result is the
// zero polynomial.
func PolyLCM(a, b *Poly) *Poly {
	v := paBinVar(a, b)
	if a.IsZero() || b.IsZero() {
		return &Poly{v: v}
	}
	g := PolyGCD(a, b)
	prod := a.MulP(b)
	quo, _, _ := prod.DivMod(g)
	return quo.Monic()
}

// --- factorization ---------------------------------------------------------

// SquareFree returns the square-free decomposition of the polynomial using Yun's
// algorithm: a list of pairwise-coprime square-free factors, each tagged with
// the multiplicity to which it divides the polynomial, so that the product of
// Base^Mult reproduces the polynomial. A non-unit leading coefficient is
// returned as a leading constant factor of multiplicity one. The zero polynomial
// yields no factors.
func (p *Poly) SquareFree() []PolyFactor {
	if p.IsZero() {
		return nil
	}
	if p.Degree() == 0 {
		return []PolyFactor{{Base: paClone(p), Mult: 1}}
	}
	var out []PolyFactor
	if lc := p.LeadingCoeff(); !isOne(lc) {
		out = append(out, PolyFactor{Base: NewPoly(p.v, lc), Mult: 1})
	}
	m := p.Monic()
	// Yun's algorithm over the monic polynomial m.
	a := PolyGCD(m, m.Derivative()) // gcd(f, f')
	b := paMustQuo(m, a)            // f / a
	c := paMustQuo(m.Derivative(), a)
	d := c.Sub(b.Derivative())
	for i := 1; b.Degree() >= 1; i++ {
		ai := PolyGCD(b, d)
		b = paMustQuo(b, ai)
		c = paMustQuo(d, ai)
		d = c.Sub(b.Derivative())
		if ai.Degree() >= 1 {
			out = append(out, PolyFactor{Base: ai, Mult: i})
		}
	}
	paSortFactors(out)
	return out
}

// Factor returns a factorization of the polynomial over the rationals: the
// product of Base^Mult reproduces the polynomial. Factoring proceeds by
// square-free decomposition followed by exact extraction of all rational roots
// (giving monic linear factors) and identification of irreducible quadratic
// factors; any remaining higher-degree factor with no rational root is returned
// as-is. A non-unit leading coefficient is returned as a leading constant
// factor. The result is deterministic.
func (p *Poly) Factor() []PolyFactor {
	if p.IsZero() {
		return nil
	}
	if p.Degree() == 0 {
		return []PolyFactor{{Base: paClone(p), Mult: 1}}
	}
	var out []PolyFactor
	if lc := p.LeadingCoeff(); !isOne(lc) {
		out = append(out, PolyFactor{Base: NewPoly(p.v, lc), Mult: 1})
	}
	for _, sf := range p.Monic().SquareFree() {
		if sf.Base.Degree() < 1 {
			continue
		}
		for _, irr := range paFactorSquareFree(sf.Base) {
			out = append(out, PolyFactor{Base: irr, Mult: sf.Mult})
		}
	}
	paSortFactors(out)
	return out
}

// --- invariants ------------------------------------------------------------

// Resultant returns the resultant of p and q, computed exactly over the
// rationals via the Euclidean remainder sequence. The resultant is zero exactly
// when p and q share a non-constant common factor.
func (p *Poly) Resultant(q *Poly) Expr {
	_ = paBinVar(p, q)
	return newRational(paResultant(p.coeffs, q.coeffs))
}

// Discriminant returns the discriminant of the polynomial, defined as
// (-1)^{n(n-1)/2} / a_n * Resultant(p, p'), where n is the degree and a_n is the
// leading coefficient. Polynomials of degree less than one yield Int(0).
func (p *Poly) Discriminant() Expr {
	n := p.Degree()
	if n < 1 {
		return Int(0)
	}
	res := paResultant(p.coeffs, p.Derivative().coeffs)
	an := p.coeffs[len(p.coeffs)-1]
	d := new(big.Rat).Quo(res, an)
	if (n*(n-1)/2)&1 == 1 {
		d.Neg(d)
	}
	return newRational(d)
}

// --- partial fractions -----------------------------------------------------

// PartialFractions computes the partial-fraction decomposition of num/den over
// the rationals. It returns the polynomial (improper) part from Euclidean
// division and a list of proper terms Coeff/(Denom^Power), one per distinct
// irreducible factor of den and each of its powers with a non-zero numerator.
// The denominator is factored with [Poly.Factor] and the residues are solved as
// an exact rational linear system. It returns an error if den is the zero
// polynomial or if den has an irreducible factor of degree three or higher that
// cannot be reduced.
func PartialFractions(num, den *Poly) (poly *Poly, terms []PFTerm, err error) {
	if den.IsZero() {
		return nil, nil, errors.New("algebra: PartialFractions with zero denominator")
	}
	v := paBinVar(num, den)
	quo, rem, e := num.DivMod(den)
	if e != nil {
		return nil, nil, e
	}
	poly = quo
	if rem.IsZero() {
		return poly, nil, nil
	}
	// Distinct irreducible factors of the denominator (skip the constant part).
	type ff struct {
		base *Poly
		mult int
	}
	var irr []ff
	for _, f := range den.Factor() {
		if f.Base.Degree() < 1 {
			continue
		}
		if f.Base.Degree() >= 3 {
			return nil, nil, errors.New("algebra: PartialFractions cannot handle an irreducible factor of degree >= 3")
		}
		irr = append(irr, ff{f.Base, f.Mult})
	}
	// Work over the monic denominator; fold the leading coefficient into the
	// numerator so the residues match the monic factorization exactly.
	denM := den.Monic()
	lcInv := paRatInv(den.coeffs[len(den.coeffs)-1])
	remM := rem.Scale(newRational(lcInv))
	n := denM.Degree()

	// Enumerate the unknown numerator coefficients and their basis polynomials.
	type unk struct{ fi, power, j int }
	var unks []unk
	var basis [][]*big.Rat
	for fi := range irr {
		df := irr[fi].base.Degree()
		for power := 1; power <= irr[fi].mult; power++ {
			cof := paCofactor(denM, irr[fi].base, power) // denM / base^power
			for j := 0; j < df; j++ {
				basis = append(basis, paShift(cof, j)) // x^j * cof
				unks = append(unks, unk{fi, power, j})
			}
		}
	}
	// Build and solve the square rational linear system M x = b.
	M := make([][]*big.Rat, n)
	for i := 0; i < n; i++ {
		M[i] = make([]*big.Rat, len(unks))
		for u := range unks {
			if i < len(basis[u]) {
				M[i][u] = new(big.Rat).Set(basis[u][i])
			} else {
				M[i][u] = big.NewRat(0, 1)
			}
		}
	}
	b := make([]*big.Rat, n)
	for i := 0; i < n; i++ {
		if i < len(remM.coeffs) {
			b[i] = new(big.Rat).Set(remM.coeffs[i])
		} else {
			b[i] = big.NewRat(0, 1)
		}
	}
	sol, err := gaussSolve(M, b)
	if err != nil {
		return nil, nil, err
	}
	// Regroup the solution into a numerator per (factor, power).
	type key struct{ fi, power int }
	numer := map[key][]*big.Rat{}
	for idx, u := range unks {
		k := key{u.fi, u.power}
		arr := numer[k]
		if arr == nil {
			arr = make([]*big.Rat, irr[u.fi].base.Degree())
			for t := range arr {
				arr[t] = big.NewRat(0, 1)
			}
			numer[k] = arr
		}
		arr[u.j] = new(big.Rat).Set(sol[idx])
	}
	for fi := range irr {
		for power := 1; power <= irr[fi].mult; power++ {
			np := &Poly{v: v, coeffs: paTrim(numer[key{fi, power}])}
			if np.IsZero() {
				continue
			}
			terms = append(terms, PFTerm{
				Coeff: np.Expr(),
				Denom: paClone(irr[fi].base),
				Power: power,
			})
		}
	}
	return poly, terms, nil
}

// ApartExpr computes the partial-fraction decomposition of the rational
// expression e in the symbol v and rebuilds it as an [Expr]. It splits e into
// numerator and denominator, forms the corresponding polynomials with
// [PolyFrom], applies [PartialFractions] and reassembles the result through
// [Add], [Mul] and [Pow]. It returns an error if e is not a rational function of
// v or if the decomposition cannot be carried out. v must be a [Symbol].
func ApartExpr(e Expr, v Expr) (Expr, error) {
	if _, ok := v.(*Symbol); !ok {
		return nil, errors.New("algebra: ApartExpr requires a symbol")
	}
	numE, denE := numDenom(Simplify(e))
	num, err := PolyFrom(numE, v)
	if err != nil {
		return nil, err
	}
	den, err := PolyFrom(denE, v)
	if err != nil {
		return nil, err
	}
	poly, terms, err := PartialFractions(num, den)
	if err != nil {
		return nil, err
	}
	parts := []Expr{poly.Expr()}
	for _, t := range terms {
		parts = append(parts, Mul(t.Coeff, Pow(t.Denom.Expr(), Int(int64(-t.Power)))))
	}
	return Add(parts...), nil
}

// --- unexported helpers ----------------------------------------------------

// paToRat converts an integer or rational expression to a fresh *big.Rat,
// reporting false for any other kind of expression.
func paToRat(e Expr) (*big.Rat, bool) {
	r, ok := toRat(e)
	if !ok {
		return nil, false
	}
	return new(big.Rat).Set(r), true
}

// paRatInv returns 1/r as a fresh *big.Rat.
func paRatInv(r *big.Rat) *big.Rat { return new(big.Rat).Inv(r) }

// paClone returns an independent copy of p.
func paClone(p *Poly) *Poly {
	out := make([]*big.Rat, len(p.coeffs))
	for i, c := range p.coeffs {
		out[i] = new(big.Rat).Set(c)
	}
	return &Poly{v: p.v, coeffs: out}
}

// paBinVar returns the variable name shared by a binary operation on p and q.
// Constant (degree <= 0) operands are variable-agnostic; it panics if both
// operands are non-constant in different variables.
func paBinVar(p, q *Poly) string {
	if p.Degree() >= 1 && q.Degree() >= 1 && p.v != q.v {
		panic("algebra: Poly variable mismatch: " + p.v + " vs " + q.v)
	}
	if p.Degree() >= 1 {
		return p.v
	}
	if q.Degree() >= 1 {
		return q.v
	}
	return p.v
}

// paTrim drops trailing (highest-degree) zero coefficients, returning the empty
// slice for the zero polynomial.
func paTrim(c []*big.Rat) []*big.Rat {
	n := len(c)
	for n > 0 && c[n-1].Sign() == 0 {
		n--
	}
	return c[:n]
}

// paTrimCopy returns a fresh trimmed copy of c.
func paTrimCopy(c []*big.Rat) []*big.Rat {
	t := paTrim(c)
	out := make([]*big.Rat, len(t))
	for i := range t {
		out[i] = new(big.Rat).Set(t[i])
	}
	return out
}

// paAddC returns the untrimmed coefficient-wise sum of a and b.
func paAddC(a, b []*big.Rat) []*big.Rat {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]*big.Rat, n)
	for i := 0; i < n; i++ {
		s := big.NewRat(0, 1)
		if i < len(a) {
			s.Add(s, a[i])
		}
		if i < len(b) {
			s.Add(s, b[i])
		}
		out[i] = s
	}
	return out
}

// paSubC returns the untrimmed coefficient-wise difference a - b.
func paSubC(a, b []*big.Rat) []*big.Rat {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	out := make([]*big.Rat, n)
	for i := 0; i < n; i++ {
		s := big.NewRat(0, 1)
		if i < len(a) {
			s.Add(s, a[i])
		}
		if i < len(b) {
			s.Sub(s, b[i])
		}
		out[i] = s
	}
	return out
}

// paMulC returns the untrimmed convolution product of a and b.
func paMulC(a, b []*big.Rat) []*big.Rat {
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	out := make([]*big.Rat, len(a)+len(b)-1)
	for i := range out {
		out[i] = big.NewRat(0, 1)
	}
	for i := range a {
		for j := range b {
			out[i+j].Add(out[i+j], new(big.Rat).Mul(a[i], b[j]))
		}
	}
	return out
}

// paDerivC returns the untrimmed derivative coefficients of c.
func paDerivC(c []*big.Rat) []*big.Rat {
	if len(c) <= 1 {
		return nil
	}
	out := make([]*big.Rat, len(c)-1)
	for i := 1; i < len(c); i++ {
		out[i-1] = new(big.Rat).Mul(c[i], new(big.Rat).SetInt64(int64(i)))
	}
	return out
}

// paMonicC returns c divided by its leading coefficient, or the empty slice for
// the zero polynomial.
func paMonicC(c []*big.Rat) []*big.Rat {
	c = paTrim(c)
	if len(c) == 0 {
		return nil
	}
	lead := c[len(c)-1]
	out := make([]*big.Rat, len(c))
	for i := range c {
		out[i] = new(big.Rat).Quo(c[i], lead)
	}
	return out
}

// paDivRat divides num by the non-zero polynomial den, returning untrimmed
// quotient and remainder coefficient slices.
func paDivRat(num, den []*big.Rat) (quo, rem []*big.Rat) {
	num = paTrimCopy(num)
	den = paTrimCopy(den)
	dd := len(den) - 1
	rem = make([]*big.Rat, len(num))
	for i := range num {
		rem[i] = new(big.Rat).Set(num[i])
	}
	if len(num) < len(den) {
		return nil, rem
	}
	quo = make([]*big.Rat, len(num)-len(den)+1)
	for i := range quo {
		quo[i] = big.NewRat(0, 1)
	}
	lead := den[dd]
	for {
		rt := paTrim(rem)
		rd := len(rt) - 1
		if rd < dd {
			break
		}
		coef := new(big.Rat).Quo(rt[rd], lead)
		shift := rd - dd
		quo[shift] = coef
		for i := 0; i <= dd; i++ {
			rem[shift+i].Sub(rem[shift+i], new(big.Rat).Mul(coef, den[i]))
		}
	}
	return quo, rem
}

// paGCDC returns the monic Euclidean GCD of a and b, or the empty slice when
// both are zero.
func paGCDC(a, b []*big.Rat) []*big.Rat {
	a = paTrimCopy(a)
	b = paTrimCopy(b)
	for len(b) > 0 {
		_, r := paDivRat(a, b)
		a = b
		b = paTrim(r)
	}
	return paMonicC(a)
}

// paMustQuo returns the exact quotient a/b, panicking if the division does not
// come out even. It is used by the square-free algorithm where divisibility is
// guaranteed.
func paMustQuo(a, b *Poly) *Poly {
	q, _, err := a.DivMod(b)
	if err != nil {
		panic("algebra: Poly square-free division by zero")
	}
	return q
}

// paFactorSquareFree factors a monic square-free polynomial s into monic linear
// factors (from its rational roots) plus at most one remaining higher-degree
// factor with no rational root. Factors are returned in a deterministic order.
func paFactorSquareFree(s *Poly) []*Poly {
	v := s.v
	core := paTrimCopy(s.coeffs)
	var out []*Poly
	for {
		roots := rationalRoots(core)
		if len(roots) == 0 {
			break
		}
		paSortRats(roots)
		for _, r := range roots {
			out = append(out, &Poly{v: v, coeffs: []*big.Rat{new(big.Rat).Neg(r), big.NewRat(1, 1)}})
			core = deflate(core, r)
		}
	}
	core = paTrim(core)
	if len(core)-1 >= 1 {
		out = append(out, &Poly{v: v, coeffs: paTrimCopy(core)})
	}
	return out
}

// paCofactor returns denM/base^power as coefficients; base^power must divide
// denM exactly.
func paCofactor(denM, base *Poly, power int) []*big.Rat {
	acc := base
	for k := 1; k < power; k++ {
		acc = acc.MulP(base)
	}
	quo, _, _ := denM.DivMod(acc)
	return quo.coeffs
}

// paShift returns x^j * c (prepending j zero coefficients).
func paShift(c []*big.Rat, j int) []*big.Rat {
	out := make([]*big.Rat, len(c)+j)
	for i := range out {
		out[i] = big.NewRat(0, 1)
	}
	for i := range c {
		out[i+j] = new(big.Rat).Set(c[i])
	}
	return out
}

// paResultant returns the resultant of the polynomials a and b via the
// Euclidean remainder sequence over the rationals.
func paResultant(a0, b0 []*big.Rat) *big.Rat {
	a := paTrimCopy(a0)
	b := paTrimCopy(b0)
	if len(a) == 0 || len(b) == 0 {
		return big.NewRat(0, 1)
	}
	res := big.NewRat(1, 1)
	for {
		m := len(a) - 1
		n := len(b) - 1
		if n == 0 {
			res.Mul(res, paRatPow(b[0], m))
			return res
		}
		if m < n {
			if (m*n)&1 == 1 {
				res.Neg(res)
			}
			a, b = b, a
			continue
		}
		_, r := paDivRat(a, b)
		r = paTrim(r)
		if len(r) == 0 {
			return big.NewRat(0, 1)
		}
		d := len(r) - 1
		if (m*n)&1 == 1 {
			res.Neg(res)
		}
		res.Mul(res, paRatPow(b[len(b)-1], m-d))
		a = b
		b = paTrimCopy(r)
	}
}

// paRatPow returns r raised to the non-negative integer power n.
func paRatPow(r *big.Rat, n int) *big.Rat {
	out := big.NewRat(1, 1)
	base := new(big.Rat).Set(r)
	for k := 0; k < n; k++ {
		out.Mul(out, base)
	}
	return out
}

// paSortRats sorts a slice of rationals in ascending order.
func paSortRats(rs []*big.Rat) {
	sort.SliceStable(rs, func(i, j int) bool { return rs[i].Cmp(rs[j]) < 0 })
}

// paSortFactors orders factors deterministically by base degree, then base
// string, then multiplicity.
func paSortFactors(fs []PolyFactor) {
	sort.SliceStable(fs, func(i, j int) bool {
		di, dj := fs[i].Base.Degree(), fs[j].Base.Degree()
		if di != dj {
			return di < dj
		}
		si, sj := fs[i].Base.String(), fs[j].Base.String()
		if si != sj {
			return si < sj
		}
		return fs[i].Mult < fs[j].Mult
	})
}
