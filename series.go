package algebra

import (
	"math"
	"math/big"
)

// This file adds the higher-level calculus operations Limit, Series, Summation
// and Product, plus the bigOp node used to carry an unevaluated summation,
// product or limit when no closed form is found.

// bigOp is an unevaluated summation, product or limit. kind is "Sum",
// "Product" or "Limit"; index is the bound variable name; for sums and products
// lo and hi are the bounds, while for limits lo holds the target and hi is
// unused.
type bigOp struct {
	builders
	kind  string
	body  Expr
	index string
	lo    Expr
	hi    Expr
}

func newBigOp(kind string, body Expr, index string, lo, hi Expr) Expr {
	e := &bigOp{kind: kind, body: body, index: index, lo: lo, hi: hi}
	e.self = e
	return e
}

func newLimit(body Expr, index string, to Expr) Expr {
	return newBigOp("Limit", body, index, to, Int(0))
}

// String renders a summation, product or limit.
func (b *bigOp) String() string {
	if b.kind == "Limit" {
		return "Limit(" + b.body.String() + ", " + b.index + " -> " + b.lo.String() + ")"
	}
	return b.kind + "(" + b.body.String() + ", " + b.index + "=" + b.lo.String() + ".." + b.hi.String() + ")"
}

// Equal reports structural equality of two bigOp nodes.
func (b *bigOp) Equal(o Expr) bool {
	x, ok := o.(*bigOp)
	return ok && b.kind == x.kind && b.index == x.index &&
		b.body.Equal(x.body) && b.lo.Equal(x.lo) && b.hi.Equal(x.hi)
}

// substBigOp substitutes into the bounds (and body, when name is not the bound
// index) of a summation/product/limit.
func substBigOp(b *bigOp, name string, val Expr) Expr {
	lo := subst(b.lo, name, val)
	hi := subst(b.hi, name, val)
	body := b.body
	if b.index != name {
		body = subst(body, name, val)
	}
	return newBigOp(b.kind, body, b.index, lo, hi)
}

// --- Series ----------------------------------------------------------------

// Series returns the Taylor expansion of e about sym = at, truncated to n terms
// (degrees 0 through n-1). With at = 0 this is the Maclaurin series. The result
// is a polynomial in (sym - at); no remainder term is attached. sym must be a
// [Symbol] and n must be positive.
func Series(e, sym, at Expr, n int) Expr {
	s, ok := sym.(*Symbol)
	if !ok || n <= 0 {
		return e
	}
	var terms []Expr
	deriv := e
	fact := big.NewInt(1)
	for k := 0; k < n; k++ {
		if k > 0 {
			fact.Mul(fact, big.NewInt(int64(k)))
		}
		coef := Simplify(Subs(deriv, sym, at))
		c := Mul(coef, Pow(newInteger(new(big.Int).Set(fact)), Int(-1)))
		var basis Expr
		if k == 0 {
			basis = Int(1)
		} else {
			basis = Pow(Add(sym, neg(at)), Int(int64(k)))
		}
		terms = append(terms, Mul(c, basis))
		deriv = Simplify(diff(deriv, s.Name))
	}
	return Simplify(Add(terms...))
}

// --- Limit -----------------------------------------------------------------

// Limit returns the limit of e as sym approaches to. It evaluates finite limits
// by substitution, applies L'Hôpital's rule to 0/0 and ∞/∞ indeterminate
// forms, and handles limits at ±∞ ([Inf], [NegInf]) of rational functions by
// comparing degrees, with a numeric fallback. When no value can be determined
// an unevaluated Limit node is returned. sym must be a [Symbol].
func Limit(e, sym, to Expr) Expr {
	s, ok := sym.(*Symbol)
	if !ok {
		return e
	}
	return limitRec(Simplify(e), s.Name, to, 0)
}

func limitRec(e Expr, name string, to Expr, depth int) Expr {
	if !containsSym(e, name) {
		return e
	}
	if depth > 12 {
		return newLimit(e, name, to)
	}
	if isInfinite(to) {
		sign := 1
		if c, ok := to.(*Constant); ok && c.Name == "-oo" {
			sign = -1
		}
		return limitInf(e, name, sign, depth)
	}
	toVal, err := Evalf(to)
	if err != nil {
		return newLimit(e, name, to)
	}
	num, den := numDenom(e)
	nv, okn := evalAtSafe(num, name, toVal)
	dv, okd := evalAtSafe(den, name, toVal)
	if !okn || !okd {
		return numericLimit(e, name, toVal, to)
	}
	zeroN := math.Abs(nv) < 1e-9
	zeroD := math.Abs(dv) < 1e-9
	switch {
	case !zeroD:
		return Simplify(Subs(e, Sym(name), to))
	case !zeroN:
		// finite nonzero over zero: diverges.
		return Inf
	default:
		// 0/0: L'Hôpital.
		np := Simplify(diff(num, name))
		dp := Simplify(diff(den, name))
		if isZero(dp) {
			return numericLimit(e, name, toVal, to)
		}
		return limitRec(Simplify(Mul(np, Pow(dp, Int(-1)))), name, to, depth+1)
	}
}

func limitInf(e Expr, name string, sign, depth int) Expr {
	num, den := numDenom(e)
	nc, ok1 := polyRatCoeffs(num, name)
	dc, ok2 := polyRatCoeffs(den, name)
	if ok1 && ok2 {
		nc = trimRat(nc)
		dc = trimRat(dc)
		dn, dd := ratDegree(nc), ratDegree(dc)
		switch {
		case dn < dd:
			return Int(0)
		case dn == dd:
			return newRational(new(big.Rat).Quo(nc[dn], dc[dd]))
		default:
			lead := new(big.Rat).Quo(nc[dn], dc[dd])
			s := lead.Sign()
			if sign < 0 && (dn-dd)%2 == 1 {
				s = -s
			}
			if s >= 0 {
				return Inf
			}
			return NegInf
		}
	}
	return numericLimitInf(e, name, sign)
}

// evalAtSafe evaluates e at name=x, reporting ok only for a finite result.
func evalAtSafe(e Expr, name string, x float64) (float64, bool) {
	v, err := Eval(e, map[string]float64{name: x})
	if err != nil || math.IsNaN(v) || math.IsInf(v, 0) {
		return 0, false
	}
	return v, true
}

// numericLimit estimates a finite limit by sampling e near toVal.
func numericLimit(e Expr, name string, toVal float64, to Expr) Expr {
	var vals []float64
	for _, h := range []float64{1e-3, 1e-5, 1e-7} {
		v, err := Eval(e, map[string]float64{name: toVal + h})
		w, err2 := Eval(e, map[string]float64{name: toVal - h})
		if err != nil || err2 != nil {
			return newLimit(e, name, to)
		}
		vals = append(vals, (v+w)/2)
	}
	if math.Abs(vals[0]-vals[2]) < 1e-4*(1+math.Abs(vals[2])) {
		return Flt(vals[2])
	}
	return newLimit(e, name, to)
}

func numericLimitInf(e Expr, name string, sign int) Expr {
	pts := []float64{1e4, 1e6, 1e8}
	var vals []float64
	for _, p := range pts {
		v, err := Eval(e, map[string]float64{name: float64(sign) * p})
		if err != nil {
			return newLimit(e, name, signInf(sign))
		}
		vals = append(vals, v)
	}
	// A finite integrand that overflows to a signed infinity at the far sample
	// diverges to that infinity (e.g. exp(x)/x as x -> oo). math.Inf here is a
	// genuine overflow of a growing expression, not an oscillation (which would
	// evaluate to NaN and be rejected by the sampling checks below).
	if last := vals[len(vals)-1]; math.IsInf(last, 0) {
		if last > 0 {
			return Inf
		}
		return NegInf
	}
	if math.Abs(vals[1]-vals[2]) < 1e-6*(1+math.Abs(vals[2])) {
		return Flt(vals[2])
	}
	if math.Abs(vals[2]) > math.Abs(vals[1]) {
		if vals[2] > 0 {
			return Inf
		}
		return NegInf
	}
	return newLimit(e, name, signInf(sign))
}

func signInf(sign int) Expr {
	if sign < 0 {
		return NegInf
	}
	return Inf
}

// --- Summation -------------------------------------------------------------

// Summation returns the sum of e over sym running from lo to hi inclusive. It
// finds closed forms for sums whose summand is constant, polynomial (via power
// sums up to degree three) or geometric, and evaluates finite numeric bounds
// directly. Unresolved sums return an unevaluated Sum node. sym must be a
// [Symbol].
func Summation(e, sym, lo, hi Expr) Expr {
	s, ok := sym.(*Symbol)
	if !ok {
		return newBigOp("Sum", e, "?", lo, hi)
	}
	name := s.Name
	e = Simplify(e)
	if !containsSym(e, name) {
		count := Simplify(Add(hi, neg(lo), Int(1)))
		return Simplify(Mul(e, count))
	}
	if r, ok := numericFinite(e, sym, lo, hi, true); ok {
		return r
	}
	if r := sumPolynomial(e, name, lo, hi); r != nil {
		return r
	}
	if r := sumGeometric(e, name, lo, hi); r != nil {
		return r
	}
	return newBigOp("Sum", e, name, lo, hi)
}

// numericFinite evaluates a sum (isSum=true) or product over concrete integer
// bounds by direct iteration.
func numericFinite(e, sym, lo, hi Expr, isSum bool) (Expr, bool) {
	li, ok1 := lo.(*Integer)
	hj, ok2 := hi.(*Integer)
	if !ok1 || !ok2 {
		return nil, false
	}
	if li.Val.Cmp(hj.Val) > 0 {
		if isSum {
			return Int(0), true
		}
		return Int(1), true
	}
	if new(big.Int).Sub(hj.Val, li.Val).Cmp(big.NewInt(100000)) > 0 {
		return nil, false
	}
	acc := Expr(Int(0))
	if !isSum {
		acc = Int(1)
	}
	for k := new(big.Int).Set(li.Val); k.Cmp(hj.Val) <= 0; k.Add(k, big.NewInt(1)) {
		term := Subs(e, sym, newInteger(new(big.Int).Set(k)))
		if isSum {
			acc = Add(acc, term)
		} else {
			acc = Mul(acc, term)
		}
	}
	return Simplify(acc), true
}

// sumPolynomial evaluates a sum whose summand is a polynomial in name using the
// power-sum (Faulhaber) formulas up to degree three.
func sumPolynomial(e Expr, name string, lo, hi Expr) Expr {
	cs, ok := polyCoeffs(e, name)
	if !ok {
		return nil
	}
	if len(cs)-1 > 3 {
		return nil
	}
	var terms []Expr
	for i, c := range cs {
		if isZero(c) {
			continue
		}
		hiPart := powerSum(i, hi)
		loPart := powerSum(i, Add(lo, Int(-1)))
		si := Add(hiPart, neg(loPart))
		terms = append(terms, Mul(c, si))
	}
	return Simplify(Expand(Add(terms...)))
}

// powerSum returns the closed form for sum_{k=1}^{m} k^i (i in 0..3).
func powerSum(i int, m Expr) Expr {
	switch i {
	case 0:
		return m
	case 1:
		return Mul(Rat(1, 2), m, Add(m, Int(1)))
	case 2:
		return Mul(Rat(1, 6), m, Add(m, Int(1)), Add(Mul(Int(2), m), Int(1)))
	case 3:
		half := Mul(Rat(1, 2), m, Add(m, Int(1)))
		return Pow(half, Int(2))
	}
	return nil
}

// sumGeometric evaluates a geometric sum C*base^k with base independent of the
// index.
func sumGeometric(e Expr, name string, lo, hi Expr) Expr {
	coeff := Expr(Int(1))
	var base Expr
	found := false
	for _, f := range factorsOf(e) {
		if !containsSym(f, name) {
			coeff = Mul(coeff, f)
			continue
		}
		b, ex := splitPow(f)
		if ex.Equal(Sym(name)) && !containsSym(b, name) {
			if found {
				return nil
			}
			base = b
			found = true
			continue
		}
		return nil
	}
	if !found {
		return nil
	}
	if isOne(base) {
		return Simplify(Mul(coeff, Add(hi, neg(lo), Int(1))))
	}
	numr := Add(Pow(base, lo), neg(Pow(base, Add(hi, Int(1)))))
	den := Add(Int(1), neg(base))
	return Simplify(Expand(Mul(coeff, numr, Pow(den, Int(-1)))))
}

// --- Product ---------------------------------------------------------------

// Product returns the product of e over sym running from lo to hi inclusive. It
// finds closed forms for constant summands (e^count) and for the summand k
// (giving hi!/(lo-1)!), and evaluates finite numeric bounds directly.
// Unresolved products return an unevaluated Product node. sym must be a
// [Symbol].
func Product(e, sym, lo, hi Expr) Expr {
	s, ok := sym.(*Symbol)
	if !ok {
		return newBigOp("Product", e, "?", lo, hi)
	}
	name := s.Name
	e = Simplify(e)
	if !containsSym(e, name) {
		count := Simplify(Add(hi, neg(lo), Int(1)))
		return Pow(e, count)
	}
	if r, ok := numericFinite(e, sym, lo, hi, false); ok {
		return r
	}
	if e.Equal(sym) {
		return Simplify(Mul(Factorial(hi), Pow(Factorial(Add(lo, Int(-1))), Int(-1))))
	}
	return newBigOp("Product", e, name, lo, hi)
}
