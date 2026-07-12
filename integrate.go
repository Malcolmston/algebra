package algebra

import "math/big"

// This file holds the extended integration strategies: inverse-square-root
// forms yielding asin/asinh/acosh, integration of rational functions by
// polynomial division plus partial fractions (distinct rational poles) or the
// arctangent form for an irreducible quadratic, and integration by parts for
// products of a polynomial with exp/sin/cos. It also defines the
// polynomial-over-the-rationals helpers shared with the cubic/quartic solvers.

// --- rational polynomial helpers -------------------------------------------

// ratCoeffs converts a coefficient slice of numeric expressions to big.Rat,
// reporting false if any entry is not an integer or rational.
func ratCoeffs(cs []Expr) ([]*big.Rat, bool) {
	out := make([]*big.Rat, len(cs))
	for i, c := range cs {
		r, ok := toRat(c)
		if !ok {
			return nil, false
		}
		out[i] = new(big.Rat).Set(r)
	}
	return out, true
}

// polyRatCoeffs returns the coefficients of e as a polynomial in name over the
// rationals (index i is the coefficient of x^i).
func polyRatCoeffs(e Expr, name string) ([]*big.Rat, bool) {
	cs, ok := polyCoeffs(e, name)
	if !ok {
		return nil, false
	}
	return ratCoeffs(cs)
}

// trimRat drops leading (highest-degree) zero coefficients, keeping at least
// one entry.
func trimRat(c []*big.Rat) []*big.Rat {
	n := len(c)
	for n > 1 && c[n-1].Sign() == 0 {
		n--
	}
	return c[:n]
}

func ratDegree(c []*big.Rat) int { return len(trimRat(c)) - 1 }

// ratPolyEval evaluates the polynomial c at x by Horner's method.
func ratPolyEval(c []*big.Rat, x *big.Rat) *big.Rat {
	acc := big.NewRat(0, 1)
	for i := len(c) - 1; i >= 0; i-- {
		acc.Mul(acc, x)
		acc.Add(acc, c[i])
	}
	return acc
}

// ratPolyDeriv returns the derivative coefficients of c.
func ratPolyDeriv(c []*big.Rat) []*big.Rat {
	if len(c) <= 1 {
		return []*big.Rat{big.NewRat(0, 1)}
	}
	out := make([]*big.Rat, len(c)-1)
	for i := 1; i < len(c); i++ {
		out[i-1] = new(big.Rat).Mul(c[i], new(big.Rat).SetInt64(int64(i)))
	}
	return out
}

// ratPolyDiv divides num by den, returning quotient and remainder coefficient
// slices (index i is the coefficient of x^i).
func ratPolyDiv(num, den []*big.Rat) (quot, rem []*big.Rat) {
	num = trimRat(num)
	den = trimRat(den)
	dd := len(den) - 1
	rem = make([]*big.Rat, len(num))
	for i := range num {
		rem[i] = new(big.Rat).Set(num[i])
	}
	qlen := len(num) - len(den) + 1
	if qlen < 1 {
		qlen = 1
	}
	quot = make([]*big.Rat, qlen)
	for i := range quot {
		quot[i] = big.NewRat(0, 1)
	}
	lead := den[dd]
	for {
		rt := trimRat(rem)
		rd := len(rt) - 1
		if rd < dd || (rd == 0 && rt[0].Sign() == 0) {
			break
		}
		coef := new(big.Rat).Quo(rt[rd], lead)
		shift := rd - dd
		if shift < len(quot) {
			quot[shift] = coef
		}
		for i := 0; i <= dd; i++ {
			rem[shift+i].Sub(rem[shift+i], new(big.Rat).Mul(coef, den[i]))
		}
	}
	return trimRat(quot), trimRat(rem)
}

// deflate divides c by (x - r) exactly (r must be a root), returning the
// quotient coefficients.
func deflate(c []*big.Rat, r *big.Rat) []*big.Rat {
	c = trimRat(c)
	d := len(c) - 1
	if d < 1 {
		return []*big.Rat{big.NewRat(0, 1)}
	}
	q := make([]*big.Rat, d)
	q[d-1] = new(big.Rat).Set(c[d])
	for i := d - 1; i >= 1; i-- {
		q[i-1] = new(big.Rat).Add(c[i], new(big.Rat).Mul(r, q[i]))
	}
	return q
}

func gcdBig(a, b *big.Int) *big.Int {
	return new(big.Int).GCD(nil, nil, new(big.Int).Abs(a), new(big.Int).Abs(b))
}

func lcmBig(a, b *big.Int) *big.Int {
	if a.Sign() == 0 || b.Sign() == 0 {
		return big.NewInt(1)
	}
	g := gcdBig(a, b)
	return new(big.Int).Abs(new(big.Int).Mul(new(big.Int).Div(a, g), b))
}

// divisorsBig returns the positive divisors of |n| (n != 0).
func divisorsBig(n *big.Int) []*big.Int {
	n = new(big.Int).Abs(n)
	if n.Sign() == 0 {
		return []*big.Int{big.NewInt(1)}
	}
	var out []*big.Int
	i := big.NewInt(1)
	for new(big.Int).Mul(i, i).Cmp(n) <= 0 {
		if new(big.Int).Rem(n, i).Sign() == 0 {
			out = append(out, new(big.Int).Set(i))
			j := new(big.Int).Div(n, i)
			if j.Cmp(i) != 0 {
				out = append(out, j)
			}
		}
		i.Add(i, big.NewInt(1))
	}
	return out
}

// rationalRoots returns the distinct rational roots of the polynomial c using
// the rational-root theorem.
func rationalRoots(c []*big.Rat) []*big.Rat {
	c = trimRat(c)
	d := len(c) - 1
	if d < 1 {
		return nil
	}
	lcm := big.NewInt(1)
	for _, x := range c {
		lcm = lcmBig(lcm, x.Denom())
	}
	ic := make([]*big.Int, len(c))
	for i, x := range c {
		t := new(big.Rat).Mul(x, new(big.Rat).SetInt(lcm))
		ic[i] = new(big.Int).Set(t.Num())
	}
	var roots []*big.Rat
	seen := map[string]bool{}
	add := func(r *big.Rat) {
		k := r.RatString()
		if !seen[k] && ratPolyEval(c, r).Sign() == 0 {
			seen[k] = true
			roots = append(roots, new(big.Rat).Set(r))
		}
	}
	low := 0
	for low <= d && ic[low].Sign() == 0 {
		low++
	}
	if low > 0 {
		add(big.NewRat(0, 1))
	}
	if low > d {
		return roots
	}
	ps := divisorsBig(ic[low])
	qs := divisorsBig(ic[d])
	for _, p := range ps {
		for _, q := range qs {
			r := new(big.Rat).SetFrac(p, q)
			add(new(big.Rat).Set(r))
			add(new(big.Rat).Neg(r))
		}
	}
	return roots
}

// --- integration strategies ------------------------------------------------

// numDenom splits e into numerator and denominator by collecting factors with
// negative integer exponents into the denominator.
func numDenom(e Expr) (num, den Expr) {
	num, den = Int(1), Int(1)
	for _, f := range factorsOf(e) {
		b, ex := splitPow(f)
		if n, ok := ex.(*Integer); ok && n.Val.Sign() < 0 {
			den = Mul(den, Pow(b, neg(ex)))
		} else {
			num = Mul(num, f)
		}
	}
	return num, den
}

// integratePolyRat integrates a polynomial given by rational coefficients.
func integratePolyRat(c []*big.Rat, v Expr) Expr {
	terms := make([]Expr, 0, len(c))
	for i, coef := range c {
		if coef.Sign() == 0 {
			continue
		}
		np1 := int64(i + 1)
		cc := new(big.Rat).Quo(coef, big.NewRat(np1, 1))
		terms = append(terms, Mul(newRational(cc), Pow(v, Int(np1))))
	}
	return Add(terms...)
}

// integrateRational integrates a rational function of the symbol name via
// polynomial division followed by partial fractions over distinct rational
// poles, or the arctangent/log form for an irreducible quadratic denominator.
// It returns nil for integrands it does not resolve.
func integrateRational(e Expr, name string, v Expr) Expr {
	num, den := numDenom(e)
	if !containsSym(den, name) {
		return nil
	}
	nc, ok := polyRatCoeffs(num, name)
	if !ok {
		return nil
	}
	dc, ok := polyRatCoeffs(den, name)
	if !ok {
		return nil
	}
	nc = trimRat(nc)
	dc = trimRat(dc)
	if ratDegree(dc) < 1 {
		return nil
	}
	var parts []Expr
	// Polynomial part for improper fractions.
	if ratDegree(nc) >= ratDegree(dc) {
		q, r := ratPolyDiv(nc, dc)
		parts = append(parts, integratePolyRat(q, v))
		nc = trimRat(r)
	}
	if len(nc) == 1 && nc[0].Sign() == 0 {
		return Simplify(Add(parts...))
	}
	// Peel off distinct rational roots, tracking multiplicity.
	roots := rationalRoots(dc)
	core := append([]*big.Rat(nil), dc...)
	mult := map[string]int{}
	for _, r := range roots {
		for ratPolyEval(core, r).Sign() == 0 && ratDegree(core) >= 1 {
			core = deflate(core, r)
			mult[r.RatString()]++
		}
	}
	switch {
	case ratDegree(core) == 0:
		// Fully split into linear factors; require simple poles.
		for _, m := range mult {
			if m != 1 {
				return nil
			}
		}
		dcDeriv := ratPolyDeriv(dc)
		for _, r := range roots {
			a := new(big.Rat).Quo(ratPolyEval(nc, r), ratPolyEval(dcDeriv, r))
			parts = append(parts, Mul(newRational(a), Log(Add(v, neg(newRational(new(big.Rat).Set(r)))))))
		}
		return Simplify(Add(parts...))
	case ratDegree(dc) == 2 && len(roots) == 0:
		if r := integrateIrreducibleQuadratic(nc, dc, v); r != nil {
			parts = append(parts, r)
			return Simplify(Add(parts...))
		}
	}
	return nil
}

// integrateIrreducibleQuadratic integrates (p*x+q)/(a*x^2+b*x+c) when the
// denominator has negative discriminant, giving a log plus an arctangent.
func integrateIrreducibleQuadratic(nc, dc []*big.Rat, v Expr) Expr {
	a, b, c := dc[2], dc[1], dc[0]
	p := big.NewRat(0, 1)
	q := big.NewRat(0, 1)
	if len(nc) >= 1 {
		q = nc[0]
	}
	if len(nc) >= 2 {
		p = nc[1]
	}
	// D = 4ac - b^2 must be positive for an irreducible quadratic.
	D := new(big.Rat).Sub(
		new(big.Rat).Mul(big.NewRat(4, 1), new(big.Rat).Mul(a, c)),
		new(big.Rat).Mul(b, b))
	if D.Sign() <= 0 {
		return nil
	}
	den := Add(Mul(newRational(new(big.Rat).Set(a)), Pow(v, Int(2))),
		Mul(newRational(new(big.Rat).Set(b)), v), newRational(new(big.Rat).Set(c)))
	twoA := new(big.Rat).Mul(big.NewRat(2, 1), a)
	p2a := new(big.Rat).Quo(p, twoA)
	term1 := Mul(newRational(p2a), Log(den))
	k := new(big.Rat).Sub(q, new(big.Rat).Mul(p2a, b))
	sqrtD := sqrtExpr(newRational(new(big.Rat).Set(D)))
	inner := Mul(Add(Mul(newRational(new(big.Rat).Set(twoA)), v), newRational(new(big.Rat).Set(b))),
		Pow(sqrtD, Int(-1)))
	atanPart := Mul(Int(2), Pow(sqrtD, Int(-1)), Atan(inner))
	term2 := Mul(newRational(k), atanPart)
	return Add(term1, term2)
}

// integrateInvSqrtQuadratic integrates 1/sqrt(quadratic) forms, recognising
// 1/sqrt(c-x^2) -> asin(x/sqrt c), 1/sqrt(x^2+c) -> asinh(x/sqrt c) and
// 1/sqrt(x^2-c) -> acosh(x/sqrt c) (no linear term in the quadratic).
func integrateInvSqrtQuadratic(e Expr, name string, v Expr) Expr {
	pw, ok := e.(*power)
	if !ok {
		return nil
	}
	var base Expr
	switch {
	case pw.exp.Equal(Int(-1)):
		f, ok := pw.base.(*fn)
		if !ok || f.name != "sqrt" {
			return nil
		}
		base = f.arg
	case pw.exp.Equal(Rat(-1, 2)):
		base = pw.base
	default:
		return nil
	}
	bc, ok := polyRatCoeffs(base, name)
	if !ok {
		return nil
	}
	bc = trimRat(bc)
	if ratDegree(bc) != 2 || bc[1].Sign() != 0 {
		return nil
	}
	c2, c0 := bc[2], bc[0]
	switch {
	case c2.Cmp(big.NewRat(-1, 1)) == 0 && c0.Sign() > 0:
		s := sqrtExpr(newRational(new(big.Rat).Set(c0)))
		return Asin(Mul(v, Pow(s, Int(-1))))
	case c2.Cmp(big.NewRat(1, 1)) == 0 && c0.Sign() > 0:
		s := sqrtExpr(newRational(new(big.Rat).Set(c0)))
		return Asinh(Mul(v, Pow(s, Int(-1))))
	case c2.Cmp(big.NewRat(1, 1)) == 0 && c0.Sign() < 0:
		s := sqrtExpr(newRational(new(big.Rat).Neg(c0)))
		return Acosh(Mul(v, Pow(s, Int(-1))))
	}
	return nil
}

// integrateByParts integrates a product of a polynomial in name with a single
// exp/sin/cos kernel whose argument is linear, using tabular integration by
// parts. It returns nil when the integrand is not of that shape.
func integrateByParts(e Expr, name string, v Expr) Expr {
	factors := factorsOf(e)
	kernelIdx := -1
	for i, f := range factors {
		fnn, ok := f.(*fn)
		if !ok {
			continue
		}
		if fnn.name != "exp" && fnn.name != "sin" && fnn.name != "cos" {
			continue
		}
		if _, _, ok := linearCoeffs(fnn.arg, name, v); ok {
			kernelIdx = i
			break
		}
	}
	if kernelIdx < 0 {
		return nil
	}
	rest := make([]Expr, 0, len(factors)-1)
	for i, f := range factors {
		if i != kernelIdx {
			rest = append(rest, f)
		}
	}
	poly := Mul(rest...)
	pc, ok := polyCoeffs(poly, name)
	if !ok || len(pc)-1 < 1 {
		return nil
	}
	kernel := factors[kernelIdx]

	deriv := poly
	anti := kernel
	sign := int64(1)
	var terms []Expr
	for i := 0; i < 100; i++ {
		anti = integrate(anti, name, v)
		if _, unresolved := anti.(*integral); unresolved {
			return nil
		}
		terms = append(terms, Mul(Int(sign), deriv, anti))
		deriv = Simplify(diff(deriv, name))
		sign = -sign
		if isZero(deriv) {
			break
		}
	}
	return Add(terms...)
}
