package algebra

import (
	"errors"
	"math/big"
)

// Solve returns the real solutions of the equation e == 0 for the symbol v.
// It handles expressions that are polynomials in v of degree 1 (linear) or 2
// (quadratic). Quadratic roots are returned exactly, using an [Integer] or
// [Rational] when the discriminant is a perfect square and a symbolic sqrt
// otherwise; a repeated root is returned once. Higher-degree or non-polynomial
// equations return an error. v must be a [Symbol].
//
// To solve lhs == rhs, pass Solve(lhs.Add(rhs.Mul(Int(-1))), v), i.e. move
// everything to one side.
func Solve(e, v Expr) ([]Expr, error) {
	s, ok := v.(*Symbol)
	if !ok {
		return nil, errors.New("algebra: Solve requires a symbol")
	}
	coeffs, ok := polyCoeffs(e, s.Name)
	if !ok {
		return nil, errors.New("algebra: cannot solve non-polynomial equation")
	}
	deg := len(coeffs) - 1
	for deg > 0 && isZero(coeffs[deg]) {
		deg--
	}
	switch deg {
	case 0:
		if isZero(coeffs[0]) {
			return nil, errors.New("algebra: equation holds for all values")
		}
		return nil, errors.New("algebra: no solution")
	case 1:
		// c1*x + c0 = 0  ->  x = -c0/c1.
		root := Simplify(Mul(neg(coeffs[0]), Pow(coeffs[1], Int(-1))))
		return []Expr{root}, nil
	case 2:
		return solveQuadratic(coeffs[2], coeffs[1], coeffs[0], v), nil
	default:
		return nil, errors.New("algebra: only linear and quadratic equations are supported")
	}
}

// solveQuadratic solves a*x^2 + b*x + c = 0 using the quadratic formula.
func solveQuadratic(a, b, c, v Expr) []Expr {
	disc := Simplify(Add(Pow(b, Int(2)), Mul(Int(-4), a, c)))
	root := sqrtExpr(disc)
	denom := Mul(Int(2), a)
	r1 := Simplify(Mul(Add(neg(b), root), Pow(denom, Int(-1))))
	r2 := Simplify(Mul(Add(neg(b), neg(root)), Pow(denom, Int(-1))))
	if r1.Equal(r2) {
		return []Expr{r1}
	}
	out := []Expr{r1, r2}
	sortExprs(out)
	return out
}

// sqrtExpr returns the square root of a numeric expression in simplest form:
// exact when it is a perfect square, and with the largest square factor pulled
// out of the radical otherwise (e.g. sqrt(8) -> 2*sqrt(2)).
func sqrtExpr(e Expr) Expr {
	switch n := e.(type) {
	case *Integer:
		if n.Val.Sign() >= 0 {
			return Sqrt(n)
		}
	case *Rational:
		if n.Val.Sign() > 0 {
			// sqrt(p/q) = sqrt(p*q)/q.
			pq := new(big.Int).Mul(n.Val.Num(), n.Val.Denom())
			out, rad := extractSquare(pq)
			return Mul(newRational(new(big.Rat).SetFrac(out, n.Val.Denom())), Sqrt(newInteger(rad)))
		}
	}
	return Sqrt(e)
}

// extractSquare factors n = out^2 * rad with rad square-free, returning out and
// rad. n must be non-negative.
func extractSquare(n *big.Int) (out, rad *big.Int) {
	out = big.NewInt(1)
	rad = new(big.Int).Set(n)
	if n.Sign() == 0 {
		return big.NewInt(0), big.NewInt(1)
	}
	f := big.NewInt(2)
	fsq := big.NewInt(4)
	for fsq.Cmp(rad) <= 0 {
		for {
			q := new(big.Int)
			m := new(big.Int)
			q.DivMod(rad, fsq, m)
			if m.Sign() != 0 {
				break
			}
			rad.Set(q)
			out.Mul(out, f)
		}
		f.Add(f, big.NewInt(1))
		fsq.Mul(f, f)
	}
	return out, rad
}

// intSqrt returns the integer square root of n and whether n is a perfect
// square. n must be non-negative.
func intSqrt(n *big.Int) (*big.Int, bool) {
	if n.Sign() < 0 {
		return nil, false
	}
	r := new(big.Int).Sqrt(n)
	sq := new(big.Int).Mul(r, r)
	return r, sq.Cmp(n) == 0
}

// polyCoeffs returns the coefficients of e viewed as a polynomial in the
// symbol named name, with index i holding the coefficient of x^i. It reports
// false if e is not polynomial in x.
func polyCoeffs(e Expr, name string) ([]Expr, bool) {
	ex := Simplify(Expand(e))
	m := map[int]Expr{}
	maxd := 0
	for _, t := range termsOf(ex) {
		k, c, ok := monomial(t, name)
		if !ok {
			return nil, false
		}
		if k > maxd {
			maxd = k
		}
		if prev, ok := m[k]; ok {
			m[k] = Add(prev, c)
		} else {
			m[k] = c
		}
	}
	out := make([]Expr, maxd+1)
	for i := range out {
		if c, ok := m[i]; ok {
			out[i] = Simplify(c)
		} else {
			out[i] = Int(0)
		}
	}
	return out, true
}

// monomial decomposes a single term into (degree in x, coefficient). It
// reports false if the term is not a monomial in x (e.g. x appears inside a
// function).
func monomial(t Expr, name string) (int, Expr, bool) {
	deg := 0
	var coeff []Expr
	for _, f := range factorsOf(t) {
		if s, ok := f.(*Symbol); ok && s.Name == name {
			deg++
			continue
		}
		if p, ok := f.(*power); ok {
			if bs, ok := p.base.(*Symbol); ok && bs.Name == name {
				if n, ok := p.exp.(*Integer); ok && n.Val.Sign() >= 0 && n.Val.IsInt64() {
					deg += int(n.Val.Int64())
					continue
				}
			}
		}
		if containsSym(f, name) {
			return 0, nil, false
		}
		coeff = append(coeff, f)
	}
	return deg, Mul(coeff...), true
}

// Factor factors a univariate polynomial e in the symbol v of degree 1 or 2
// into a leading coefficient times monic linear factors, when its roots are
// rational or expressible with a single sqrt. If e cannot be factored this
// way it is returned unchanged (canonicalized). This is a convenience helper,
// not a general factorization algorithm.
func Factor(e, v Expr) Expr {
	s, ok := v.(*Symbol)
	if !ok {
		return Simplify(e)
	}
	coeffs, ok := polyCoeffs(e, s.Name)
	if !ok {
		return Simplify(e)
	}
	deg := len(coeffs) - 1
	for deg > 0 && isZero(coeffs[deg]) {
		deg--
	}
	if deg != 2 {
		return Simplify(e)
	}
	a := coeffs[2]
	roots := solveQuadratic(a, coeffs[1], coeffs[0], v)
	factors := []Expr{a}
	for _, r := range roots {
		factors = append(factors, Add(v, neg(r)))
	}
	if len(roots) == 1 {
		// (x-r) appears squared.
		factors = append(factors, Add(v, neg(roots[0])))
	}
	return Mul(factors...)
}

// Collect rewrites e as a canonical polynomial in the symbol v, grouping like
// powers of v. If e is not polynomial in v it is simplified and returned.
func Collect(e, v Expr) Expr {
	s, ok := v.(*Symbol)
	if !ok {
		return Simplify(e)
	}
	coeffs, ok := polyCoeffs(e, s.Name)
	if !ok {
		return Simplify(e)
	}
	terms := make([]Expr, 0, len(coeffs))
	for i, c := range coeffs {
		terms = append(terms, Mul(c, Pow(v, Int(int64(i)))))
	}
	return Add(terms...)
}

// neg returns -e in canonical form.
func neg(e Expr) Expr { return Mul(Int(-1), e) }
