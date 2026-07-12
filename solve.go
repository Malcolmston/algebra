package algebra

import (
	"errors"
	"math"
	"math/big"
	"math/cmplx"
)

// Solve returns the solutions of the equation e == 0 for the symbol v.
//
// It handles polynomials of any degree by first extracting all rational roots
// exactly (via the rational-root theorem and synthetic division) and then
// resolving the remaining factor: a linear or quadratic factor is solved
// exactly with the quadratic formula (complex conjugate roots a±b*I are
// returned when the discriminant is negative), and an irreducible factor of
// degree three or higher is solved numerically with the Durand–Kerner method,
// yielding correct real and complex roots. Repeated roots are returned once.
// Non-polynomial equations return an error. v must be a [Symbol].
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
	if deg == 0 {
		if isZero(coeffs[0]) {
			return nil, errors.New("algebra: equation holds for all values")
		}
		return nil, errors.New("algebra: no solution")
	}
	rc, ok := ratCoeffs(coeffs[:deg+1])
	if !ok {
		// Non-rational coefficients: fall back to the exact low-degree formulas.
		switch deg {
		case 1:
			return []Expr{Simplify(Mul(neg(coeffs[0]), Pow(coeffs[1], Int(-1))))}, nil
		case 2:
			return solveQuadratic(coeffs[2], coeffs[1], coeffs[0], v), nil
		}
		return nil, errors.New("algebra: cannot solve this equation exactly")
	}
	return solveRatPoly(rc, v)
}

// solveRatPoly solves a polynomial with rational coefficients.
func solveRatPoly(rc []*big.Rat, v Expr) ([]Expr, error) {
	rc = trimRat(rc)
	var roots []Expr
	add := func(r Expr) {
		for _, ex := range roots {
			if ex.Equal(r) {
				return
			}
		}
		roots = append(roots, r)
	}
	// Extract rational roots exactly.
	core := append([]*big.Rat(nil), rc...)
	for {
		rr := rationalRoots(core)
		if len(rr) == 0 {
			break
		}
		for _, r := range rr {
			add(newRational(new(big.Rat).Set(r)))
			for ratPolyEval(core, r).Sign() == 0 && ratDegree(core) >= 1 {
				core = deflate(core, r)
			}
		}
	}
	switch ratDegree(core) {
	case 0:
		// Fully solved by rational roots (or a nonzero constant remains).
	case 1:
		add(newRational(new(big.Rat).Neg(new(big.Rat).Quo(core[0], core[1]))))
	case 2:
		for _, r := range solveQuadratic(
			newRational(new(big.Rat).Set(core[2])),
			newRational(new(big.Rat).Set(core[1])),
			newRational(new(big.Rat).Set(core[0])), v) {
			add(r)
		}
	default:
		for _, r := range durandKerner(core) {
			add(complexToExpr(r))
		}
	}
	if len(roots) == 0 {
		return nil, errors.New("algebra: no solution")
	}
	sortExprs(roots)
	return roots, nil
}

// solveQuadratic solves a*x^2 + b*x + c = 0 using the quadratic formula,
// returning complex conjugate roots a±b*I when the discriminant is negative.
func solveQuadratic(a, b, c, v Expr) []Expr {
	disc := Simplify(Add(Pow(b, Int(2)), Mul(Int(-4), a, c)))
	var root Expr
	if isNum(disc) && numSign(disc) < 0 {
		root = Mul(I, sqrtExpr(neg(disc)))
	} else {
		root = sqrtExpr(disc)
	}
	denom := Mul(Int(2), a)
	r1 := Simplify(Expand(Mul(Add(neg(b), root), Pow(denom, Int(-1)))))
	r2 := Simplify(Expand(Mul(Add(neg(b), neg(root)), Pow(denom, Int(-1)))))
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

// durandKerner returns numeric approximations to all roots (real and complex)
// of the polynomial with rational coefficients c using the Durand–Kerner
// iteration.
func durandKerner(c []*big.Rat) []complex128 {
	c = trimRat(c)
	n := ratDegree(c)
	if n < 1 {
		return nil
	}
	// Monic complex coefficients, index i = x^i.
	coef := make([]complex128, n+1)
	lead, _ := c[n].Float64()
	for i := 0; i <= n; i++ {
		f, _ := c[i].Float64()
		coef[i] = complex(f/lead, 0)
	}
	eval := func(z complex128) complex128 {
		acc := complex(0, 0)
		for i := n; i >= 0; i-- {
			acc = acc*z + coef[i]
		}
		return acc
	}
	// Distinct initial guesses on a spiral.
	roots := make([]complex128, n)
	seed := complex(0.4, 0.9)
	cur := complex(1, 0)
	for i := 0; i < n; i++ {
		cur *= seed
		roots[i] = cur
	}
	for iter := 0; iter < 500; iter++ {
		maxDelta := 0.0
		for i := 0; i < n; i++ {
			denom := complex(1, 0)
			for j := 0; j < n; j++ {
				if j != i {
					denom *= roots[i] - roots[j]
				}
			}
			if denom == 0 {
				continue
			}
			delta := eval(roots[i]) / denom
			roots[i] -= delta
			if d := cmplx.Abs(delta); d > maxDelta {
				maxDelta = d
			}
		}
		if maxDelta < 1e-14 {
			break
		}
	}
	return roots
}

// complexToExpr converts a numeric root to an expression, snapping values very
// close to integers, and dropping negligible real or imaginary parts.
func complexToExpr(z complex128) Expr {
	re := snap(real(z))
	im := snap(imag(z))
	if im == 0 {
		return realToExpr(re)
	}
	return Add(realToExpr(re), Mul(realToExpr(im), I))
}

func realToExpr(x float64) Expr {
	if x == math.Trunc(x) && math.Abs(x) < 1e15 {
		return Int(int64(x))
	}
	return Flt(x)
}

// snap rounds x to the nearest integer when it is within a small tolerance.
func snap(x float64) float64 {
	r := math.Round(x)
	if math.Abs(x-r) < 1e-9 {
		return r
	}
	if math.Abs(x) < 1e-12 {
		return 0
	}
	return x
}

// SolveSystem solves a system of linear equations eqs (each read as eq == 0)
// for the given syms, using Gaussian elimination over the rationals. It returns
// the solution values aligned with syms. The number of equations must equal the
// number of unknowns and the system must be linear with a unique solution;
// otherwise an error is returned. Each entry of syms must be a [Symbol].
func SolveSystem(eqs []Expr, syms []Expr) ([]Expr, error) {
	n := len(syms)
	if n == 0 {
		return nil, errors.New("algebra: no unknowns")
	}
	if len(eqs) != n {
		return nil, errors.New("algebra: need as many equations as unknowns")
	}
	names := make([]string, n)
	for j, sy := range syms {
		s, ok := sy.(*Symbol)
		if !ok {
			return nil, errors.New("algebra: SolveSystem requires symbols")
		}
		names[j] = s.Name
	}
	// Build augmented matrix A|b with A x = b, where b = -constant term.
	a := make([][]*big.Rat, n)
	b := make([]*big.Rat, n)
	for i, eq := range eqs {
		a[i] = make([]*big.Rat, n)
		linComb := Expr(Int(0))
		for j := range names {
			cij := Simplify(Diff(eq, syms[j]))
			for _, nm := range names {
				if containsSym(cij, nm) {
					return nil, errors.New("algebra: system is not linear")
				}
			}
			r, ok := toRat(cij)
			if !ok {
				return nil, errors.New("algebra: non-rational coefficient")
			}
			a[i][j] = new(big.Rat).Set(r)
			linComb = Add(linComb, Mul(cij, syms[j]))
		}
		// Constant term: equation with all unknowns set to zero.
		constExpr := eq
		for _, sy := range syms {
			constExpr = Subs(constExpr, sy, Int(0))
		}
		constExpr = Simplify(constExpr)
		// Verify linearity: eq == linComb + const.
		if !Simplify(Expand(Add(eq, neg(Add(linComb, constExpr))))).Equal(Int(0)) {
			return nil, errors.New("algebra: system is not linear")
		}
		r, ok := toRat(constExpr)
		if !ok {
			return nil, errors.New("algebra: non-rational constant")
		}
		b[i] = new(big.Rat).Neg(r)
	}
	sol, err := gaussSolve(a, b)
	if err != nil {
		return nil, err
	}
	out := make([]Expr, n)
	for j := range sol {
		out[j] = newRational(sol[j])
	}
	return out, nil
}

// gaussSolve solves the n×n rational linear system a x = b by Gaussian
// elimination with partial pivoting.
func gaussSolve(a [][]*big.Rat, b []*big.Rat) ([]*big.Rat, error) {
	n := len(a)
	m := make([][]*big.Rat, n)
	for i := range a {
		m[i] = make([]*big.Rat, n+1)
		for j := 0; j < n; j++ {
			m[i][j] = new(big.Rat).Set(a[i][j])
		}
		m[i][n] = new(big.Rat).Set(b[i])
	}
	for col := 0; col < n; col++ {
		piv := -1
		for r := col; r < n; r++ {
			if m[r][col].Sign() != 0 {
				piv = r
				break
			}
		}
		if piv < 0 {
			return nil, errors.New("algebra: singular or underdetermined system")
		}
		m[col], m[piv] = m[piv], m[col]
		inv := new(big.Rat).Inv(m[col][col])
		for j := col; j <= n; j++ {
			m[col][j].Mul(m[col][j], inv)
		}
		for r := 0; r < n; r++ {
			if r == col || m[r][col].Sign() == 0 {
				continue
			}
			factor := new(big.Rat).Set(m[r][col])
			for j := col; j <= n; j++ {
				m[r][j].Sub(m[r][j], new(big.Rat).Mul(factor, m[col][j]))
			}
		}
	}
	sol := make([]*big.Rat, n)
	for i := 0; i < n; i++ {
		sol[i] = new(big.Rat).Set(m[i][n])
	}
	return sol, nil
}
