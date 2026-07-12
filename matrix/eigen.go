package matrix

import (
	"math"
	"math/big"
	"sort"

	"github.com/malcolmston/algebra"
)

// Lambda is the symbol used as the eigenvalue variable in [Matrix.CharPoly].
var Lambda = algebra.Sym("lambda")

// CharPoly returns the characteristic polynomial det(A - λI) as an algebra.Expr
// in the symbol [Lambda]. It is defined for any square matrix (numeric or
// symbolic) and returns [ErrNotSquare] otherwise. The result is expanded and
// simplified, so it reads as a polynomial in lambda.
//
// The eigenvalues of A are the roots of this polynomial. Note the sign
// convention: for an n×n matrix the leading term is (-1)^n·λ^n; the roots are
// unaffected.
func (m *Matrix) CharPoly() (algebra.Expr, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	b := m.Clone()
	for i := 0; i < m.rows; i++ {
		b.data[i][i] = algebra.Add(m.data[i][i], algebra.Mul(algebra.Int(-1), Lambda))
	}
	d, err := b.Det()
	if err != nil {
		return nil, err
	}
	return simp(algebra.Expand(d)), nil
}

// Eigenvalues returns the eigenvalues of a square matrix.
//
// Coverage:
//   - 1×1 and 2×2 matrices (numeric or symbolic) are solved exactly by handing
//     the characteristic polynomial to the parent package's Solve, so symbolic
//     entries yield closed-form eigenvalues.
//   - 3×3 numeric matrices are solved exactly when a rational eigenvalue exists:
//     it is peeled off and the remaining quadratic is solved exactly (its roots
//     may be irrational, expressed with a symbolic sqrt). If no rational
//     eigenvalue exists the real eigenvalues are computed numerically and
//     returned as algebra.Flt values; complex eigenvalues are omitted.
//
// Limits: a 3×3 matrix with free symbols, and any matrix of size 4×4 or larger,
// return [ErrUnsupported]. Callers can still obtain [Matrix.CharPoly] and solve
// it by other means. Repeated eigenvalues are returned once.
func (m *Matrix) Eigenvalues() ([]algebra.Expr, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.rows
	switch {
	case n == 0:
		return nil, nil
	case n == 1:
		return []algebra.Expr{simp(m.data[0][0])}, nil
	case n == 2:
		p, err := m.CharPoly()
		if err != nil {
			return nil, err
		}
		roots, err := algebra.Solve(p, Lambda)
		if err != nil {
			return nil, ErrUnsupported
		}
		return roots, nil
	case n == 3:
		return m.eig3()
	default:
		return nil, ErrUnsupported
	}
}

// eig3 solves the 3×3 numeric case.
func (m *Matrix) eig3() ([]algebra.Expr, error) {
	p, err := m.CharPoly()
	if err != nil {
		return nil, err
	}
	coeffs, ok := cubicCoeffs(p)
	if !ok {
		return nil, ErrUnsupported // symbolic entries: no numeric eigenvalue path
	}
	// Try to peel a rational root, then solve the remaining quadratic exactly.
	if r, found := rationalRoot(coeffs); found {
		b2, b1, b0 := deflate(coeffs, r)
		roots := []algebra.Expr{ratToExpr(r)}
		q, err := algebra.Solve(buildQuadratic(b2, b1, b0), Lambda)
		if err == nil {
			roots = append(roots, q...)
		}
		return dedupSort(roots), nil
	}
	// No rational root: numeric real roots of the cubic.
	fc := make([]float64, 4)
	for i := 0; i < 4; i++ {
		f, _ := coeffs[i].Float64()
		fc[i] = f
	}
	reals := realCubicRoots(fc[3], fc[2], fc[1], fc[0])
	out := make([]algebra.Expr, 0, len(reals))
	for _, x := range reals {
		out = append(out, algebra.Flt(x))
	}
	return dedupSort(out), nil
}

// cubicCoeffs extracts the exact rational coefficients c0..c3 (ascending powers
// of lambda) of the degree-3 polynomial p by interpolating it at lambda =
// 0,1,2,3 and solving the resulting Vandermonde system exactly. It reports false
// if any interpolated value is not numeric (i.e. p has free symbols).
func cubicCoeffs(p algebra.Expr) ([]*big.Rat, bool) {
	xs := []int64{0, 1, 2, 3}
	v := New(4, 4)
	rhs := make([]algebra.Expr, 4)
	for i, x := range xs {
		for j := 0; j < 4; j++ {
			v.data[i][j] = algebra.Int(ipow(x, j))
		}
		rhs[i] = simp(algebra.Subs(p, Lambda, algebra.Int(x)))
	}
	sol, err := Solve(v, &Vector{data: rhs})
	if err != nil {
		return nil, false
	}
	out := make([]*big.Rat, 4)
	for i := 0; i < 4; i++ {
		r, ok := exprToRat(sol.data[i])
		if !ok {
			return nil, false
		}
		out[i] = r
	}
	return out, true
}

// exprToRat returns the exact rational value of a numeric expression, using its
// canonical string form. It reports false for anything that is not a plain
// integer, rational or float.
func exprToRat(e algebra.Expr) (*big.Rat, bool) {
	r, ok := new(big.Rat).SetString(simp(e).String())
	return r, ok
}

// ratToExpr converts a rational to an exact algebra.Expr (Integer or Rational).
func ratToExpr(r *big.Rat) algebra.Expr {
	if r.IsInt() {
		return algebra.IntBig(r.Num())
	}
	num := algebra.IntBig(r.Num())
	den := algebra.IntBig(r.Denom())
	return simp(algebra.Mul(num, algebra.Pow(den, algebra.Int(-1))))
}

// buildQuadratic constructs the expression a2·λ² + a1·λ + a0.
func buildQuadratic(a2, a1, a0 *big.Rat) algebra.Expr {
	return algebra.Add(
		algebra.Mul(ratToExpr(a2), algebra.Pow(Lambda, algebra.Int(2))),
		algebra.Mul(ratToExpr(a1), Lambda),
		ratToExpr(a0),
	)
}

// rationalRoot searches for a rational root of the cubic with the given
// ascending rational coefficients using the rational-root theorem. A zero
// constant term yields the root 0 directly. It gives up (returns false) when the
// integer numerators are too large to factor cheaply.
func rationalRoot(c []*big.Rat) (*big.Rat, bool) {
	if c[0].Sign() == 0 {
		return new(big.Rat), true // root 0
	}
	a := clearDenoms(c) // integer coefficients, same roots
	a0 := new(big.Int).Abs(a[0])
	an := new(big.Int).Abs(a[3])
	const cap = 1 << 20
	if !a0.IsInt64() || !an.IsInt64() || a0.Int64() > cap || an.Int64() > cap {
		return nil, false
	}
	ps := divisorsInt64(a0.Int64())
	qs := divisorsInt64(an.Int64())
	for _, dp := range ps {
		for _, dq := range qs {
			for _, s := range []int64{1, -1} {
				cand := big.NewRat(s*dp, dq)
				if evalRatPoly(c, cand).Sign() == 0 {
					return cand, true
				}
			}
		}
	}
	return nil, false
}

// clearDenoms scales the rational coefficients by the least common multiple of
// their denominators, returning proportional integer coefficients.
func clearDenoms(c []*big.Rat) []*big.Int {
	l := big.NewInt(1)
	for _, r := range c {
		l = lcm(l, r.Denom())
	}
	out := make([]*big.Int, len(c))
	lr := new(big.Rat).SetInt(l)
	for i, r := range c {
		t := new(big.Rat).Mul(r, lr) // integer-valued
		out[i] = new(big.Int).Set(t.Num())
	}
	return out
}

// evalRatPoly evaluates the polynomial with ascending rational coefficients at x
// using Horner's method.
func evalRatPoly(c []*big.Rat, x *big.Rat) *big.Rat {
	acc := new(big.Rat)
	for i := len(c) - 1; i >= 0; i-- {
		acc.Mul(acc, x)
		acc.Add(acc, c[i])
	}
	return acc
}

// deflate divides the cubic c3·x³+c2·x²+c1·x+c0 by (x - r) via synthetic
// division, returning the quotient quadratic coefficients (leading first).
func deflate(c []*big.Rat, r *big.Rat) (b2, b1, b0 *big.Rat) {
	b2 = new(big.Rat).Set(c[3])
	b1 = new(big.Rat).Add(c[2], new(big.Rat).Mul(r, b2))
	b0 = new(big.Rat).Add(c[1], new(big.Rat).Mul(r, b1))
	return b2, b1, b0
}

// dedupSort removes structurally duplicate expressions and returns them in a
// deterministic order (by canonical string form).
func dedupSort(es []algebra.Expr) []algebra.Expr {
	out := make([]algebra.Expr, 0, len(es))
	for _, e := range es {
		dup := false
		for _, o := range out {
			if simp(e).Equal(simp(o)) {
				dup = true
				break
			}
		}
		if !dup {
			out = append(out, e)
		}
	}
	sort.SliceStable(out, func(i, j int) bool { return out[i].String() < out[j].String() })
	return out
}

// --- small integer helpers ---------------------------------------------------

// ipow returns base**exp for small non-negative exponents.
func ipow(base int64, exp int) int64 {
	r := int64(1)
	for i := 0; i < exp; i++ {
		r *= base
	}
	return r
}

// divisorsInt64 returns the positive divisors of n (n > 0).
func divisorsInt64(n int64) []int64 {
	var ds []int64
	for i := int64(1); i*i <= n; i++ {
		if n%i == 0 {
			ds = append(ds, i)
			if i != n/i {
				ds = append(ds, n/i)
			}
		}
	}
	return ds
}

// lcm returns the least common multiple of two positive big integers.
func lcm(a, b *big.Int) *big.Int {
	g := new(big.Int).GCD(nil, nil, a, b)
	if g.Sign() == 0 {
		return new(big.Int).Set(a)
	}
	l := new(big.Int).Div(a, g)
	l.Mul(l, b)
	return l.Abs(l)
}

// realCubicRoots returns the distinct real roots of a·x³+b·x²+c·x+d (a != 0),
// used only as the numeric fallback when no rational eigenvalue exists.
func realCubicRoots(a, b, c, d float64) []float64 {
	// Normalize and depress: x = t - b/(3a).
	b, c, d = b/a, c/a, d/a
	off := -b / 3
	p := c - b*b/3
	q := 2*b*b*b/27 - b*c/3 + d
	const eps = 1e-9
	disc := q*q/4 + p*p*p/27
	switch {
	case disc > eps:
		sq := math.Sqrt(disc)
		u := math.Cbrt(-q/2 + sq)
		v := math.Cbrt(-q/2 - sq)
		return []float64{u + v + off}
	case disc < -eps:
		// Three distinct real roots via the trigonometric method (p < 0 here).
		m := 2 * math.Sqrt(-p/3)
		arg := 3 * q / (p * m) // = (3q)/(2p) * sqrt(-3/p)
		if arg > 1 {
			arg = 1
		} else if arg < -1 {
			arg = -1
		}
		phi := math.Acos(arg)
		return []float64{
			m*math.Cos(phi/3) + off,
			m*math.Cos((phi+2*math.Pi)/3) + off,
			m*math.Cos((phi+4*math.Pi)/3) + off,
		}
	default:
		if math.Abs(p) < eps {
			return []float64{off} // triple root
		}
		// Double root plus a simple root of the depressed cubic.
		return []float64{3*q/p + off, -3*q/(2*p) + off}
	}
}
