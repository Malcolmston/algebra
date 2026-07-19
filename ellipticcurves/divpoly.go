package ellipticcurves

import (
	"errors"
	"math/big"
)

// ErrDivisionPolyAtTwoTorsion indicates that a division-polynomial value was
// requested at a point with y = 0, where the numeric recursion (which divides
// by 2y for even indices) does not apply.
var ErrDivisionPolyAtTwoTorsion = errors.New("ellipticcurves: division polynomial undefined at y=0 for even index")

// DivisionPolynomial2 returns psi_2(x, y) = 2*y evaluated modulo p.
func (c *CurveFp) DivisionPolynomial2(x, y *big.Int) *big.Int {
	return ModDouble(y, c.P)
}

// DivisionPolynomial3 returns psi_3(x, y) = 3*x^4 + 6*A*x^2 + 12*B*x - A^2
// evaluated modulo p. It depends only on x.
func (c *CurveFp) DivisionPolynomial3(x, y *big.Int) *big.Int {
	p := c.P
	x2 := ModSquare(x, p)
	x4 := ModSquare(x2, p)
	t := ModMul(bigThree, x4, p)
	t = ModAdd(t, ModMul(big.NewInt(6), ModMul(c.A, x2, p), p), p)
	t = ModAdd(t, ModMul(big.NewInt(12), ModMul(c.B, x, p), p), p)
	t = ModSub(t, ModSquare(c.A, p), p)
	return t
}

// DivisionPolynomial4 returns psi_4(x, y) modulo p, equal to
// 4*y*(x^6 + 5*A*x^4 + 20*B*x^3 - 5*A^2*x^2 - 4*A*B*x - 8*B^2 - A^3).
func (c *CurveFp) DivisionPolynomial4(x, y *big.Int) *big.Int {
	p := c.P
	x2 := ModSquare(x, p)
	x3 := ModMul(x2, x, p)
	x4 := ModSquare(x2, p)
	x6 := ModMul(x4, x2, p)
	a2 := ModSquare(c.A, p)
	inner := x6
	inner = ModAdd(inner, ModMul(big.NewInt(5), ModMul(c.A, x4, p), p), p)
	inner = ModAdd(inner, ModMul(big.NewInt(20), ModMul(c.B, x3, p), p), p)
	inner = ModSub(inner, ModMul(big.NewInt(5), ModMul(a2, x2, p), p), p)
	inner = ModSub(inner, ModMul(bigFour, ModMul(ModMul(c.A, c.B, p), x, p), p), p)
	inner = ModSub(inner, ModMul(big.NewInt(8), ModSquare(c.B, p), p), p)
	inner = ModSub(inner, ModMul(c.A, a2, p), p)
	return ModMul(ModMul(bigFour, y, p), inner, p)
}

// DivisionPolynomialValue returns psi_n(x, y) evaluated at the point (x, y)
// modulo p, using the standard recursion. The point must lie on the curve. For
// even n the recursion divides by 2*y, so a point with y = 0 yields
// ErrDivisionPolyAtTwoTorsion. Negative n uses psi_{-n} = -psi_n.
func (c *CurveFp) DivisionPolynomialValue(n int, x, y *big.Int) (*big.Int, error) {
	if n < 0 {
		v, err := c.DivisionPolynomialValue(-n, x, y)
		if err != nil {
			return nil, err
		}
		return ModNeg(v, c.P), nil
	}
	p := c.P
	psi := make([]*big.Int, n+5)
	psi[0] = big.NewInt(0)
	if n+5 > 1 {
		psi[1] = big.NewInt(1)
	}
	if n+5 > 2 {
		psi[2] = c.DivisionPolynomial2(x, y)
	}
	if n+5 > 3 {
		psi[3] = c.DivisionPolynomial3(x, y)
	}
	if n+5 > 4 {
		psi[4] = c.DivisionPolynomial4(x, y)
	}
	if n <= 4 {
		return new(big.Int).Set(psi[n]), nil
	}
	twoY := ModDouble(y, p)
	var twoYinv *big.Int
	if twoY.Sign() != 0 {
		twoYinv = new(big.Int).ModInverse(twoY, p)
	}
	for k := 5; k <= n; k++ {
		if k%2 == 1 {
			m := (k - 1) / 2
			// psi_{2m+1} = psi_{m+2} psi_m^3 - psi_{m-1} psi_{m+1}^3
			a := ModMul(psi[m+2], powMod(psi[m], 3, p), p)
			b := ModMul(psi[m-1], powMod(psi[m+1], 3, p), p)
			psi[k] = ModSub(a, b, p)
		} else {
			m := k / 2
			// psi_{2m} = psi_m (psi_{m+2} psi_{m-1}^2 - psi_{m-2} psi_{m+1}^2) / (2y)
			if twoYinv == nil {
				return nil, ErrDivisionPolyAtTwoTorsion
			}
			t1 := ModMul(psi[m+2], ModSquare(psi[m-1], p), p)
			t2 := ModMul(psi[m-2], ModSquare(psi[m+1], p), p)
			inner := ModSub(t1, t2, p)
			num := ModMul(psi[m], inner, p)
			psi[k] = ModMul(num, twoYinv, p)
		}
	}
	return new(big.Int).Set(psi[n]), nil
}

// powMod returns v^e mod p for a small non-negative exponent.
func powMod(v *big.Int, e int, p *big.Int) *big.Int {
	return new(big.Int).Exp(v, big.NewInt(int64(e)), p)
}

// IsNTorsionByDivisionPoly reports whether pt is an n-torsion point (n*pt = O)
// by testing whether psi_n vanishes at pt. For a 2-torsion point it falls back
// to the direct group-law test. The point at infinity is n-torsion for every n.
func (c *CurveFp) IsNTorsionByDivisionPoly(n int, pt PointFp) (bool, error) {
	if pt.Infinity {
		return true, nil
	}
	if pt.Y.Sign() == 0 {
		// 2-torsion: n*pt = O iff n even.
		return n%2 == 0, nil
	}
	v, err := c.DivisionPolynomialValue(n, pt.X, pt.Y)
	if err != nil {
		return false, err
	}
	return v.Sign() == 0, nil
}

// IsNTorsionPoint reports whether n*pt equals the point at infinity, computed
// directly with the group law. It is the reference against which the
// division-polynomial test is checked.
func (c *CurveFp) IsNTorsionPoint(n *big.Int, pt PointFp) bool {
	return c.ScalarMul(n, pt).Infinity
}

// DivisionPolynomialSequence returns psi_0, psi_1, ..., psi_n evaluated at
// (x, y) modulo p. It reuses the recursion of DivisionPolynomialValue.
func (c *CurveFp) DivisionPolynomialSequence(n int, x, y *big.Int) ([]*big.Int, error) {
	out := make([]*big.Int, n+1)
	for k := 0; k <= n; k++ {
		v, err := c.DivisionPolynomialValue(k, x, y)
		if err != nil {
			return nil, err
		}
		out[k] = v
	}
	return out, nil
}
