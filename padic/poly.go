package padic

import "math/big"

// PolyTrim removes trailing (high-degree) zero coefficients, returning the
// polynomial in canonical low-to-high form. A polynomial equal to zero becomes
// the empty slice.
func PolyTrim(coeffs []*big.Int) []*big.Int {
	n := len(coeffs)
	for n > 0 && coeffs[n-1].Sign() == 0 {
		n--
	}
	out := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		out[i] = new(big.Int).Set(coeffs[i])
	}
	return out
}

// PolyDegree returns the degree of the integer polynomial coeffs (low-to-high),
// or -1 for the zero polynomial.
func PolyDegree(coeffs []*big.Int) int {
	for i := len(coeffs) - 1; i >= 0; i-- {
		if coeffs[i].Sign() != 0 {
			return i
		}
	}
	return -1
}

// PolyIsZero reports whether every coefficient is zero.
func PolyIsZero(coeffs []*big.Int) bool { return PolyDegree(coeffs) == -1 }

// PolyAdd returns the sum of two integer polynomials given low-to-high.
func PolyAdd(a, b []*big.Int) []*big.Int {
	n := maxInt(len(a), len(b))
	out := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		s := new(big.Int)
		if i < len(a) {
			s.Add(s, a[i])
		}
		if i < len(b) {
			s.Add(s, b[i])
		}
		out[i] = s
	}
	return PolyTrim(out)
}

// PolyNeg returns the negation of an integer polynomial.
func PolyNeg(a []*big.Int) []*big.Int {
	out := make([]*big.Int, len(a))
	for i := range a {
		out[i] = new(big.Int).Neg(a[i])
	}
	return out
}

// PolySub returns a - b for integer polynomials given low-to-high.
func PolySub(a, b []*big.Int) []*big.Int { return PolyAdd(a, PolyNeg(b)) }

// PolyScale multiplies every coefficient of a by the scalar c.
func PolyScale(a []*big.Int, c *big.Int) []*big.Int {
	out := make([]*big.Int, len(a))
	for i := range a {
		out[i] = new(big.Int).Mul(a[i], c)
	}
	return PolyTrim(out)
}

// PolyMul returns the product of two integer polynomials given low-to-high.
func PolyMul(a, b []*big.Int) []*big.Int {
	if PolyIsZero(a) || PolyIsZero(b) {
		return []*big.Int{}
	}
	out := make([]*big.Int, len(a)+len(b)-1)
	for i := range out {
		out[i] = new(big.Int)
	}
	for i := range a {
		if a[i].Sign() == 0 {
			continue
		}
		for j := range b {
			out[i+j].Add(out[i+j], new(big.Int).Mul(a[i], b[j]))
		}
	}
	return PolyTrim(out)
}

// PolyModPk reduces every coefficient of a modulo p^k into [0, p^k).
func PolyModPk(a []*big.Int, p *big.Int, k int) []*big.Int {
	mod := PPow(p, k)
	out := make([]*big.Int, len(a))
	for i := range a {
		out[i] = new(big.Int).Mod(a[i], mod)
	}
	return out
}

// PolyValuation returns the minimum p-adic valuation over the non-zero
// coefficients of a, and false if a is the zero polynomial.
func PolyValuation(p *big.Int, a []*big.Int) (int, bool) {
	best := 0
	found := false
	for _, c := range a {
		if c.Sign() == 0 {
			continue
		}
		v := ValuationInt(p, c)
		if !found || v < best {
			best = v
			found = true
		}
	}
	return best, found
}

// PolyContent returns the p-adic content of a: p raised to the minimum
// coefficient valuation, i.e. the largest power of p dividing every
// coefficient. The zero polynomial yields nil.
func PolyContent(p *big.Int, a []*big.Int) *big.Int {
	v, ok := PolyValuation(p, a)
	if !ok {
		return nil
	}
	return PPow(p, v)
}

// PolyEvalRat evaluates an integer polynomial at a rational point x = n/d.
func PolyEvalRat(coeffs []*big.Int, x *big.Rat) *big.Rat {
	acc := new(big.Rat)
	for i := len(coeffs) - 1; i >= 0; i-- {
		acc.Mul(acc, x)
		acc.Add(acc, new(big.Rat).SetInt(coeffs[i]))
	}
	return acc
}

// IsEisenstein reports whether the integer polynomial coeffs (low-to-high) is
// Eisenstein at p: the leading coefficient is a unit, every lower coefficient
// is divisible by p, and the constant term is divisible by p but not p^2. Such
// a polynomial is irreducible over Q_p and generates a totally ramified
// extension.
func IsEisenstein(p *big.Int, coeffs []*big.Int) bool {
	deg := PolyDegree(coeffs)
	if deg < 1 {
		return false
	}
	if ValuationInt(p, coeffs[deg]) != 0 {
		return false
	}
	for i := 0; i < deg; i++ {
		if ValuationInt(p, coeffs[i]) < 1 {
			return false
		}
	}
	return ValuationInt(p, coeffs[0]) == 1
}
