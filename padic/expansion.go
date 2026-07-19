package padic

import "math/big"

// ExpandRational returns the p-adic expansion of the rational a/b: the
// valuation v and the first count digits d[0], d[1], ... of the unit, each in
// [0, p), such that a/b = p^v * (d[0] + d[1]*p + d[2]*p^2 + ...). b must be
// non-zero and count positive. A zero numerator yields v = count and all-zero
// digits.
func ExpandRational(p *big.Int, a, b *big.Int, count int) (int, []*big.Int, error) {
	if b.Sign() == 0 {
		return 0, nil, ErrZeroDivision
	}
	if count <= 0 {
		return 0, nil, ErrPrecision
	}
	if a.Sign() == 0 {
		zeros := make([]*big.Int, count)
		for i := range zeros {
			zeros[i] = big.NewInt(0)
		}
		return count, zeros, nil
	}
	v := ValuationRat(p, a, b)
	num := UnitPartInt(p, a)
	den := UnitPartInt(p, b)
	digits := digitsOfUnit(p, num, den, count)
	return v, digits, nil
}

// digitsOfUnit returns the first count p-adic digits of the unit num/den, where
// den is coprime to p. Each digit is in [0, p).
func digitsOfUnit(p, num, den *big.Int, count int) []*big.Int {
	out := make([]*big.Int, count)
	n := new(big.Int).Set(num)
	pinv := new(big.Int).ModInverse(new(big.Int).Mod(den, p), p)
	for i := 0; i < count; i++ {
		nm := new(big.Int).Mod(n, p)
		d := new(big.Int).Mul(nm, pinv)
		d.Mod(d, p)
		out[i] = d
		// n <- (n - d*den) / p
		n.Sub(n, new(big.Int).Mul(d, den))
		n.Div(n, p)
	}
	return out
}

// Digits returns the p-adic digits of the unit part of x: a slice of length
// RelativePrecision, each entry in [0, p), least significant first. For the
// zero element it returns an empty slice.
func (x *Padic) Digits() ([]*big.Int, error) {
	if x.IsZero() {
		return []*big.Int{}, nil
	}
	return digitsOfUnit(x.p, x.unit, bigOne, x.prec), nil
}

// FullDigits returns the p-adic digits of x from the constant term up to (but
// not including) the absolute precision. Positions below the valuation are
// zero. The slice has length AbsolutePrecision (clamped at 0).
func (x *Padic) FullDigits() ([]*big.Int, error) {
	n := x.AbsolutePrecision()
	if n <= 0 || x.IsZero() {
		out := make([]*big.Int, maxInt(n, 0))
		for i := range out {
			out[i] = big.NewInt(0)
		}
		return out, nil
	}
	out := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		out[i] = big.NewInt(0)
	}
	unitDigits := digitsOfUnit(x.p, x.unit, bigOne, x.prec)
	for i, d := range unitDigits {
		pos := x.val + i
		if pos >= 0 && pos < n {
			out[pos] = d
		}
	}
	return out, nil
}

// DigitsToRat reconstructs the rational p^val * sum(digits[i]*p^i) from a
// valuation and a slice of digits. It is the inverse of ExpandRational up to
// the finite number of digits provided.
func DigitsToRat(p *big.Int, val int, digits []*big.Int) *big.Rat {
	sum := new(big.Int)
	pw := big.NewInt(1)
	for _, d := range digits {
		sum.Add(sum, new(big.Int).Mul(d, pw))
		pw.Mul(pw, p)
	}
	r := new(big.Rat).SetInt(sum)
	if val >= 0 {
		r.Mul(r, new(big.Rat).SetInt(PPow(p, val)))
	} else {
		r.Quo(r, new(big.Rat).SetInt(PPow(p, -val)))
	}
	return r
}

// DigitsToPadic builds a p-adic number from a valuation and a slice of digits
// (least significant first), to relative precision len(digits).
func DigitsToPadic(p *big.Int, val int, digits []*big.Int) *Padic {
	if len(digits) == 0 {
		return newZero(p, val)
	}
	sum := new(big.Int)
	pw := big.NewInt(1)
	for _, d := range digits {
		sum.Add(sum, new(big.Int).Mul(d, pw))
		pw.Mul(pw, p)
	}
	return makeScaled(p, val, sum, val+len(digits))
}

// RationalReconstruction attempts to recover a rational a/b congruent to x
// modulo p^absprec, with |a|, |b| bounded by sqrt(p^absprec / 2). It returns
// numerator and denominator, or an error if no such rational is found. This is
// the standard p-adic (Wang) rational reconstruction.
func (x *Padic) RationalReconstruction() (*big.Int, *big.Int, error) {
	m := PPow(x.p, x.AbsolutePrecision())
	var r *big.Int
	if bi, err := x.BigInt(); err == nil {
		r = new(big.Int).Mod(bi, m)
	} else {
		// Non-integral: scale by p^(-val) using modular inverse of p.
		mod := PPow(x.p, x.prec)
		pinv := new(big.Int).ModInverse(new(big.Int).Exp(x.p, big.NewInt(int64(-x.val)), mod), mod)
		r = new(big.Int).Mul(x.unit, pinv)
		r.Mod(r, mod)
		m = mod
	}
	bound := new(big.Int).Sqrt(new(big.Int).Div(m, bigTwo))
	num, den, ok := rationalReconstruct(r, m, bound)
	if !ok {
		return nil, nil, ErrNoRoot
	}
	return num, den, nil
}

// rationalReconstruct runs the extended-Euclidean rational reconstruction of r
// modulo m with numerator/denominator bounded by bound.
func rationalReconstruct(r, m, bound *big.Int) (*big.Int, *big.Int, bool) {
	r0 := new(big.Int).Set(m)
	r1 := new(big.Int).Mod(r, m)
	s0 := big.NewInt(0)
	s1 := big.NewInt(1)
	for r1.Sign() != 0 && r1.Cmp(bound) > 0 {
		q := new(big.Int).Div(r0, r1)
		r0, r1 = r1, new(big.Int).Sub(r0, new(big.Int).Mul(q, r1))
		s0, s1 = s1, new(big.Int).Sub(s0, new(big.Int).Mul(q, s1))
	}
	if s1.Sign() == 0 {
		return nil, nil, false
	}
	if new(big.Int).Abs(s1).Cmp(bound) > 0 {
		return nil, nil, false
	}
	num := new(big.Int).Set(r1)
	den := new(big.Int).Set(s1)
	if den.Sign() < 0 {
		num.Neg(num)
		den.Neg(den)
	}
	if !Coprime(num, den) {
		g := GCD(num, den)
		num.Div(num, g)
		den.Div(den, g)
	}
	return num, den, true
}
