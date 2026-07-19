package ellipticcurves

import (
	"errors"
	"math/big"
)

// ErrBadReduction indicates that a rational curve cannot be reduced to a smooth
// curve at the given prime, either because the prime divides a coefficient
// denominator or because the reduced curve is singular (bad reduction).
var ErrBadReduction = errors.New("ellipticcurves: bad reduction at prime")

// DiscriminantInt returns the discriminant of an integral curve as an integer,
// together with true. For a non-integral curve it returns (nil, false).
func (c *CurveQ) DiscriminantInt() (*big.Int, bool) {
	if !c.IsIntegralCurve() {
		return nil, false
	}
	d := c.Discriminant()
	if !d.IsInt() {
		return nil, false
	}
	return new(big.Int).Set(d.Num()), true
}

// BadPrimes returns the primes of bad reduction of an integral curve, namely the
// prime divisors of its discriminant. It returns ErrNonIntegralCurve for a
// non-integral curve.
func (c *CurveQ) BadPrimes() ([]*big.Int, error) {
	d, ok := c.DiscriminantInt()
	if !ok {
		return nil, ErrNonIntegralCurve
	}
	return PrimeFactors(d), nil
}

// HasGoodReduction reports whether an integral curve has good reduction at the
// prime p, i.e. p does not divide the discriminant. It returns false when p is
// not prime or the curve is non-integral.
func (c *CurveQ) HasGoodReduction(p *big.Int) bool {
	if !p.ProbablyPrime(20) {
		return false
	}
	d, ok := c.DiscriminantInt()
	if !ok {
		return false
	}
	return new(big.Int).Mod(d, p).Sign() != 0
}

// ratModP reduces a rational number modulo the prime p, returning nil when p
// divides its denominator.
func ratModP(r *big.Rat, p *big.Int) *big.Int {
	den := new(big.Int).Mod(r.Denom(), p)
	if den.Sign() == 0 {
		return nil
	}
	inv := new(big.Int).ModInverse(den, p)
	if inv == nil {
		return nil
	}
	num := new(big.Int).Mod(r.Num(), p)
	return ModMul(num, inv, p)
}

// ReduceModP returns the reduction of the rational curve modulo the prime p as a
// CurveFp. It returns ErrBadReduction when p divides a coefficient denominator
// or when the reduced curve is singular.
func (c *CurveQ) ReduceModP(p *big.Int) (*CurveFp, error) {
	if p.Cmp(bigThree) <= 0 || !p.ProbablyPrime(20) {
		return nil, ErrNotPrime
	}
	a := ratModP(c.A, p)
	b := ratModP(c.B, p)
	if a == nil || b == nil {
		return nil, ErrBadReduction
	}
	curve, err := NewCurveFp(a, b, p)
	if err != nil {
		return nil, ErrBadReduction
	}
	return curve, nil
}

// ReducePointModP reduces a rational point of the curve modulo the prime p to a
// point on the reduced curve, returning ErrBadReduction when a coordinate has a
// denominator divisible by p. The point at infinity maps to the point at
// infinity.
func (c *CurveQ) ReducePointModP(p *big.Int, pt PointQ) (PointFp, error) {
	if pt.Infinity {
		return PointAtInfinityFp(), nil
	}
	x := ratModP(pt.X, p)
	y := ratModP(pt.Y, p)
	if x == nil || y == nil {
		return PointFp{}, ErrBadReduction
	}
	return PointFp{X: x, Y: y}, nil
}

// APAtPrime returns the trace of Frobenius a_p = p + 1 - #E~(F_p) of the good
// reduction of an integral curve at the prime p. It returns ErrBadReduction at a
// prime of bad reduction.
func (c *CurveQ) APAtPrime(p *big.Int) (*big.Int, error) {
	reduced, err := c.ReduceModP(p)
	if err != nil {
		return nil, err
	}
	return reduced.TraceOfFrobenius(), nil
}
