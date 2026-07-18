package ecc

import (
	"errors"
	"io"
	"math/big"
)

// DiscriminantFp returns the discriminant -16*(4*a^3 + 27*b^2) of the curve
// y^2 = x^3 + a*x + b modulo p, without constructing a CurveFp. The modulus p
// must be positive.
func DiscriminantFp(a, b, p *big.Int) *big.Int {
	a3 := ModMul(ModMul(a, a, p), a, p)
	b2 := ModMul(b, b, p)
	term := ModAdd(ModMul(big.NewInt(4), a3, p), ModMul(big.NewInt(27), b2, p), p)
	return ModMul(big.NewInt(-16), term, p)
}

// JInvariantFp returns the j-invariant 1728 * 4*a^3 / (4*a^3 + 27*b^2) of the
// curve y^2 = x^3 + a*x + b modulo p, without constructing a CurveFp, together
// with a boolean that is false when the curve is singular (its discriminant is
// zero, making the j-invariant undefined).
func JInvariantFp(a, b, p *big.Int) (*big.Int, bool) {
	a3 := ModMul(ModMul(a, a, p), a, p)
	b2 := ModMul(b, b, p)
	fourA3 := ModMul(big.NewInt(4), a3, p)
	denom := ModAdd(fourA3, ModMul(big.NewInt(27), b2, p), p)
	num := ModMul(big.NewInt(1728), fourA3, p)
	return ModDiv(num, denom, p)
}

// SameJInvariant reports whether this curve and other share the same
// j-invariant over the same base field, a necessary condition for the two
// curves to be isomorphic over the algebraic closure of GF(p). It returns false
// if the curves have different moduli.
func (c *CurveFp) SameJInvariant(other *CurveFp) bool {
	if c.P.Cmp(other.P) != 0 {
		return false
	}
	return c.JInvariant().Cmp(other.JInvariant()) == 0
}

// OnCurveX reports whether the field element x is the abscissa of some point on
// the curve, i.e. whether x^3 + A*x + B is a square modulo P.
func (c *CurveFp) OnCurveX(x *big.Int) bool {
	return Legendre(c.RHS(x), c.P) >= 0
}

// Cofactor returns the cofactor #E(GF(p)) / n, the index of the order-n subgroup
// in the full group, together with a boolean that is false when n does not
// divide the group order. The full order is obtained by naive point counting,
// so this is intended for small primes.
func (c *CurveFp) Cofactor(n *big.Int) (*big.Int, bool) {
	order := c.CountPointsNaive()
	q, r := new(big.Int).QuoRem(order, n, new(big.Int))
	if r.Sign() != 0 {
		return nil, false
	}
	return q, true
}

// IsGenerator reports whether pt has order exactly n in the curve group, i.e.
// whether n*pt is the identity and pt generates a subgroup of size n. It uses
// the baby-step/giant-step order routine.
func (c *CurveFp) IsGenerator(pt PointFp, n *big.Int) bool {
	if pt.Infinity {
		return n.Cmp(big.NewInt(1)) == 0
	}
	ord := c.PointOrderBSGS(pt)
	return ord != nil && ord.Cmp(n) == 0
}

// eccRandInt returns a uniformly distributed integer in [0, max) drawn from the
// supplied reader by rejection sampling. max must be positive.
func eccRandInt(reader io.Reader, max *big.Int) (*big.Int, error) {
	if max.Sign() <= 0 {
		return nil, errors.New("ecc: random bound must be positive")
	}
	bitLen := max.BitLen()
	byteLen := (bitLen + 7) / 8
	// Number of high bits to mask off so the candidate stays close to max.
	excess := uint(byteLen*8 - bitLen)
	buf := make([]byte, byteLen)
	for {
		if _, err := io.ReadFull(reader, buf); err != nil {
			return nil, err
		}
		if excess > 0 {
			buf[0] &= 0xff >> excess
		}
		candidate := new(big.Int).SetBytes(buf)
		if candidate.Cmp(max) < 0 {
			return candidate, nil
		}
	}
}

// GenerateKeyPair draws a private scalar uniformly from [1, n) using the
// supplied reader and returns it together with the public point priv*G. It is
// deterministic given a deterministic reader, which makes it testable; callers
// wanting real keys should pass crypto/rand.Reader. It returns an error if the
// reader fails.
func (c *CurveFp) GenerateKeyPair(g PointFp, n *big.Int, reader io.Reader) (priv *big.Int, pub PointFp, err error) {
	if n.Cmp(big.NewInt(2)) < 0 {
		return nil, PointFp{}, errors.New("ecc: GenerateKeyPair requires n >= 2")
	}
	// Sample from [0, n-1) then shift into [1, n).
	nm1 := new(big.Int).Sub(n, big.NewInt(1))
	r, err := eccRandInt(reader, nm1)
	if err != nil {
		return nil, PointFp{}, err
	}
	priv = r.Add(r, big.NewInt(1))
	pub = c.ScalarMul(priv, g)
	return priv, pub, nil
}

// GenerateKeyPair draws a private key for the named curve from the supplied
// reader and returns it together with the corresponding public point. It is a
// convenience wrapper over CurveFp.GenerateKeyPair using the curve's standard
// base point and subgroup order.
func (nc NamedCurve) GenerateKeyPair(reader io.Reader) (priv *big.Int, pub PointFp, err error) {
	return nc.Curve.GenerateKeyPair(nc.G, nc.N, reader)
}
