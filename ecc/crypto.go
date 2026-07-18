package ecc

import (
	"errors"
	"math/big"
)

// PublicKeyFp derives the public key d*G from a private scalar d and a base
// point G on the curve. The scalar d must be non-zero; the returned point is
// the caller's public key.
func (c *CurveFp) PublicKeyFp(g PointFp, d *big.Int) PointFp {
	return c.ScalarMul(d, g)
}

// ECDHSharedSecret computes the Elliptic-Curve Diffie-Hellman shared secret
// d*Q, where d is the local private scalar and Q is the remote party's public
// point. It returns the x-coordinate of the shared point, the conventional
// shared secret value, together with the full shared point. The boolean result
// is false if the shared point is the point at infinity, which indicates an
// invalid or small-subgroup public key.
func (c *CurveFp) ECDHSharedSecret(d *big.Int, q PointFp) (secret *big.Int, shared PointFp, ok bool) {
	shared = c.ScalarMul(d, q)
	if shared.Infinity {
		return nil, shared, false
	}
	return new(big.Int).Set(shared.X), shared, true
}

// ECDSASignature holds the pair (R, S) produced by the ECDSA signing equation.
type ECDSASignature struct {
	// R is the first signature component, the x-coordinate of k*G reduced
	// modulo the subgroup order n.
	R *big.Int
	// S is the second signature component.
	S *big.Int
}

// ECDSASign produces an ECDSA signature (r, s) for the integer message digest e
// under the private key d, using the base point G of prime subgroup order n and
// the caller-supplied per-message secret k. Supplying k explicitly keeps the
// routine deterministic and testable; production callers must draw a fresh,
// uniformly random, secret k in [1, n) for every signature. It returns an error
// when k is out of range or when the resulting r or s would be zero, in which
// case the caller should retry with a different k.
func (c *CurveFp) ECDSASign(g PointFp, n, d, e, k *big.Int) (ECDSASignature, error) {
	if n.Sign() <= 0 {
		return ECDSASignature{}, errors.New("ecc: ECDSASign requires a positive order n")
	}
	km := eccMod(k, n)
	if km.Sign() == 0 {
		return ECDSASignature{}, errors.New("ecc: ECDSASign requires k in [1, n)")
	}
	// r = (x-coordinate of k*G) mod n
	rp := c.ScalarMul(km, g)
	if rp.Infinity {
		return ECDSASignature{}, errors.New("ecc: ECDSASign produced point at infinity, retry with new k")
	}
	r := eccMod(rp.X, n)
	if r.Sign() == 0 {
		return ECDSASignature{}, errors.New("ecc: ECDSASign produced r = 0, retry with new k")
	}
	// s = k^-1 * (e + r*d) mod n
	kInv, ok := ModInverse(km, n)
	if !ok {
		return ECDSASignature{}, errors.New("ecc: ECDSASign could not invert k modulo n")
	}
	rd := ModMul(r, eccMod(d, n), n)
	s := ModMul(kInv, ModAdd(eccMod(e, n), rd, n), n)
	if s.Sign() == 0 {
		return ECDSASignature{}, errors.New("ecc: ECDSASign produced s = 0, retry with new k")
	}
	return ECDSASignature{R: r, S: s}, nil
}

// ECDSAVerify checks an ECDSA signature (r, s) against the message digest e and
// the public key Q, using the base point G of prime subgroup order n. It returns
// true exactly when the signature is valid. The signature components must lie in
// [1, n); values outside that range are rejected.
func (c *CurveFp) ECDSAVerify(g, q PointFp, n, e *big.Int, sig ECDSASignature) bool {
	if n.Sign() <= 0 || sig.R == nil || sig.S == nil {
		return false
	}
	if sig.R.Sign() <= 0 || sig.R.Cmp(n) >= 0 || sig.S.Sign() <= 0 || sig.S.Cmp(n) >= 0 {
		return false
	}
	if q.Infinity || !c.IsOnCurve(q) {
		return false
	}
	sInv, ok := ModInverse(sig.S, n)
	if !ok {
		return false
	}
	u1 := ModMul(eccMod(e, n), sInv, n)
	u2 := ModMul(sig.R, sInv, n)
	point := c.Add(c.ScalarMul(u1, g), c.ScalarMul(u2, q))
	if point.Infinity {
		return false
	}
	v := eccMod(point.X, n)
	return v.Cmp(eccMod(sig.R, n)) == 0
}
