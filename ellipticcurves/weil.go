package ellipticcurves

import (
	"errors"
	"math/big"
	"math/rand"
)

// ErrPairingDegenerate indicates that a Miller-function evaluation hit a zero of
// a line or vertical function at the evaluation point, so a fresh auxiliary
// point is required. WeilPairing retries automatically; this error surfaces only
// from the low-level MillerFunction.
var ErrPairingDegenerate = errors.New("ellipticcurves: degenerate pairing evaluation")

// ErrPairingInput indicates that the pairing inputs are not valid n-torsion
// points, or n is not positive.
var ErrPairingInput = errors.New("ellipticcurves: invalid pairing input")

// verticalValue returns the value at Z of the vertical line through R, namely
// x_Z - x_R, or 1 when R is the point at infinity.
func (c *CurveFp) verticalValue(r, z PointFp) *big.Int {
	if r.Infinity {
		return big.NewInt(1)
	}
	return ModSub(z.X, r.X, c.P)
}

// tangentValue returns the value at Z of the tangent line to the curve at T.
// For a 2-torsion point the tangent is vertical and the value is x_Z - x_T.
func (c *CurveFp) tangentValue(t, z PointFp) (*big.Int, error) {
	if t.Y.Sign() == 0 {
		return ModSub(z.X, t.X, c.P), nil
	}
	num := ModAdd(ModMul(bigThree, ModSquare(t.X, c.P), c.P), c.A, c.P)
	den := ModDouble(t.Y, c.P)
	lambda, err := ModDiv(num, den, c.P)
	if err != nil {
		return nil, err
	}
	// (y_Z - y_T) - lambda*(x_Z - x_T)
	v := ModSub(ModSub(z.Y, t.Y, c.P), ModMul(lambda, ModSub(z.X, t.X, c.P), c.P), c.P)
	return v, nil
}

// chordValue returns the value at Z of the line through T and S. When T and S
// share an x-coordinate the line is vertical and the value is x_Z - x_T.
func (c *CurveFp) chordValue(t, s, z PointFp) (*big.Int, error) {
	if t.X.Cmp(s.X) == 0 {
		return ModSub(z.X, t.X, c.P), nil
	}
	lambda, err := ModDiv(ModSub(s.Y, t.Y, c.P), ModSub(s.X, t.X, c.P), c.P)
	if err != nil {
		return nil, err
	}
	v := ModSub(ModSub(z.Y, t.Y, c.P), ModMul(lambda, ModSub(z.X, t.X, c.P), c.P), c.P)
	return v, nil
}

// MillerFunction evaluates the Miller function f_{n,P} at the affine point Z,
// where div(f_{n,P}) = n*(P) - n*(O) when n*P = O. Z must avoid the support of
// that divisor and all intermediate lines; otherwise ErrPairingDegenerate is
// returned so the caller can retry with a different evaluation point.
func (c *CurveFp) MillerFunction(n *big.Int, p, z PointFp) (*big.Int, error) {
	if z.Infinity {
		return nil, ErrPairingDegenerate
	}
	f := big.NewInt(1)
	t := clonePointFp(p)
	for i := n.BitLen() - 2; i >= 0; i-- {
		ln, err := c.tangentValue(t, z)
		if err != nil {
			return nil, ErrPairingDegenerate
		}
		t2 := c.Double(t)
		vv := c.verticalValue(t2, z)
		f = ModSquare(f, c.P)
		f = ModMul(f, ln, c.P)
		fq, err := ModDiv(f, vv, c.P)
		if err != nil {
			return nil, ErrPairingDegenerate
		}
		f = fq
		t = t2
		if n.Bit(i) == 1 {
			ln, err := c.chordValue(t, p, z)
			if err != nil {
				return nil, ErrPairingDegenerate
			}
			tsum := c.Add(t, p)
			vv := c.verticalValue(tsum, z)
			f = ModMul(f, ln, c.P)
			fq, err := ModDiv(f, vv, c.P)
			if err != nil {
				return nil, ErrPairingDegenerate
			}
			f = fq
			t = tsum
		}
	}
	if f.Sign() == 0 {
		return nil, ErrPairingDegenerate
	}
	return f, nil
}

// millerRatio evaluates f_{n,P} on the degree-zero divisor (a) - (b), returning
// f_{n,P}(a) / f_{n,P}(b).
func (c *CurveFp) millerRatio(n *big.Int, p, a, b PointFp) (*big.Int, error) {
	fa, err := c.MillerFunction(n, p, a)
	if err != nil {
		return nil, err
	}
	fb, err := c.MillerFunction(n, p, b)
	if err != nil {
		return nil, err
	}
	return ModDiv(fa, fb, c.P)
}

// WeilPairing returns the Weil pairing e_n(P, Q) of two n-torsion points as an
// element of F_p, using Miller's algorithm with random divisor shifts drawn
// from rng to avoid degeneracies. The result is an n-th root of unity in F_p.
// It returns ErrPairingInput when P or Q is not n-torsion, and
// ErrPairingDegenerate if repeated random shifts all fail (which is
// astronomically unlikely for valid independent inputs).
func (c *CurveFp) WeilPairing(n *big.Int, pPt, qPt PointFp, rng *rand.Rand) (*big.Int, error) {
	if n.Sign() <= 0 {
		return nil, ErrPairingInput
	}
	if !c.ScalarMul(n, pPt).Infinity || !c.ScalarMul(n, qPt).Infinity {
		return nil, ErrPairingInput
	}
	if pPt.Infinity || qPt.Infinity {
		return big.NewInt(1), nil
	}
	if c.PointEqual(pPt, qPt) {
		return big.NewInt(1), nil
	}
	// Weil pairing via Miller with a single auxiliary point S:
	//   e_n(P,Q) = [ f_P(Q+S)/f_P(S) ] / [ f_Q(P-S)/f_Q(-S) ].
	// The shared shift S is what makes the value independent of the choice of S
	// (Weil reciprocity); random retries only avoid the measure-zero set of
	// shifts hitting the support of the divisors.
	for attempt := 0; attempt < 200; attempt++ {
		s := c.RandomPointFp(rng64(rng, attempt))
		negS := c.Neg(s)
		fp, err := c.millerRatio(n, pPt, c.Add(qPt, s), s)
		if err != nil {
			continue
		}
		fq, err := c.millerRatio(n, qPt, c.Add(pPt, negS), negS)
		if err != nil {
			continue
		}
		res, err := ModDiv(fp, fq, c.P)
		if err != nil {
			continue
		}
		if res.Sign() == 0 {
			continue
		}
		return res, nil
	}
	return nil, ErrPairingDegenerate
}

// rng64 derives a deterministic sub-source from rng and a counter so that
// retries use fresh but reproducible randomness.
func rng64(rng *rand.Rand, counter int) *rand.Rand {
	var golden uint64 = 0x9E3779B97F4A7C15
	seed := rng.Int63() ^ (int64(counter) * int64(golden))
	return rand.New(rand.NewSource(seed))
}

// WeilPairingRootOfUnity reports whether v is an n-th root of unity in F_p, a
// sanity check on a Weil pairing value: v^n = 1.
func (c *CurveFp) WeilPairingRootOfUnity(v, n *big.Int) bool {
	return new(big.Int).Exp(v, n, c.P).Cmp(bigOne) == 0
}
