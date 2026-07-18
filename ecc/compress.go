package ecc

import (
	"fmt"
	"math/big"
)

// CompressedPoint holds the SEC1-style compressed representation of an affine
// point: its x-coordinate together with a single parity bit selecting one of
// the two possible y values.
type CompressedPoint struct {
	// X is the affine x-coordinate reduced into [0, P).
	X *big.Int
	// YOdd is the least-significant bit of the affine y-coordinate, used to
	// recover the correct square root during decompression.
	YOdd bool
	// Infinity marks the compressed encoding of the point at infinity.
	Infinity bool
}

// Compress returns the compressed representation of pt, storing its
// x-coordinate and the parity of its y-coordinate. The point at infinity is
// encoded with the Infinity flag set.
func (c *CurveFp) Compress(pt PointFp) CompressedPoint {
	if pt.Infinity {
		return CompressedPoint{Infinity: true}
	}
	y := eccMod(pt.Y, c.P)
	return CompressedPoint{
		X:    new(big.Int).Set(eccMod(pt.X, c.P)),
		YOdd: y.Bit(0) == 1,
	}
}

// Decompress reconstructs the affine point from a compressed representation by
// solving y^2 = x^3 + A*x + B for y and selecting the root matching the stored
// parity bit. It returns an error if the x-coordinate is not the abscissa of
// any curve point (the right-hand side is a non-residue).
func (c *CurveFp) Decompress(cp CompressedPoint) (PointFp, error) {
	if cp.Infinity {
		return PointFp{Infinity: true}, nil
	}
	rhs := c.RHS(cp.X)
	y, ok := ModSqrt(rhs, c.P)
	if !ok {
		return PointFp{}, fmt.Errorf("ecc: x = %s is not on the curve", cp.X)
	}
	if (y.Bit(0) == 1) != cp.YOdd {
		y = ModNeg(y, c.P)
	}
	return PointFp{X: new(big.Int).Set(eccMod(cp.X, c.P)), Y: y}, nil
}

// CompressBytes returns the SEC1 compressed octet string of pt: a single byte
// 0x00 for the point at infinity, or the prefix 0x02 (even y) or 0x03 (odd y)
// followed by the x-coordinate encoded as a big-endian fixed-length field
// element. The field element length is ceil(bitlen(P)/8).
func (c *CurveFp) CompressBytes(pt PointFp) []byte {
	if pt.Infinity {
		return []byte{0x00}
	}
	size := (c.P.BitLen() + 7) / 8
	cp := c.Compress(pt)
	prefix := byte(0x02)
	if cp.YOdd {
		prefix = 0x03
	}
	out := make([]byte, 1+size)
	out[0] = prefix
	cp.X.FillBytes(out[1:])
	return out
}

// DecompressBytes parses a SEC1 compressed octet string produced by
// CompressBytes and returns the corresponding affine point. It returns an error
// for a malformed prefix, an incorrect length, or an x-coordinate that is not
// on the curve.
func (c *CurveFp) DecompressBytes(b []byte) (PointFp, error) {
	if len(b) == 1 && b[0] == 0x00 {
		return PointFp{Infinity: true}, nil
	}
	size := (c.P.BitLen() + 7) / 8
	if len(b) != 1+size {
		return PointFp{}, fmt.Errorf("ecc: compressed point has length %d, want %d", len(b), 1+size)
	}
	if b[0] != 0x02 && b[0] != 0x03 {
		return PointFp{}, fmt.Errorf("ecc: invalid compressed point prefix 0x%02x", b[0])
	}
	cp := CompressedPoint{
		X:    new(big.Int).SetBytes(b[1:]),
		YOdd: b[0] == 0x03,
	}
	return c.Decompress(cp)
}
