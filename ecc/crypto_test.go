package ecc

import (
	"math/big"
	"testing"
)

// eccByteStream is a deterministic io.Reader yielding a fixed repeating pattern,
// used to make key generation reproducible in tests.
type eccByteStream struct {
	pat []byte
	i   int
}

func (s *eccByteStream) Read(b []byte) (int, error) {
	for i := range b {
		b[i] = s.pat[s.i%len(s.pat)]
		s.i++
	}
	return len(b), nil
}

func TestECDHSharedSecret(t *testing.T) {
	c, g := eccTestCurve(t)
	// Alice d=11, Bob d=7 on the order-19 subgroup.
	aPriv := bi(11)
	bPriv := bi(7)
	aPub := c.PublicKeyFp(g, aPriv) // 11G = (13,10)
	bPub := c.PublicKeyFp(g, bPriv) // 7G  = (0,6)
	if aPub.X.Cmp(bi(13)) != 0 || aPub.Y.Cmp(bi(10)) != 0 {
		t.Fatalf("Alice public = (%v,%v), want (13,10)", aPub.X, aPub.Y)
	}
	sA, ptA, okA := c.ECDHSharedSecret(aPriv, bPub)
	sB, _, okB := c.ECDHSharedSecret(bPriv, aPub)
	if !okA || !okB {
		t.Fatalf("ECDH failed okA=%v okB=%v", okA, okB)
	}
	if sA.Cmp(sB) != 0 {
		t.Errorf("shared secrets differ: %s vs %s", sA, sB)
	}
	// 11*7 = 77 ≡ 1 (mod 19), so the shared point is 1G = (5,1), secret x = 5.
	if sA.Cmp(bi(5)) != 0 {
		t.Errorf("shared secret = %s, want 5", sA)
	}
	if !c.IsOnCurve(ptA) {
		t.Errorf("shared point not on curve")
	}
}

func TestECDSAKnownAnswer(t *testing.T) {
	c, g := eccTestCurve(t)
	n := bi(19)
	d := bi(11)
	e := bi(5)
	k := bi(3)
	sig, err := c.ECDSASign(g, n, d, e, k)
	if err != nil {
		t.Fatalf("ECDSASign: %v", err)
	}
	// Hand-computed: r = x(3G) mod 19 = 10, s = k^-1 (e + r d) mod 19 = 13.
	if sig.R.Cmp(bi(10)) != 0 || sig.S.Cmp(bi(13)) != 0 {
		t.Errorf("signature = (%s,%s), want (10,13)", sig.R, sig.S)
	}
	q := c.PublicKeyFp(g, d)
	if !c.ECDSAVerify(g, q, n, e, sig) {
		t.Errorf("valid signature failed to verify")
	}
	// Wrong message must fail.
	if c.ECDSAVerify(g, q, n, bi(6), sig) {
		t.Errorf("signature verified against wrong message")
	}
	// Tampered s must fail.
	bad := ECDSASignature{R: sig.R, S: new(big.Int).Add(sig.S, bi(1))}
	if c.ECDSAVerify(g, q, n, e, bad) {
		t.Errorf("tampered signature verified")
	}
	// Out-of-range r rejected.
	if c.ECDSAVerify(g, q, n, e, ECDSASignature{R: bi(0), S: sig.S}) {
		t.Errorf("r=0 should be rejected")
	}
}

func TestECDSARoundTripAllNonces(t *testing.T) {
	c, g := eccTestCurve(t)
	n := bi(19)
	d := bi(9)
	q := c.PublicKeyFp(g, d)
	for _, e := range []int64{1, 2, 8, 15, 18} {
		for k := int64(1); k < 19; k++ {
			sig, err := c.ECDSASign(g, n, d, bi(e), bi(k))
			if err != nil {
				continue // r=0 or s=0 cases legitimately skipped
			}
			if !c.ECDSAVerify(g, q, n, bi(e), sig) {
				t.Errorf("round-trip verify failed e=%d k=%d", e, k)
			}
		}
	}
}

func TestGenerateKeyPair(t *testing.T) {
	c, g := eccTestCurve(t)
	n := bi(19)
	reader := &eccByteStream{pat: []byte{0x01, 0x04, 0x07, 0x0a, 0x0d}}
	for i := 0; i < 20; i++ {
		priv, pub, err := c.GenerateKeyPair(g, n, reader)
		if err != nil {
			t.Fatalf("GenerateKeyPair: %v", err)
		}
		if priv.Sign() <= 0 || priv.Cmp(n) >= 0 {
			t.Errorf("private key %s out of range [1,19)", priv)
		}
		if !c.Equal(pub, c.ScalarMul(priv, g)) {
			t.Errorf("public key does not match priv*G")
		}
		if !pub.Infinity && !c.IsOnCurve(pub) {
			t.Errorf("public key not on curve")
		}
	}
}

func TestNamedCurvesConsistent(t *testing.T) {
	for _, nc := range []NamedCurve{Secp256k1(), P256()} {
		if !nc.Curve.IsOnCurve(nc.G) {
			t.Errorf("named curve generator not on curve")
		}
		// N*G must be the identity for a subgroup of order N.
		if !nc.Curve.ScalarMul(nc.N, nc.G).Infinity {
			t.Errorf("N*G is not the identity")
		}
		// (N+1)*G = G.
		np1 := new(big.Int).Add(nc.N, big.NewInt(1))
		if !nc.Curve.Equal(nc.Curve.ScalarMul(np1, nc.G), nc.G) {
			t.Errorf("(N+1)*G != G")
		}
	}
}

func TestNamedCurveECDSA(t *testing.T) {
	nc := Secp256k1()
	d := eccHex("1122334455667788990011223344556677889900112233445566778899001122")
	q := nc.Curve.PublicKeyFp(nc.G, d)
	e := eccHex("abcdef0123456789abcdef0123456789abcdef0123456789abcdef0123456789")
	k := eccHex("00fedcba9876543210fedcba9876543210fedcba9876543210fedcba98765432")
	sig, err := nc.Curve.ECDSASign(nc.G, nc.N, d, e, k)
	if err != nil {
		t.Fatalf("ECDSASign secp256k1: %v", err)
	}
	if !nc.Curve.ECDSAVerify(nc.G, q, nc.N, e, sig) {
		t.Errorf("secp256k1 signature failed to verify")
	}
	if nc.Curve.ECDSAVerify(nc.G, q, nc.N, new(big.Int).Add(e, big.NewInt(1)), sig) {
		t.Errorf("secp256k1 signature verified wrong digest")
	}
}

func BenchmarkScalarMulSecp256k1(b *testing.B) {
	nc := Secp256k1()
	k := eccHex("7fedcba9876543210123456789abcdef7fedcba9876543210123456789abcdef")
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = nc.Curve.ScalarMul(k, nc.G)
	}
}
