package ecc

import "testing"

func TestCompressRoundTrip(t *testing.T) {
	c, g := eccTestCurve(t)
	for i := 1; i <= 19; i++ {
		p := c.ScalarMul(bi(int64(i)), g)
		cp := c.Compress(p)
		got, err := c.Decompress(cp)
		if err != nil {
			t.Fatalf("Decompress %dG: %v", i, err)
		}
		if !c.Equal(got, p) {
			t.Errorf("compress round-trip %dG: got (%v,%v,inf=%v)", i, got.X, got.Y, got.Infinity)
		}
	}
}

func TestCompressParity(t *testing.T) {
	c, g := eccTestCurve(t)
	// G = (5,1): y = 1 is odd.
	if cp := c.Compress(g); !cp.YOdd {
		t.Errorf("G should compress with odd parity")
	}
	// 7G = (0,6): y = 6 is even.
	p7 := c.ScalarMul(bi(7), g)
	if cp := c.Compress(p7); cp.YOdd {
		t.Errorf("7G should compress with even parity")
	}
}

func TestCompressBytes(t *testing.T) {
	c, g := eccTestCurve(t)
	b := c.CompressBytes(g)
	// prefix 0x03 (odd y) plus one x byte for p = 17.
	if len(b) != 2 || b[0] != 0x03 || b[1] != 5 {
		t.Errorf("CompressBytes(G) = %x, want 0305", b)
	}
	got, err := c.DecompressBytes(b)
	if err != nil {
		t.Fatalf("DecompressBytes: %v", err)
	}
	if !c.Equal(got, g) {
		t.Errorf("DecompressBytes round-trip mismatch")
	}
	// Point at infinity encodes as a single zero byte.
	inf := c.CompressBytes(c.Identity())
	if len(inf) != 1 || inf[0] != 0x00 {
		t.Errorf("infinity encoding = %x, want 00", inf)
	}
	infPt, err := c.DecompressBytes(inf)
	if err != nil || !infPt.Infinity {
		t.Errorf("infinity round-trip failed: %v", err)
	}
	// Malformed inputs are rejected.
	if _, err := c.DecompressBytes([]byte{0x05, 0x05}); err == nil {
		t.Errorf("bad prefix should error")
	}
	if _, err := c.DecompressBytes([]byte{0x02, 0x05, 0x05}); err == nil {
		t.Errorf("bad length should error")
	}
}
