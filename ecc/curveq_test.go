package ecc

import (
	"math/big"
	"testing"
)

func rat(a, b int64) *big.Rat { return big.NewRat(a, b) }

func TestCurveQBasics(t *testing.T) {
	// E: y^2 = x^3 - 2 over Q, with the integral point P = (3, 5).
	c, err := NewCurveQ(rat(0, 1), rat(-2, 1))
	if err != nil {
		t.Fatalf("NewCurveQ: %v", err)
	}
	// Discriminant = -16*(27*4) = -1728, j-invariant = 0 (since A = 0).
	if c.Discriminant().Cmp(rat(-1728, 1)) != 0 {
		t.Errorf("Discriminant = %s, want -1728", c.Discriminant().RatString())
	}
	if c.JInvariant().Sign() != 0 {
		t.Errorf("JInvariant = %s, want 0", c.JInvariant().RatString())
	}
	p, err := c.NewPoint(rat(3, 1), rat(5, 1))
	if err != nil {
		t.Fatalf("NewPoint: %v", err)
	}
	if !c.IsOnCurve(p) {
		t.Errorf("(3,5) should be on the curve")
	}
}

func TestCurveQDouble(t *testing.T) {
	c, _ := NewCurveQ(rat(0, 1), rat(-2, 1))
	p, _ := c.NewPoint(rat(3, 1), rat(5, 1))
	// Hand-computed: 2P = (129/100, -383/1000).
	two := c.Double(p)
	if two.X.Cmp(rat(129, 100)) != 0 || two.Y.Cmp(rat(-383, 1000)) != 0 {
		t.Errorf("2P = (%s,%s), want (129/100,-383/1000)", two.X.RatString(), two.Y.RatString())
	}
	if !c.IsOnCurve(two) {
		t.Errorf("2P should be on the curve")
	}
	// ScalarMul(2) must agree with Double.
	if !c.Equal(c.ScalarMul(big.NewInt(2), p), two) {
		t.Errorf("ScalarMul(2,P) != Double(P)")
	}
}

func TestCurveQGroupLaw(t *testing.T) {
	c, _ := NewCurveQ(rat(0, 1), rat(-2, 1))
	p, _ := c.NewPoint(rat(3, 1), rat(5, 1))
	// P + (-P) = O.
	if !c.Add(p, c.Negate(p)).Infinity {
		t.Errorf("P + (-P) should be infinity")
	}
	// P + O = P.
	if !c.Equal(c.Add(p, c.Identity()), p) {
		t.Errorf("P + O should be P")
	}
	// (-2)P = -(2P).
	neg2 := c.ScalarMul(big.NewInt(-2), p)
	if !c.Equal(neg2, c.Negate(c.ScalarMul(big.NewInt(2), p))) {
		t.Errorf("(-2)P should equal -(2P)")
	}
	// 3P = 2P + P must stay on the curve (exact rational arithmetic).
	three := c.ScalarMul(big.NewInt(3), p)
	if three.Infinity || !c.IsOnCurve(three) {
		t.Errorf("3P should be a finite on-curve point")
	}
	if !c.Equal(three, c.Add(c.Double(p), p)) {
		t.Errorf("3P != 2P + P")
	}
}

func TestCurveQSingularRejected(t *testing.T) {
	if _, err := NewCurveQ(rat(0, 1), rat(0, 1)); err == nil {
		t.Errorf("y^2 = x^3 is singular and should be rejected")
	}
}
