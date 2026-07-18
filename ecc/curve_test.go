package ecc

import (
	"math/big"
	"testing"
)

// eccTestCurve is the standard textbook curve E: y^2 = x^3 + 2x + 2 over
// GF(17). It is cyclic of prime order 19 with generator G = (5, 1).
func eccTestCurve(t *testing.T) (*CurveFp, PointFp) {
	t.Helper()
	c, err := NewCurveFp(bi(2), bi(2), bi(17))
	if err != nil {
		t.Fatalf("NewCurveFp: %v", err)
	}
	g, err := c.NewPoint(bi(5), bi(1))
	if err != nil {
		t.Fatalf("NewPoint generator: %v", err)
	}
	return c, g
}

// eccMultiples is the known table k -> k*G for the test curve, k = 1..19.
var eccMultiples = [][2]int64{
	{5, 1}, {6, 3}, {10, 6}, {3, 1}, {9, 16}, {16, 13}, {0, 6}, {13, 7},
	{7, 6}, {7, 11}, {13, 10}, {0, 11}, {16, 4}, {9, 1}, {3, 16}, {10, 11},
	{6, 14}, {5, 16},
}

func TestInvariants(t *testing.T) {
	c, _ := eccTestCurve(t)
	if got := c.Discriminant(); got.Cmp(bi(4)) != 0 {
		t.Errorf("Discriminant = %s, want 4", got)
	}
	if got := c.JInvariant(); got.Cmp(bi(3)) != 0 {
		t.Errorf("JInvariant = %s, want 3", got)
	}
	if !c.IsSmooth() {
		t.Errorf("curve should be smooth")
	}
	// Free-function forms must agree.
	if got := DiscriminantFp(bi(2), bi(2), bi(17)); got.Cmp(bi(4)) != 0 {
		t.Errorf("DiscriminantFp = %s, want 4", got)
	}
	j, ok := JInvariantFp(bi(2), bi(2), bi(17))
	if !ok || j.Cmp(bi(3)) != 0 {
		t.Errorf("JInvariantFp = %s ok=%v, want 3 true", j, ok)
	}
}

func TestSingularRejected(t *testing.T) {
	// y^2 = x^3 (a=0,b=0) is singular.
	if _, err := NewCurveFp(bi(0), bi(0), bi(17)); err == nil {
		t.Errorf("expected singular curve to be rejected")
	}
	// p must exceed 3.
	if _, err := NewCurveFp(bi(1), bi(1), bi(3)); err == nil {
		t.Errorf("expected p<=3 to be rejected")
	}
}

func TestScalarMulTable(t *testing.T) {
	c, g := eccTestCurve(t)
	for i, want := range eccMultiples {
		k := int64(i + 1)
		got := c.ScalarMul(bi(k), g)
		if got.Infinity || got.X.Cmp(bi(want[0])) != 0 || got.Y.Cmp(bi(want[1])) != 0 {
			t.Errorf("%dG = (%v,%v,inf=%v), want (%d,%d)", k, got.X, got.Y, got.Infinity, want[0], want[1])
		}
		if !c.IsOnCurve(got) {
			t.Errorf("%dG is not on the curve", k)
		}
	}
	// 19G must be the identity, and 20G = G.
	if !c.ScalarMul(bi(19), g).Infinity {
		t.Errorf("19G should be the point at infinity")
	}
	if !c.Equal(c.ScalarMul(bi(20), g), g) {
		t.Errorf("20G should equal G")
	}
}

func TestAddDoubleNegate(t *testing.T) {
	c, g := eccTestCurve(t)
	// G + G via Add must equal Double(G) = 2G = (6,3).
	if !c.Equal(c.Add(g, g), c.Double(g)) {
		t.Errorf("Add(G,G) != Double(G)")
	}
	two := c.Double(g)
	if two.X.Cmp(bi(6)) != 0 || two.Y.Cmp(bi(3)) != 0 {
		t.Errorf("2G = (%v,%v), want (6,3)", two.X, two.Y)
	}
	// G + 2G = 3G = (10,6).
	three := c.Add(g, two)
	if three.X.Cmp(bi(10)) != 0 || three.Y.Cmp(bi(6)) != 0 {
		t.Errorf("3G = (%v,%v), want (10,6)", three.X, three.Y)
	}
	// G + (-G) = infinity, and Negate matches 18G.
	neg := c.Negate(g)
	if !c.Add(g, neg).Infinity {
		t.Errorf("G + (-G) should be infinity")
	}
	if !c.Equal(neg, c.ScalarMul(bi(18), g)) {
		t.Errorf("-G should equal 18G")
	}
	// Identity behaviour.
	if !c.Equal(c.Add(g, c.Identity()), g) {
		t.Errorf("G + O should be G")
	}
	// Negative scalar: (-3)G = -(3G) = 16G.
	if !c.Equal(c.ScalarMul(bi(-3), g), c.ScalarMul(bi(16), g)) {
		t.Errorf("(-3)G should equal 16G")
	}
}

func TestNoInputMutation(t *testing.T) {
	c, g := eccTestCurve(t)
	gx := new(big.Int).Set(g.X)
	gy := new(big.Int).Set(g.Y)
	_ = c.ScalarMul(bi(7), g)
	_ = c.Add(g, g)
	if g.X.Cmp(gx) != 0 || g.Y.Cmp(gy) != 0 {
		t.Errorf("generator coordinates mutated by operations")
	}
}
