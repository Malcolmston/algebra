package ecc

import "testing"

func TestPointCounting(t *testing.T) {
	c, _ := eccTestCurve(t)
	if got := c.CountPointsNaive(); got.Cmp(bi(19)) != 0 {
		t.Errorf("CountPointsNaive = %s, want 19", got)
	}
	// trace t = p + 1 - N = 17 + 1 - 19 = -1.
	if got := c.TraceOfFrobenius(); got.Cmp(bi(-1)) != 0 {
		t.Errorf("TraceOfFrobenius = %s, want -1", got)
	}
	lo, hi := c.HasseInterval()
	if lo.Cmp(bi(19)) > 0 || hi.Cmp(bi(19)) < 0 {
		t.Errorf("Hasse interval [%s,%s] should contain 19", lo, hi)
	}
}

func TestPointOrder(t *testing.T) {
	c, g := eccTestCurve(t)
	// The group has prime order 19, so every non-identity point has order 19.
	if got := c.PointOrderNaive(g); got == nil || got.Cmp(bi(19)) != 0 {
		t.Errorf("PointOrderNaive(G) = %s, want 19", got)
	}
	if got := c.PointOrderBSGS(g); got == nil || got.Cmp(bi(19)) != 0 {
		t.Errorf("PointOrderBSGS(G) = %s, want 19", got)
	}
	// Naive and BSGS must agree for every multiple of G.
	for i := 1; i <= 19; i++ {
		p := c.ScalarMul(bi(int64(i)), g)
		a := c.PointOrderNaive(p)
		b := c.PointOrderBSGS(p)
		if a == nil || b == nil || a.Cmp(b) != 0 {
			t.Errorf("order mismatch at %dG: naive=%s bsgs=%s", i, a, b)
		}
	}
	// Identity has order 1.
	if got := c.PointOrderBSGS(c.Identity()); got.Cmp(bi(1)) != 0 {
		t.Errorf("PointOrderBSGS(O) = %s, want 1", got)
	}
}

func TestGeneratorAndCofactor(t *testing.T) {
	c, g := eccTestCurve(t)
	if !c.IsGenerator(g, bi(19)) {
		t.Errorf("G should be a generator of order 19")
	}
	if c.IsGenerator(g, bi(7)) {
		t.Errorf("G does not have order 7")
	}
	cof, ok := c.Cofactor(bi(19))
	if !ok || cof.Cmp(bi(1)) != 0 {
		t.Errorf("Cofactor(19) = %s ok=%v, want 1 true", cof, ok)
	}
	if _, ok := c.Cofactor(bi(7)); ok {
		t.Errorf("7 does not divide the group order")
	}
}
