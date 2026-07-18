package linprog

import "testing"

func TestStrongDuality(t *testing.T) {
	// Primal: min 2x1+3x2 s.t. x1+x2>=4, x1+2x2>=5, x>=0.
	// Optimum (3,1), obj 9. Dual optimum also 9.
	c := []float64{2, 3}
	a := [][]float64{{1, 1}, {1, 2}}
	b := []float64{4, 5}
	primal := NewLP(Minimize, c,
		a,
		[]Relation{GreaterEqual, GreaterEqual},
		b,
	)
	ps := primal.Solve()
	if ps.Status != StatusOptimal {
		t.Fatalf("primal status = %v", ps.Status)
	}
	if !linprogClose(ps.Objective, 9, tol) {
		t.Errorf("primal obj = %v, want 9", ps.Objective)
	}
	if !linprogVecClose(ps.X, []float64{3, 1}, tol) {
		t.Errorf("primal x = %v, want [3 1]", ps.X)
	}

	dual := DualCanonical(c, a, b)
	ds := dual.Solve()
	if ds.Status != StatusOptimal {
		t.Fatalf("dual status = %v", ds.Status)
	}
	if !linprogClose(ds.Objective, 9, tol) {
		t.Errorf("dual obj = %v, want 9", ds.Objective)
	}
	if gap := DualityGap(ps.Objective, ds.Objective); !linprogClose(gap, 0, tol) {
		t.Errorf("duality gap = %v, want 0", gap)
	}
	if !WeakDualityHolds(ps.Objective, ds.Objective, tol) {
		t.Errorf("weak duality violated")
	}
}

func TestComplementarySlacknessAndKKT(t *testing.T) {
	c := []float64{2, 3}
	a := [][]float64{{1, 1}, {1, 2}}
	b := []float64{4, 5}
	x := []float64{3, 1} // primal optimum
	y := []float64{1, 1} // dual optimum
	if !ComplementarySlackness(c, a, b, x, y, tol) {
		t.Errorf("complementary slackness should hold at optimum")
	}
	k := CheckKKTLP(c, a, b, x, y)
	if !k.Satisfied(tol) {
		t.Errorf("KKT residuals not satisfied: %+v (max=%v)", k, k.Max())
	}

	// A clearly non-optimal pair should violate complementary slackness.
	if ComplementarySlackness(c, a, b, []float64{5, 0}, []float64{1, 1}, tol) {
		t.Errorf("expected complementary slackness to fail for non-optimal pair")
	}
}
