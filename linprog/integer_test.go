package linprog

import "testing"

func TestSolveIntegerMaximize(t *testing.T) {
	// maximize 5x1+4x2 s.t. 6x1+4x2<=24, x1+2x2<=6, integer.
	// LP relaxation optimum 21 at (3,1.5); integer optimum 20 at (4,0).
	lp := NewLP(Maximize,
		[]float64{5, 4},
		[][]float64{{6, 4}, {1, 2}},
		[]Relation{LessEqual, LessEqual},
		[]float64{24, 6},
	)
	relax := lp.Solve()
	if !linprogClose(relax.Objective, 21, tol) {
		t.Errorf("relaxation obj = %v, want 21", relax.Objective)
	}
	sol := lp.SolveInteger([]bool{true, true})
	if sol.Status != StatusOptimal {
		t.Fatalf("integer status = %v", sol.Status)
	}
	if !linprogClose(sol.Objective, 20, tol) {
		t.Errorf("integer obj = %v, want 20", sol.Objective)
	}
	if !linprogVecClose(sol.X, []float64{4, 0}, tol) {
		t.Errorf("integer x = %v, want [4 0]", sol.X)
	}
	if !IsIntegral(sol.X, tol) {
		t.Errorf("solution not integral: %v", sol.X)
	}
}

func TestSolveBinary(t *testing.T) {
	// maximize 3x1+2x2+4x3 s.t. x1+x2+x3<=2, x binary.
	// Pick the two most valuable items: x1 and x3, obj 7.
	lp := NewLP(Maximize,
		[]float64{3, 2, 4},
		[][]float64{{1, 1, 1}},
		[]Relation{LessEqual},
		[]float64{2},
	)
	sol := lp.SolveBinary([]bool{true, true, true})
	if sol.Status != StatusOptimal {
		t.Fatalf("status = %v", sol.Status)
	}
	if !linprogClose(sol.Objective, 7, tol) {
		t.Errorf("obj = %v, want 7", sol.Objective)
	}
	if !linprogVecClose(sol.X, []float64{1, 0, 1}, tol) {
		t.Errorf("x = %v, want [1 0 1]", sol.X)
	}
}

func TestSolveIntegerMinimize(t *testing.T) {
	// minimize x1+x2 s.t. 2x1+2x2>=7 (i.e. x1+x2>=3.5), integer -> x1+x2>=4.
	lp := NewLP(Minimize,
		[]float64{1, 1},
		[][]float64{{2, 2}},
		[]Relation{GreaterEqual},
		[]float64{7},
	)
	sol := lp.SolveInteger([]bool{true, true})
	if sol.Status != StatusOptimal {
		t.Fatalf("status = %v", sol.Status)
	}
	if !linprogClose(sol.Objective, 4, tol) {
		t.Errorf("obj = %v, want 4", sol.Objective)
	}
}
