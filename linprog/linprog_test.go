package linprog

import (
	"math"
	"testing"
)

// linprogClose reports whether a and b agree to within tol.
func linprogClose(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

// linprogVecClose reports whether every component of a and b agrees within tol.
func linprogVecClose(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if math.Abs(a[i]-b[i]) > tol {
			return false
		}
	}
	return true
}

const tol = 1e-6

func TestSimplexMaximize(t *testing.T) {
	// maximize 3x + 5y s.t. x<=4, 2y<=12, 3x+2y<=18, x,y>=0.
	// Classic result: x=2, y=6, objective 36.
	lp := NewLP(Maximize,
		[]float64{3, 5},
		[][]float64{{1, 0}, {0, 2}, {3, 2}},
		[]Relation{LessEqual, LessEqual, LessEqual},
		[]float64{4, 12, 18},
	)
	sol := lp.Solve()
	if sol.Status != StatusOptimal {
		t.Fatalf("status = %v, want Optimal", sol.Status)
	}
	if !linprogClose(sol.Objective, 36, tol) {
		t.Errorf("objective = %v, want 36", sol.Objective)
	}
	if !linprogVecClose(sol.X, []float64{2, 6}, tol) {
		t.Errorf("x = %v, want [2 6]", sol.X)
	}
	if !lp.Feasible(sol.X, tol) {
		t.Errorf("solution not feasible: %v", sol.X)
	}
}

func TestSimplexMinimize(t *testing.T) {
	// minimize x+y s.t. x+2y>=3, 2x+y>=3, x,y>=0. Optimum (1,1), obj 2.
	lp := NewLP(Minimize,
		[]float64{1, 1},
		[][]float64{{1, 2}, {2, 1}},
		[]Relation{GreaterEqual, GreaterEqual},
		[]float64{3, 3},
	)
	sol := lp.Solve()
	if sol.Status != StatusOptimal {
		t.Fatalf("status = %v", sol.Status)
	}
	if !linprogClose(sol.Objective, 2, tol) {
		t.Errorf("objective = %v, want 2", sol.Objective)
	}
	if !linprogVecClose(sol.X, []float64{1, 1}, tol) {
		t.Errorf("x = %v, want [1 1]", sol.X)
	}
}

func TestSimplexUnbounded(t *testing.T) {
	// maximize x s.t. x - y <= 1, x,y>=0 is unbounded.
	lp := NewLP(Maximize,
		[]float64{1, 0},
		[][]float64{{1, -1}},
		[]Relation{LessEqual},
		[]float64{1},
	)
	sol := lp.Solve()
	if sol.Status != StatusUnbounded {
		t.Fatalf("status = %v, want Unbounded", sol.Status)
	}
}

func TestSimplexInfeasible(t *testing.T) {
	// x <= 1 and x >= 2 has empty feasible region.
	lp := NewLP(Minimize,
		[]float64{1},
		[][]float64{{1}, {1}},
		[]Relation{LessEqual, GreaterEqual},
		[]float64{1, 2},
	)
	sol := lp.Solve()
	if sol.Status != StatusInfeasible {
		t.Fatalf("status = %v, want Infeasible", sol.Status)
	}
}

func TestSimplexEqualityConstraint(t *testing.T) {
	// minimize 2x+3y s.t. x+y=10, x<=6, x,y>=0.
	// Cheaper to maximize x (coeff 2<3): x=6, y=4, obj=12+12=24.
	lp := NewLP(Minimize,
		[]float64{2, 3},
		[][]float64{{1, 1}, {1, 0}},
		[]Relation{Equal, LessEqual},
		[]float64{10, 6},
	)
	sol := lp.Solve()
	if sol.Status != StatusOptimal {
		t.Fatalf("status = %v", sol.Status)
	}
	if !linprogClose(sol.Objective, 24, tol) {
		t.Errorf("objective = %v, want 24", sol.Objective)
	}
	if !linprogVecClose(sol.X, []float64{6, 4}, tol) {
		t.Errorf("x = %v, want [6 4]", sol.X)
	}
}

func TestStandardFormRoundTrip(t *testing.T) {
	lp := NewLP(Maximize,
		[]float64{3, 5},
		[][]float64{{1, 0}, {0, 2}, {3, 2}},
		[]Relation{LessEqual, LessEqual, LessEqual},
		[]float64{4, 12, 18},
	)
	std := lp.Standard()
	// 2 structural + 3 slack columns.
	if std.NumVars() != 5 {
		t.Errorf("NumVars = %d, want 5", std.NumVars())
	}
	if std.NumConstraints() != 3 {
		t.Errorf("NumConstraints = %d, want 3", std.NumConstraints())
	}
	// Feed the optimal point plus its slacks and check equality feasibility.
	x := []float64{2, 6, 2, 0, 0} // slacks: 4-2=2, 12-12=0, 18-18=0
	if !std.Feasible(x, tol) {
		t.Errorf("standard point not feasible")
	}
}
