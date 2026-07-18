package linprog

import "testing"

func TestSolveQPEqualitySimple(t *testing.T) {
	// minimize 1/2(x1^2+x2^2) s.t. x1+x2=1. Optimum (0.5,0.5), obj 0.25.
	q := [][]float64{{1, 0}, {0, 1}}
	c := []float64{0, 0}
	a := [][]float64{{1, 1}}
	b := []float64{1}
	x, lambda, err := SolveQPEquality(q, c, a, b)
	if err != nil {
		t.Fatalf("err = %v", err)
	}
	if !linprogVecClose(x, []float64{0.5, 0.5}, tol) {
		t.Errorf("x = %v, want [0.5 0.5]", x)
	}
	if !linprogClose(lambda[0], -0.5, tol) {
		t.Errorf("lambda = %v, want -0.5", lambda[0])
	}
}

func TestSolveQPEqualityShifted(t *testing.T) {
	// minimize x1^2+x2^2-2x1-5x2 s.t. x1+x2=3. Optimum (0.75,2.25).
	q := [][]float64{{2, 0}, {0, 2}}
	c := []float64{-2, -5}
	a := [][]float64{{1, 1}}
	b := []float64{3}
	x, _, err := SolveQPEquality(q, c, a, b)
	if err != nil {
		t.Fatalf("err = %v", err)
	}
	if !linprogVecClose(x, []float64{0.75, 2.25}, tol) {
		t.Errorf("x = %v, want [0.75 2.25]", x)
	}
	qp := QP{Q: q, C: c, A: a, B: b}
	if !linprogClose(qp.Objective(x), -7.125, tol) {
		t.Errorf("objective = %v, want -7.125", qp.Objective(x))
	}
}

func TestSolveQPActiveSet(t *testing.T) {
	// Nocedal & Wright example 16.4:
	// minimize (x1-1)^2 + (x2-2.5)^2 with inequality constraints.
	// In 1/2 x^T Q x + c x form (constant dropped): Q=2I, c=[-2,-5].
	// Solution x* = (1.4, 1.7).
	qp := QP{
		Q: [][]float64{{2, 0}, {0, 2}},
		C: []float64{-2, -5},
		G: [][]float64{
			{-1, 2}, // -x1+2x2 <= 2
			{1, 2},  //  x1+2x2 <= 6
			{1, -2}, //  x1-2x2 <= 2
			{-1, 0}, // -x1     <= 0
			{0, -1}, //     -x2 <= 0
		},
		H: []float64{2, 6, 2, 0, 0},
	}
	x0 := []float64{2, 0} // feasible starting vertex
	if !qp.Feasible(x0, 1e-9) {
		t.Fatalf("starting point not feasible")
	}
	sol := SolveQP(qp, x0)
	if sol.Status != StatusOptimal {
		t.Fatalf("status = %v", sol.Status)
	}
	if !linprogVecClose(sol.X, []float64{1.4, 1.7}, 1e-5) {
		t.Errorf("x = %v, want [1.4 1.7]", sol.X)
	}
	// Objective in the reduced form equals 0.8 - 7.25 = -6.45.
	if !linprogClose(sol.Objective, -6.45, 1e-5) {
		t.Errorf("objective = %v, want -6.45", sol.Objective)
	}
	// KKT residuals must vanish at the reported solution.
	k := ComputeKKTResidual(qp, sol.X, sol.LambdaEq, sol.MuIneq)
	if !k.Satisfied(1e-5) {
		t.Errorf("KKT residuals not satisfied: %+v (max=%v)", k, k.Max())
	}
}

func TestKKTResidualNonzero(t *testing.T) {
	qp := QP{
		Q: [][]float64{{2, 0}, {0, 2}},
		C: []float64{-2, -5},
		G: [][]float64{{-1, 2}},
		H: []float64{2},
	}
	// A clearly non-stationary point should show a large residual.
	k := ComputeKKTResidual(qp, []float64{0, 0}, nil, []float64{0})
	if k.Satisfied(1e-3) {
		t.Errorf("expected non-trivial residual at non-optimal point, got %+v", k)
	}
}

func TestLinearSolver(t *testing.T) {
	// Solve [[2,1],[1,3]] x = [3,5] -> x = [0.8, 1.4].
	x, err := SolveLinear([][]float64{{2, 1}, {1, 3}}, []float64{3, 5})
	if err != nil {
		t.Fatalf("err = %v", err)
	}
	if !linprogVecClose(x, []float64{0.8, 1.4}, tol) {
		t.Errorf("x = %v, want [0.8 1.4]", x)
	}
	// Singular matrix must be reported.
	if _, err := SolveLinear([][]float64{{1, 1}, {1, 1}}, []float64{1, 2}); err != ErrSingular {
		t.Errorf("err = %v, want ErrSingular", err)
	}
}

func TestLinearAlgebraHelpers(t *testing.T) {
	if got := Dot([]float64{1, 2, 3}, []float64{4, 5, 6}); !linprogClose(got, 32, tol) {
		t.Errorf("Dot = %v, want 32", got)
	}
	a := [][]float64{{1, 2}, {3, 4}}
	if got := MatVec(a, []float64{1, 1}); !linprogVecClose(got, []float64{3, 7}, tol) {
		t.Errorf("MatVec = %v, want [3 7]", got)
	}
	if got := MatTVec(a, []float64{1, 1}); !linprogVecClose(got, []float64{4, 6}, tol) {
		t.Errorf("MatTVec = %v, want [4 6]", got)
	}
	tr := Transpose(a)
	if tr[0][1] != 3 || tr[1][0] != 2 {
		t.Errorf("Transpose wrong: %v", tr)
	}
}

// BenchmarkTransportation solves a deterministic 8x8 balanced transportation
// problem, which is the heaviest routine (it builds and runs the two-phase
// simplex on a 16-constraint, 64-variable standard program).
func BenchmarkTransportation(b *testing.B) {
	const n = 8
	supply := make([]float64, n)
	demand := make([]float64, n)
	cost := make([][]float64, n)
	for i := 0; i < n; i++ {
		supply[i] = 10
		demand[i] = 10
		cost[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			// Deterministic cost pattern.
			cost[i][j] = float64((i*7+j*3)%11 + 1)
		}
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, _, err := Transportation(supply, demand, cost); err != nil {
			b.Fatalf("err = %v", err)
		}
	}
}
