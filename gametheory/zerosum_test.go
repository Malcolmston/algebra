package gametheory

import (
	"testing"
)

func TestSolveZeroSumMatchingPennies(t *testing.T) {
	sol, err := SolveZeroSum(MatchingPennies().Row)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !gametheoryApprox(sol.Value, 0, 1e-9) {
		t.Fatalf("value = %v, want 0", sol.Value)
	}
	if !gametheoryVecApprox(sol.RowStrategy, []float64{0.5, 0.5}, 1e-9) {
		t.Fatalf("row strategy = %v, want [0.5 0.5]", sol.RowStrategy)
	}
	if !gametheoryVecApprox(sol.ColStrategy, []float64{0.5, 0.5}, 1e-9) {
		t.Fatalf("col strategy = %v, want [0.5 0.5]", sol.ColStrategy)
	}
}

func TestSolveZeroSumRPS(t *testing.T) {
	// Rock-paper-scissors: value 0, uniform optimal strategies.
	rps := [][]float64{
		{0, -1, 1},
		{1, 0, -1},
		{-1, 1, 0},
	}
	sol, err := SolveZeroSum(rps)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !gametheoryApprox(sol.Value, 0, 1e-9) {
		t.Fatalf("RPS value = %v, want 0", sol.Value)
	}
	third := []float64{1.0 / 3, 1.0 / 3, 1.0 / 3}
	if !gametheoryVecApprox(sol.RowStrategy, third, 1e-9) {
		t.Fatalf("RPS row strategy = %v, want uniform", sol.RowStrategy)
	}
}

func TestSolveZeroSumSaddle(t *testing.T) {
	// Game with a pure saddle point. Row mins: {3,1} -> maximin row 0 = 3.
	// Col maxes: {4,3} -> minimax col 1 = 3. Saddle at (0,1), value 3.
	payoff := [][]float64{
		{4, 3},
		{2, 1},
	}
	i, j, ok := SaddlePoint(payoff)
	if !ok || i != 0 || j != 1 {
		t.Fatalf("saddle = (%d,%d,%v), want (0,1,true)", i, j, ok)
	}
	sol, err := SolveZeroSum(payoff)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !gametheoryApprox(sol.Value, 3, 1e-9) {
		t.Fatalf("value = %v, want 3", sol.Value)
	}
}

func TestMaximinMinimax(t *testing.T) {
	payoff := [][]float64{
		{4, 3},
		{2, 1},
	}
	if row, v := Maximin(payoff); row != 0 || !gametheoryApprox(v, 3, 1e-12) {
		t.Fatalf("maximin = (%d,%v), want (0,3)", row, v)
	}
	if col, v := Minimax(payoff); col != 1 || !gametheoryApprox(v, 3, 1e-12) {
		t.Fatalf("minimax = (%d,%v), want (1,3)", col, v)
	}
}

func TestSolveZeroSumMethodRejectsNonZeroSum(t *testing.T) {
	pd := StandardPrisonersDilemma().Game()
	if _, err := pd.SolveZeroSum(); err != ErrNotZeroSum {
		t.Fatalf("expected ErrNotZeroSum, got %v", err)
	}
}

func TestSolveZeroSumMixed(t *testing.T) {
	// A 2x2 zero-sum game with a known mixed value. Row payoff:
	//   [[ 3, -1], [-1, 1]]. Optimal value = 1/3, with row mixing (1/3, 2/3)?
	// Compute value via formula: for 2x2 with no saddle, v = (ad - bc)/(a-b-c+d)
	// where matrix = [[a,b],[c,d]] = [[3,-1],[-1,1]].
	// v = (3*1 - (-1)(-1)) / (3 - (-1) - (-1) + 1) = (3-1)/(6) = 1/3.
	payoff := [][]float64{{3, -1}, {-1, 1}}
	sol, err := SolveZeroSum(payoff)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !gametheoryApprox(sol.Value, 1.0/3, 1e-9) {
		t.Fatalf("value = %v, want 1/3", sol.Value)
	}
	// Verify optimality: row strategy guarantees at least the value against both
	// pure columns.
	for j := 0; j < 2; j++ {
		var v float64
		for i := 0; i < 2; i++ {
			v += sol.RowStrategy[i] * payoff[i][j]
		}
		if v < sol.Value-1e-9 {
			t.Fatalf("row strategy under-performs against column %d: %v < %v", j, v, sol.Value)
		}
	}
}
