package linprog

import "testing"

func TestHungarianMin(t *testing.T) {
	cost := [][]float64{
		{4, 1, 3},
		{2, 0, 5},
		{3, 2, 2},
	}
	assign, total, err := Hungarian(cost)
	if err != nil {
		t.Fatalf("err = %v", err)
	}
	// Optimal assignment row0->1, row1->0, row2->2, total 1+2+2 = 5.
	if !linprogClose(total, 5, tol) {
		t.Errorf("total = %v, want 5", total)
	}
	if AssignmentCost(cost, assign) != total {
		t.Errorf("AssignmentCost mismatch: %v vs %v", AssignmentCost(cost, assign), total)
	}
	// Verify it is a valid permutation.
	seen := make([]bool, 3)
	for _, j := range assign {
		if seen[j] {
			t.Errorf("assignment not a permutation: %v", assign)
		}
		seen[j] = true
	}
}

func TestHungarianMax(t *testing.T) {
	value := [][]float64{
		{4, 1, 3},
		{2, 0, 5},
		{3, 2, 2},
	}
	_, total, err := HungarianMax(value)
	if err != nil {
		t.Fatalf("err = %v", err)
	}
	// Best assignment row0->0, row1->2, row2->1, total 4+5+2 = 11.
	if !linprogClose(total, 11, tol) {
		t.Errorf("total = %v, want 11", total)
	}
}

func TestHungarianNotSquare(t *testing.T) {
	_, _, err := Hungarian([][]float64{{1, 2, 3}, {4, 5, 6}})
	if err != ErrNotSquare {
		t.Errorf("err = %v, want ErrNotSquare", err)
	}
}

func TestTransportation(t *testing.T) {
	supply := []float64{30, 40}
	demand := []float64{20, 50}
	cost := [][]float64{
		{2, 3},
		{4, 1},
	}
	flow, total, err := Transportation(supply, demand, cost)
	if err != nil {
		t.Fatalf("err = %v", err)
	}
	// Optimal cost 110: x11=20, x12=10, x22=40.
	if !linprogClose(total, 110, tol) {
		t.Errorf("total = %v, want 110", total)
	}
	// Verify supply and demand are met exactly.
	for i := range supply {
		var s float64
		for j := range demand {
			s += flow[i][j]
		}
		if !linprogClose(s, supply[i], tol) {
			t.Errorf("row %d sum = %v, want %v", i, s, supply[i])
		}
	}
	for j := range demand {
		var d float64
		for i := range supply {
			d += flow[i][j]
		}
		if !linprogClose(d, demand[j], tol) {
			t.Errorf("col %d sum = %v, want %v", j, d, demand[j])
		}
	}
}

func TestTransportationUnbalanced(t *testing.T) {
	_, _, err := Transportation([]float64{10}, []float64{5}, [][]float64{{1}})
	if err != ErrUnbalanced {
		t.Errorf("err = %v, want ErrUnbalanced", err)
	}
}

func TestNorthWestCorner(t *testing.T) {
	flow := NorthWestCorner([]float64{30, 40}, []float64{20, 50})
	// Row sums must equal supply, column sums equal demand.
	if !linprogClose(flow[0][0]+flow[0][1], 30, tol) {
		t.Errorf("row 0 sum wrong: %v", flow[0])
	}
	if !linprogClose(flow[1][0]+flow[1][1], 40, tol) {
		t.Errorf("row 1 sum wrong: %v", flow[1])
	}
	if !IsBalanced([]float64{30, 40}, []float64{20, 50}, tol) {
		t.Errorf("problem should be balanced")
	}
}
