package linprog

import (
	"errors"
	"math"
)

// ErrUnbalanced is returned by [Transportation] when total supply does not
// equal total demand.
var ErrUnbalanced = errors.New("linprog: unbalanced transportation problem")

// ErrNotSquare is returned by [Hungarian] when the cost matrix is not square.
var ErrNotSquare = errors.New("linprog: cost matrix must be square")

// Hungarian solves the assignment problem: given an n-by-n cost matrix, it
// finds a permutation assignment of rows to columns that minimizes the total
// cost. It runs the O(n^3) Kuhn-Munkres algorithm with dual potentials.
//
// It returns assign, where assign[i] is the column assigned to row i, and the
// minimized total cost. A non-square matrix yields [ErrNotSquare].
func Hungarian(cost [][]float64) (assign []int, total float64, err error) {
	n := len(cost)
	for _, row := range cost {
		if len(row) != n {
			return nil, 0, ErrNotSquare
		}
	}
	if n == 0 {
		return []int{}, 0, nil
	}
	const inf = math.MaxFloat64
	u := make([]float64, n+1)
	v := make([]float64, n+1)
	p := make([]int, n+1) // p[j] = row matched to column j (1-indexed)
	way := make([]int, n+1)
	for i := 1; i <= n; i++ {
		p[0] = i
		j0 := 0
		minv := make([]float64, n+1)
		used := make([]bool, n+1)
		for j := range minv {
			minv[j] = inf
		}
		for {
			used[j0] = true
			i0 := p[j0]
			delta := inf
			j1 := -1
			for j := 1; j <= n; j++ {
				if used[j] {
					continue
				}
				cur := cost[i0-1][j-1] - u[i0] - v[j]
				if cur < minv[j] {
					minv[j] = cur
					way[j] = j0
				}
				if minv[j] < delta {
					delta = minv[j]
					j1 = j
				}
			}
			for j := 0; j <= n; j++ {
				if used[j] {
					u[p[j]] += delta
					v[j] -= delta
				} else {
					minv[j] -= delta
				}
			}
			j0 = j1
			if p[j0] == 0 {
				break
			}
		}
		for j0 != 0 {
			j1 := way[j0]
			p[j0] = p[j1]
			j0 = j1
		}
	}
	assign = make([]int, n)
	for j := 1; j <= n; j++ {
		if p[j] != 0 {
			assign[p[j]-1] = j - 1
		}
	}
	total = 0
	for i := 0; i < n; i++ {
		total += cost[i][assign[i]]
	}
	return assign, total, nil
}

// HungarianMax solves the maximization assignment problem by negating the cost
// matrix, running [Hungarian], and returning the assignment together with the
// maximized total value.
func HungarianMax(value [][]float64) (assign []int, total float64, err error) {
	n := len(value)
	neg := make([][]float64, n)
	for i, row := range value {
		if len(row) != n {
			return nil, 0, ErrNotSquare
		}
		neg[i] = make([]float64, n)
		for j, v := range row {
			neg[i][j] = -v
		}
	}
	assign, negTotal, err := Hungarian(neg)
	if err != nil {
		return nil, 0, err
	}
	return assign, -negTotal, nil
}

// AssignmentCost returns the total cost of the assignment, summing
// cost[i][assign[i]] over all rows.
func AssignmentCost(cost [][]float64, assign []int) float64 {
	var t float64
	for i, j := range assign {
		t += cost[i][j]
	}
	return t
}

// IsBalanced reports whether total supply equals total demand within tol.
func IsBalanced(supply, demand []float64, tol float64) bool {
	var s, d float64
	for _, v := range supply {
		s += v
	}
	for _, v := range demand {
		d += v
	}
	return math.Abs(s-d) <= tol
}

// NorthWestCorner returns an initial basic feasible flow for a balanced
// transportation problem using the north-west-corner rule. The result is an
// m-by-n matrix whose row sums equal supply and whose column sums equal
// demand.
func NorthWestCorner(supply, demand []float64) [][]float64 {
	m := len(supply)
	n := len(demand)
	s := append([]float64(nil), supply...)
	d := append([]float64(nil), demand...)
	flow := make([][]float64, m)
	for i := range flow {
		flow[i] = make([]float64, n)
	}
	i, j := 0, 0
	for i < m && j < n {
		q := math.Min(s[i], d[j])
		flow[i][j] = q
		s[i] -= q
		d[j] -= q
		if s[i] <= d[j] {
			i++
		} else {
			j++
		}
	}
	return flow
}

// Transportation solves the balanced transportation problem: ship from m
// sources with the given supplies to n sinks with the given demands at the
// given per-unit costs, minimizing total cost while meeting every supply and
// demand exactly. It reduces the problem to a [StandardLP] and solves it with
// [Simplex].
//
// It returns the optimal flow matrix (m-by-n) and total cost. Total supply
// must equal total demand or [ErrUnbalanced] is returned.
func Transportation(supply, demand []float64, cost [][]float64) (flow [][]float64, total float64, err error) {
	m := len(supply)
	n := len(demand)
	if len(cost) != m {
		return nil, 0, ErrDimension
	}
	for _, row := range cost {
		if len(row) != n {
			return nil, 0, ErrDimension
		}
	}
	if !IsBalanced(supply, demand, 1e-9) {
		return nil, 0, ErrUnbalanced
	}
	// Variable index for (i,j) is i*n + j.
	nv := m * n
	c := make([]float64, nv)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			c[i*n+j] = cost[i][j]
		}
	}
	rows := make([][]float64, 0, m+n)
	b := make([]float64, 0, m+n)
	// Supply equalities: sum_j x_ij = supply_i.
	for i := 0; i < m; i++ {
		row := make([]float64, nv)
		for j := 0; j < n; j++ {
			row[i*n+j] = 1
		}
		rows = append(rows, row)
		b = append(b, supply[i])
	}
	// Demand equalities: sum_i x_ij = demand_j.
	for j := 0; j < n; j++ {
		row := make([]float64, nv)
		for i := 0; i < m; i++ {
			row[i*n+j] = 1
		}
		rows = append(rows, row)
		b = append(b, demand[j])
	}
	res := Simplex(StandardLP{C: c, A: rows, B: b})
	if res.Status != StatusOptimal {
		return nil, 0, errors.New("linprog: transportation solve failed: " + res.Status.String())
	}
	flow = make([][]float64, m)
	for i := 0; i < m; i++ {
		flow[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			flow[i][j] = res.X[i*n+j]
		}
	}
	return flow, res.Objective, nil
}
