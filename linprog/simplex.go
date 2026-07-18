package linprog

import "math"

// SimplexResult is the outcome of running the [Simplex] method on a
// [StandardLP].
type SimplexResult struct {
	// Status classifies the outcome.
	Status Status
	// X is the optimal primal vector (length NumVars) when Status is
	// StatusOptimal; otherwise it holds the last computed point.
	X []float64
	// Objective is C · X.
	Objective float64
	// Basis holds the column indices of the basic variables, one per row.
	Basis []int
	// Iterations counts the total pivots performed across both phases.
	Iterations int
}

// LPSolution is the outcome of solving a general [LP].
type LPSolution struct {
	// Status classifies the outcome.
	Status Status
	// X is the optimal setting of the original variables (length NumVars).
	X []float64
	// Objective is the value of the original objective C · X (respecting the
	// program's Sense).
	Objective float64
	// Iterations counts the simplex pivots performed.
	Iterations int
}

// linprogTableau is the working state of the primal simplex: the current
// B^{-1}A tableau, right-hand side and basis.
type linprogTableau struct {
	t     [][]float64 // m x n
	rhs   []float64   // m
	basis []int       // m
	m, n  int
}

// linprogPivot performs a Gauss-Jordan pivot on entry (row, col), updating the
// tableau, right-hand side and basis in place.
func (tab *linprogTableau) linprogPivot(row, col int) {
	p := tab.t[row][col]
	for j := 0; j < tab.n; j++ {
		tab.t[row][j] /= p
	}
	tab.rhs[row] /= p
	for i := 0; i < tab.m; i++ {
		if i == row {
			continue
		}
		f := tab.t[i][col]
		if f == 0 {
			continue
		}
		for j := 0; j < tab.n; j++ {
			tab.t[i][j] -= f * tab.t[row][j]
		}
		tab.rhs[i] -= f * tab.rhs[row]
	}
	tab.basis[row] = col
}

// linprogInBasis reports whether column j is currently basic.
func (tab *linprogTableau) linprogInBasis(j int) bool {
	for _, b := range tab.basis {
		if b == j {
			return true
		}
	}
	return false
}

// linprogOptimize runs the primal simplex on the tableau minimizing cost·x,
// using Bland's rule for both entering and leaving choices so it cannot cycle.
// Columns flagged in forbidden are never allowed to enter the basis. It
// returns the resulting Status and the number of pivots performed.
func (tab *linprogTableau) linprogOptimize(cost []float64, forbidden []bool, maxIter int) (Status, int) {
	iters := 0
	for iters < maxIter {
		// Reduced costs r[j] = cost[j] - cost_B · column_j.
		entering := -1
		for j := 0; j < tab.n; j++ {
			if forbidden != nil && forbidden[j] {
				continue
			}
			if tab.linprogInBasis(j) {
				continue
			}
			r := cost[j]
			for i := 0; i < tab.m; i++ {
				r -= cost[tab.basis[i]] * tab.t[i][j]
			}
			if r < -Eps {
				entering = j // Bland: smallest eligible index.
				break
			}
		}
		if entering == -1 {
			return StatusOptimal, iters
		}
		// Ratio test with Bland tie-break on smallest basic index.
		leave := -1
		best := math.Inf(1)
		bestBasis := math.MaxInt
		for i := 0; i < tab.m; i++ {
			if tab.t[i][entering] > Eps {
				ratio := tab.rhs[i] / tab.t[i][entering]
				if ratio < best-Eps || (math.Abs(ratio-best) <= Eps && tab.basis[i] < bestBasis) {
					best = ratio
					leave = i
					bestBasis = tab.basis[i]
				}
			}
		}
		if leave == -1 {
			return StatusUnbounded, iters
		}
		tab.linprogPivot(leave, entering)
		iters++
	}
	return StatusIterations, iters
}

// Simplex solves the standard-form program by the two-phase primal simplex
// method. Phase 1 minimizes the sum of artificial variables to find a basic
// feasible solution (reporting [StatusInfeasible] if none exists); phase 2
// minimizes the true objective (reporting [StatusUnbounded] if the objective
// decreases without bound). Bland's anti-cycling rule guarantees termination.
func Simplex(std StandardLP) SimplexResult {
	m := std.NumConstraints()
	n0 := std.NumVars()
	maxIter := 50*(m+n0) + 2000

	// Copy A and b, forcing b >= 0 by row negation.
	a := make([][]float64, m)
	b := make([]float64, m)
	for i := 0; i < m; i++ {
		row := append([]float64(nil), std.A[i]...)
		bi := std.B[i]
		if bi < 0 {
			for j := range row {
				row[j] = -row[j]
			}
			bi = -bi
		}
		a[i] = row
		b[i] = bi
	}

	// Extend with m artificial columns forming an identity.
	N := n0 + m
	tab := &linprogTableau{
		t:     make([][]float64, m),
		rhs:   b,
		basis: make([]int, m),
		m:     m,
		n:     N,
	}
	forbidden := make([]bool, N)
	for i := 0; i < m; i++ {
		row := make([]float64, N)
		copy(row, a[i])
		row[n0+i] = 1
		tab.t[i] = row
		tab.basis[i] = n0 + i
		forbidden[n0+i] = true // artificials may not re-enter in phase 2
	}

	// Phase 1: minimize sum of artificials.
	phase1Cost := make([]float64, N)
	for j := n0; j < N; j++ {
		phase1Cost[j] = 1
	}
	_, it1 := tab.linprogOptimize(phase1Cost, nil, maxIter)

	// Phase-1 objective value = sum of artificial values.
	var infeas float64
	for i := 0; i < m; i++ {
		if tab.basis[i] >= n0 {
			infeas += tab.rhs[i]
		}
	}
	if infeas > 1e-7 {
		return SimplexResult{Status: StatusInfeasible, Iterations: it1}
	}

	// Drive any artificial that remains basic (at value ~0) out of the basis.
	for i := 0; i < m; i++ {
		if tab.basis[i] < n0 {
			continue
		}
		for k := 0; k < n0; k++ {
			if !tab.linprogInBasis(k) && math.Abs(tab.t[i][k]) > Eps {
				tab.linprogPivot(i, k)
				break
			}
		}
	}

	// Phase 2: minimize the real objective, artificials forbidden.
	cost := make([]float64, N)
	copy(cost, std.C)
	status, it2 := tab.linprogOptimize(cost, forbidden, maxIter)

	x := make([]float64, n0)
	for i := 0; i < m; i++ {
		if tab.basis[i] < n0 {
			x[tab.basis[i]] = tab.rhs[i]
		}
	}
	res := SimplexResult{
		Status:     status,
		X:          x,
		Objective:  std.Objective(x),
		Basis:      append([]int(nil), tab.basis...),
		Iterations: it1 + it2,
	}
	return res
}

// Solve solves the general [LP] by reducing it to standard form and running
// the two-phase [Simplex]. The returned [LPSolution] reports the original
// variables and the value of the original objective (with the program's
// [Sense] respected).
func (lp LP) Solve() LPSolution {
	std := lp.Standard()
	res := Simplex(std)
	sol := LPSolution{
		Status:     res.Status,
		Iterations: res.Iterations,
	}
	if res.Status == StatusOptimal {
		x := make([]float64, lp.NumVars())
		copy(x, res.X[:lp.NumVars()])
		sol.X = x
		sol.Objective = lp.Objective(x)
	}
	return sol
}
