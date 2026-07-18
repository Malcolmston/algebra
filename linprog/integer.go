package linprog

import "math"

// IntTol is the tolerance within which a relaxed variable value is treated as
// integral by the branch-and-bound solver.
const IntTol = 1e-6

// MaxBranchNodes caps the number of branch-and-bound nodes explored before the
// solver gives up and returns [StatusIterations] with the best incumbent found
// so far.
const MaxBranchNodes = 200000

// SolveInteger solves the mixed-integer linear program obtained from lp by
// requiring the variables flagged in integer to take integer values. It uses
// branch and bound on the LP relaxation, branching on the first fractional
// integer-constrained variable and pruning nodes that cannot improve on the
// incumbent.
//
// integer must have length NumVars; a true entry constrains that variable to
// integers. The returned [LPSolution] reports the best integer-feasible point,
// or [StatusInfeasible] if none exists.
func (lp LP) SolveInteger(integer []bool) LPSolution {
	if len(integer) != lp.NumVars() {
		panic(ErrDimension)
	}
	best := LPSolution{Status: StatusInfeasible}
	haveBest := false
	nodes := 0
	hitLimit := false

	var branch func(cur LP)
	branch = func(cur LP) {
		if hitLimit {
			return
		}
		if nodes >= MaxBranchNodes {
			hitLimit = true
			return
		}
		nodes++

		sol := cur.Solve()
		if sol.Status != StatusOptimal {
			return // infeasible or unbounded relaxation prunes this node
		}
		// Bound: relaxation cannot beat the incumbent.
		if haveBest {
			if lp.Sense == Minimize && sol.Objective >= best.Objective-IntTol {
				return
			}
			if lp.Sense == Maximize && sol.Objective <= best.Objective+IntTol {
				return
			}
		}
		// Find the first fractional integer variable.
		frac := -1
		for j := range integer {
			if !integer[j] {
				continue
			}
			f := sol.X[j] - math.Floor(sol.X[j])
			if math.Min(f, 1-f) > IntTol {
				frac = j
				break
			}
		}
		if frac == -1 {
			// Integer feasible: update incumbent.
			better := !haveBest ||
				(lp.Sense == Minimize && sol.Objective < best.Objective) ||
				(lp.Sense == Maximize && sol.Objective > best.Objective)
			if better {
				// Round near-integer components to clean integers.
				for j := range integer {
					if integer[j] {
						sol.X[j] = math.Round(sol.X[j])
					}
				}
				sol.Objective = lp.Objective(sol.X)
				best = sol
				haveBest = true
			}
			return
		}
		val := sol.X[frac]
		row := make([]float64, lp.NumVars())
		row[frac] = 1
		// Branch down: x[frac] <= floor(val).
		branch(cur.AddConstraint(row, LessEqual, math.Floor(val)))
		// Branch up: x[frac] >= ceil(val).
		branch(cur.AddConstraint(row, GreaterEqual, math.Ceil(val)))
	}

	branch(lp)
	if hitLimit && haveBest {
		best.Status = StatusIterations
	}
	if hitLimit && !haveBest {
		return LPSolution{Status: StatusIterations}
	}
	return best
}

// SolveBinary solves the integer program in which the variables flagged in
// binary are additionally restricted to {0,1} by appending x[j] <= 1 upper
// bounds before branch and bound. Variables with a false flag are left
// continuous and nonnegative.
func (lp LP) SolveBinary(binary []bool) LPSolution {
	if len(binary) != lp.NumVars() {
		panic(ErrDimension)
	}
	cur := lp
	for j := range binary {
		if binary[j] {
			row := make([]float64, lp.NumVars())
			row[j] = 1
			cur = cur.AddConstraint(row, LessEqual, 1)
		}
	}
	return cur.SolveInteger(binary)
}

// IsIntegral reports whether every entry of x is within tol of an integer.
func IsIntegral(x []float64, tol float64) bool {
	for _, v := range x {
		if math.Abs(v-math.Round(v)) > tol {
			return false
		}
	}
	return true
}
