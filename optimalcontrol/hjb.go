package optimalcontrol

import (
	"math"
	"sort"
)

// LinearInterp1D linearly interpolates the tabulated function (grid, values) at
// x. The grid must be sorted ascending. Queries outside the grid are clamped to
// the endpoint values.
func LinearInterp1D(grid, values []float64, x float64) float64 {
	n := len(grid)
	if n == 0 {
		return 0
	}
	if x <= grid[0] {
		return values[0]
	}
	if x >= grid[n-1] {
		return values[n-1]
	}
	i := sort.SearchFloat64s(grid, x)
	// grid[i-1] < x <= grid[i]
	x0, x1 := grid[i-1], grid[i]
	y0, y1 := values[i-1], values[i]
	t := (x - x0) / (x1 - x0)
	return y0 + t*(y1-y0)
}

// HJBGrid1D describes a one-dimensional infinite-horizon optimal-control problem
// discretized for solution by value iteration (a semi-Lagrangian scheme for the
// Hamilton–Jacobi–Bellman equation ρ V = min_u [ L(x,u) + V'(x) f(x,u) ]).
type HJBGrid1D struct {
	// Grid is the ascending state grid.
	Grid []float64
	// Controls is the finite set of admissible control values searched at each
	// grid point.
	Controls []float64
	// Dynamics returns f(x, u), the state velocity.
	Dynamics func(x, u float64) float64
	// RunningCost returns L(x, u), the instantaneous cost.
	RunningCost func(x, u float64) float64
	// Dt is the time step of the semi-Lagrangian discretization.
	Dt float64
	// Rho is the (non-negative) discount rate.
	Rho float64
}

// HJBResult holds the value function and greedy control at each grid point.
type HJBResult struct {
	// Value is the optimal cost-to-go at each grid point.
	Value []float64
	// Control is the optimal control at each grid point.
	Control []float64
	// Iterations is the number of value-iteration sweeps performed.
	Iterations int
	// Converged reports whether the sweep residual fell below the tolerance.
	Converged bool
}

// Solve runs value iteration for the infinite-horizon discounted problem until
// the max-norm change is below tol or maxIter sweeps elapse. The per-step
// discount is exp(−ρ Δt) and the successor value is obtained by linear
// interpolation of the current value function at x + Δt f(x, u).
func (h *HJBGrid1D) Solve(tol float64, maxIter int) *HJBResult {
	n := len(h.Grid)
	v := make([]float64, n)
	disc := math.Exp(-h.Rho * h.Dt)
	res := &HJBResult{Value: v, Control: make([]float64, n)}
	for iter := 1; iter <= maxIter; iter++ {
		newV := make([]float64, n)
		ctrl := make([]float64, n)
		var delta float64
		for i, x := range h.Grid {
			best := math.Inf(1)
			var bestU float64
			for _, u := range h.Controls {
				xn := x + h.Dt*h.Dynamics(x, u)
				cost := h.Dt*h.RunningCost(x, u) + disc*LinearInterp1D(h.Grid, v, xn)
				if cost < best {
					best = cost
					bestU = u
				}
			}
			newV[i] = best
			ctrl[i] = bestU
			if d := math.Abs(best - v[i]); d > delta {
				delta = d
			}
		}
		v = newV
		res.Value = newV
		res.Control = ctrl
		res.Iterations = iter
		if delta < tol {
			res.Converged = true
			break
		}
	}
	return res
}

// SolveFiniteHorizon solves the finite-horizon problem with terminal cost g(x)
// by backward value iteration over the given number of stages, returning the
// value function and control at the initial time (stage 0).
func (h *HJBGrid1D) SolveFiniteHorizon(stages int, terminal func(x float64) float64) *HJBResult {
	n := len(h.Grid)
	v := make([]float64, n)
	for i, x := range h.Grid {
		v[i] = terminal(x)
	}
	disc := math.Exp(-h.Rho * h.Dt)
	ctrl := make([]float64, n)
	for k := 0; k < stages; k++ {
		newV := make([]float64, n)
		for i, x := range h.Grid {
			best := math.Inf(1)
			var bestU float64
			for _, u := range h.Controls {
				xn := x + h.Dt*h.Dynamics(x, u)
				cost := h.Dt*h.RunningCost(x, u) + disc*LinearInterp1D(h.Grid, v, xn)
				if cost < best {
					best = cost
					bestU = u
				}
			}
			newV[i] = best
			ctrl[i] = bestU
		}
		v = newV
	}
	return &HJBResult{Value: v, Control: ctrl, Iterations: stages, Converged: true}
}

// ControlAt returns the greedy control at state x by interpolating the tabulated
// policy of a solved HJBResult onto the grid.
func (h *HJBGrid1D) ControlAt(res *HJBResult, x float64) float64 {
	// Nearest grid point is adequate for a piecewise-constant control law.
	best := 0
	bestDist := math.Inf(1)
	for i, g := range h.Grid {
		if d := math.Abs(g - x); d < bestDist {
			bestDist = d
			best = i
		}
	}
	return res.Control[best]
}
