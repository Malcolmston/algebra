package odesolvers

// BoundaryFunc encodes the boundary conditions of a two-point boundary-value
// problem. Given the state ya at the left endpoint and yb at the right
// endpoint, it returns a residual vector that must vanish at a solution. For a
// well-posed problem the residual has the same length as the state.
type BoundaryFunc func(ya, yb []float64) []float64

// ShootingOptions configures the shooting solvers.
type ShootingOptions struct {
	// SubSteps is the number of RK4 substeps used to propagate the state across
	// each integration interval; larger values increase propagation accuracy.
	SubSteps int
	// Tol is the Newton convergence tolerance on the residual, and MaxIter the
	// Newton iteration budget.
	Tol     float64
	MaxIter int
}

// DefaultShootingOptions returns SubSteps=200, Tol=1e-8, MaxIter=50.
func DefaultShootingOptions() ShootingOptions {
	return ShootingOptions{SubSteps: 200, Tol: 1e-8, MaxIter: 50}
}

func (o ShootingOptions) withDefaults() ShootingOptions {
	if o.SubSteps <= 0 {
		o.SubSteps = 200
	}
	if o.Tol <= 0 {
		o.Tol = 1e-8
	}
	if o.MaxIter <= 0 {
		o.MaxIter = 50
	}
	return o
}

// ShootingResidual returns, for a candidate initial state s at t0, the boundary
// residual bc(s, y(tEnd; s)) obtained by propagating s to tEnd with RK4 using
// the given number of substeps. It is the function whose root single shooting
// finds.
func ShootingResidual(f Field, bc BoundaryFunc, t0 float64, s []float64, tEnd float64, subSteps int) []float64 {
	yb := rk4FinalState(f, t0, s, tEnd, subSteps)
	return bc(s, yb)
}

// SingleShooting solves the two-point boundary-value problem y' = f(t, y) on
// [t0, tEnd] with boundary conditions bc, using single shooting: it applies
// Newton's method to find an initial state s such that bc(s, y(tEnd; s)) = 0,
// starting from guess. On success it returns the converged initial state
// together with a fine RK4 Solution of the trajectory. A non-nil error reports
// a failure of the Newton iteration.
func SingleShooting(f Field, bc BoundaryFunc, t0 float64, tEnd float64, guess []float64, opts ShootingOptions) ([]float64, *Solution, error) {
	o := opts.withDefaults()
	residual := func(s []float64) []float64 {
		return ShootingResidual(f, bc, t0, s, tEnd, o.SubSteps)
	}
	s0, err := NewtonSolve(residual, guess, o.Tol, o.MaxIter)
	if err != nil {
		return s0, nil, err
	}
	sol := IntegrateFixedDense(f, RK4Tableau(), t0, s0, tEnd, (tEnd-t0)/float64(o.SubSteps))
	return s0, sol, nil
}

// MultipleShooting solves the boundary-value problem with multiple shooting.
// The interval [t0, tEnd] is split into nIntervals subintervals at equally
// spaced nodes. The unknowns are the states at every node; interior continuity
// (the propagated state from one node must equal the next node's state) plus
// the boundary conditions bc are enforced simultaneously by Newton's method.
// guess is the initial state estimate at t0, which is propagated to seed all
// node values. It returns the node times, the converged node states and a fine
// Solution, or a non-nil error if Newton fails.
func MultipleShooting(f Field, bc BoundaryFunc, t0, tEnd float64, nIntervals int, guess []float64, opts ShootingOptions) ([]float64, [][]float64, *Solution, error) {
	o := opts.withDefaults()
	if nIntervals < 1 {
		nIntervals = 1
	}
	n := len(guess)
	nodes := Linspace(t0, tEnd, nIntervals+1)
	subPer := o.SubSteps / nIntervals
	if subPer < 1 {
		subPer = 1
	}

	// Seed node states by propagating the guess forward.
	seed := make([][]float64, nIntervals+1)
	seed[0] = Clone(guess)
	for k := 0; k < nIntervals; k++ {
		seed[k+1] = rk4FinalState(f, nodes[k], seed[k], nodes[k+1], subPer)
	}
	// Pack node states into a single unknown vector.
	pack := func(states [][]float64) []float64 {
		out := make([]float64, (nIntervals+1)*n)
		for k := 0; k <= nIntervals; k++ {
			copy(out[k*n:k*n+n], states[k])
		}
		return out
	}
	unpack := func(flat []float64) [][]float64 {
		states := make([][]float64, nIntervals+1)
		for k := 0; k <= nIntervals; k++ {
			states[k] = flat[k*n : k*n+n]
		}
		return states
	}

	residual := func(flat []float64) []float64 {
		states := unpack(flat)
		out := make([]float64, (nIntervals+1)*n)
		// Continuity residuals for each interval.
		for k := 0; k < nIntervals; k++ {
			prop := rk4FinalState(f, nodes[k], states[k], nodes[k+1], subPer)
			for m := 0; m < n; m++ {
				out[k*n+m] = prop[m] - states[k+1][m]
			}
		}
		// Boundary residuals in the final block.
		br := bc(states[0], states[nIntervals])
		for m := 0; m < n; m++ {
			out[nIntervals*n+m] = br[m]
		}
		return out
	}

	solFlat, err := NewtonSolve(residual, pack(seed), o.Tol, o.MaxIter)
	states := unpack(solFlat)
	// Copy node states out (unpack aliases solFlat storage).
	nodeStates := make([][]float64, nIntervals+1)
	for k := range states {
		nodeStates[k] = Clone(states[k])
	}
	if err != nil {
		return nodes, nodeStates, nil, err
	}
	sol := IntegrateFixedDense(f, RK4Tableau(), t0, nodeStates[0], tEnd, (tEnd-t0)/float64(o.SubSteps))
	return nodes, nodeStates, sol, nil
}
