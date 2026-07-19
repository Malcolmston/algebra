package odesolvers

import "math"

// ImplicitRKStep advances one step of a general (possibly implicit) Runge-Kutta
// method from (t, y) by h and returns the new state. The stage system
//
//	K_i = f(t + C[i]*h, y + h * sum_j A[i][j] K_j)
//
// is solved for the stage derivatives K_i by Newton's method (finite-difference
// Jacobian). It works for any tableau, including the explicit ones, but is
// intended for the implicit / stiff methods. A non-nil error indicates the
// Newton iteration failed to converge.
func ImplicitRKStep(f Field, bt *ButcherTableau, t float64, y []float64, h float64) ([]float64, error) {
	s := bt.Stages()
	n := len(y)

	// Unknown vector: K flattened as [k_0..., k_1..., ...] of length s*n.
	// Residual R_i = K_i - f(t + c_i h, y + h sum_j a_ij K_j).
	residual := func(flat []float64) []float64 {
		out := make([]float64, s*n)
		for i := 0; i < s; i++ {
			yi := Clone(y)
			for j := 0; j < s; j++ {
				a := bt.A[i][j]
				if a == 0 {
					continue
				}
				kj := flat[j*n : j*n+n]
				for m := 0; m < n; m++ {
					yi[m] += h * a * kj[m]
				}
			}
			fi := f(t+bt.C[i]*h, yi)
			ki := flat[i*n : i*n+n]
			for m := 0; m < n; m++ {
				out[i*n+m] = ki[m] - fi[m]
			}
		}
		return out
	}

	// Initial guess: all stages equal to f(t, y).
	f0 := f(t, y)
	guess := make([]float64, s*n)
	for i := 0; i < s; i++ {
		copy(guess[i*n:i*n+n], f0)
	}
	sol, err := NewtonSolve(residual, guess, 1e-11, 60)
	if err != nil {
		return nil, err
	}
	out := Clone(y)
	for i := 0; i < s; i++ {
		b := bt.B[i]
		if b == 0 {
			continue
		}
		ki := sol[i*n : i*n+n]
		for m := 0; m < n; m++ {
			out[m] += h * b * ki[m]
		}
	}
	return out, nil
}

// IntegrateImplicit integrates y' = f(t, y) from t0 to tEnd with a fixed step h
// using the implicit Runge-Kutta method bt, solving each stage system with
// Newton's method. A non-nil error aborts integration, returning the partial
// Solution reached so far.
func IntegrateImplicit(f Field, bt *ButcherTableau, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	sol := newSolution(bt.Name, t0, y0)
	nSteps, step := stepCount(t0, tEnd, h)
	y := Clone(y0)
	t := t0
	for i := 0; i < nSteps; i++ {
		ynext, err := ImplicitRKStep(f, bt, t, y, step)
		if err != nil {
			return sol, err
		}
		y = ynext
		t = t0 + float64(i+1)*step
		sol.push(t, y)
	}
	return sol, nil
}

// --- Implicit tableau constructors -----------------------------------------

// BackwardEulerTableau returns the tableau of the backward (implicit) Euler
// method, order 1, L-stable.
func BackwardEulerTableau() *ButcherTableau {
	return &ButcherTableau{
		A:     [][]float64{{1}},
		B:     []float64{1},
		C:     []float64{1},
		Order: 1,
		Name:  "Backward Euler",
	}
}

// TrapezoidalTableau returns the tableau of the implicit trapezoidal rule
// (Crank-Nicolson), order 2, A-stable.
func TrapezoidalTableau() *ButcherTableau {
	return &ButcherTableau{
		A:     [][]float64{{0, 0}, {0.5, 0.5}},
		B:     []float64{0.5, 0.5},
		C:     []float64{0, 1},
		Order: 2,
		Name:  "Trapezoidal",
	}
}

// ImplicitMidpointTableau returns the tableau of the implicit midpoint method
// (the one-stage Gauss method), order 2, A-stable and symplectic.
func ImplicitMidpointTableau() *ButcherTableau {
	return &ButcherTableau{
		A:     [][]float64{{0.5}},
		B:     []float64{1},
		C:     []float64{0.5},
		Order: 2,
		Name:  "Implicit Midpoint",
	}
}

// GaussLegendre4Tableau returns the two-stage, fourth-order Gauss-Legendre
// tableau, A-stable and symplectic.
func GaussLegendre4Tableau() *ButcherTableau {
	r3 := math.Sqrt(3)
	return &ButcherTableau{
		A: [][]float64{
			{0.25, 0.25 - r3/6},
			{0.25 + r3/6, 0.25},
		},
		B:     []float64{0.5, 0.5},
		C:     []float64{0.5 - r3/6, 0.5 + r3/6},
		Order: 4,
		Name:  "Gauss-Legendre 4",
	}
}

// RadauIIA3Tableau returns the two-stage Radau IIA tableau of order 3,
// L-stable and stiffly accurate.
func RadauIIA3Tableau() *ButcherTableau {
	return &ButcherTableau{
		A: [][]float64{
			{5.0 / 12.0, -1.0 / 12.0},
			{3.0 / 4.0, 1.0 / 4.0},
		},
		B:     []float64{3.0 / 4.0, 1.0 / 4.0},
		C:     []float64{1.0 / 3.0, 1},
		Order: 3,
		Name:  "Radau IIA (order 3)",
	}
}

// RadauIIA5Tableau returns the three-stage Radau IIA tableau of order 5,
// L-stable and stiffly accurate, the workhorse implicit method for stiff
// systems.
func RadauIIA5Tableau() *ButcherTableau {
	s6 := math.Sqrt(6)
	return &ButcherTableau{
		A: [][]float64{
			{(88 - 7*s6) / 360, (296 - 169*s6) / 1800, (-2 + 3*s6) / 225},
			{(296 + 169*s6) / 1800, (88 + 7*s6) / 360, (-2 - 3*s6) / 225},
			{(16 - s6) / 36, (16 + s6) / 36, 1.0 / 9.0},
		},
		B:     []float64{(16 - s6) / 36, (16 + s6) / 36, 1.0 / 9.0},
		C:     []float64{(4 - s6) / 10, (4 + s6) / 10, 1},
		Order: 5,
		Name:  "Radau IIA (order 5)",
	}
}

// --- Single-step implicit convenience functions ----------------------------

// BackwardEulerStep advances one backward-Euler step and returns the new state.
func BackwardEulerStep(f Field, t float64, y []float64, h float64) ([]float64, error) {
	return ImplicitRKStep(f, BackwardEulerTableau(), t, y, h)
}

// TrapezoidalStep advances one implicit-trapezoidal step and returns the new
// state.
func TrapezoidalStep(f Field, t float64, y []float64, h float64) ([]float64, error) {
	return ImplicitRKStep(f, TrapezoidalTableau(), t, y, h)
}

// ImplicitMidpointStep advances one implicit-midpoint step and returns the new
// state.
func ImplicitMidpointStep(f Field, t float64, y []float64, h float64) ([]float64, error) {
	return ImplicitRKStep(f, ImplicitMidpointTableau(), t, y, h)
}

// RadauIIA5Step advances one order-5 Radau IIA step and returns the new state.
func RadauIIA5Step(f Field, t float64, y []float64, h float64) ([]float64, error) {
	return ImplicitRKStep(f, RadauIIA5Tableau(), t, y, h)
}

// --- Full-interval implicit convenience solvers ----------------------------

// SolveBackwardEuler integrates with the backward Euler method and a fixed
// step h.
func SolveBackwardEuler(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	return IntegrateImplicit(f, BackwardEulerTableau(), t0, y0, tEnd, h)
}

// SolveTrapezoidal integrates with the implicit trapezoidal rule and a fixed
// step h.
func SolveTrapezoidal(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	return IntegrateImplicit(f, TrapezoidalTableau(), t0, y0, tEnd, h)
}

// SolveImplicitMidpoint integrates with the implicit midpoint method and a
// fixed step h.
func SolveImplicitMidpoint(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	return IntegrateImplicit(f, ImplicitMidpointTableau(), t0, y0, tEnd, h)
}

// SolveGaussLegendre4 integrates with the fourth-order Gauss-Legendre method
// and a fixed step h.
func SolveGaussLegendre4(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	return IntegrateImplicit(f, GaussLegendre4Tableau(), t0, y0, tEnd, h)
}

// SolveRadauIIA3 integrates with the order-3 Radau IIA method and a fixed
// step h.
func SolveRadauIIA3(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	return IntegrateImplicit(f, RadauIIA3Tableau(), t0, y0, tEnd, h)
}

// SolveRadauIIA5 integrates with the order-5 Radau IIA method and a fixed
// step h.
func SolveRadauIIA5(f Field, t0 float64, y0 []float64, tEnd, h float64) (*Solution, error) {
	return IntegrateImplicit(f, RadauIIA5Tableau(), t0, y0, tEnd, h)
}
