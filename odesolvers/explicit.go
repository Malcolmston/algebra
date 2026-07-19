package odesolvers

import "math"

// stepCount returns the number of uniform steps of nominal size h needed to
// advance from t0 to tEnd, together with the signed step actually used so that
// an integral number of steps lands exactly on tEnd.
func stepCount(t0, tEnd, h float64) (int, float64) {
	span := tEnd - t0
	if span == 0 || h == 0 {
		return 0, 0
	}
	h = math.Abs(h)
	n := int(math.Ceil(math.Abs(span)/h - 1e-12))
	if n < 1 {
		n = 1
	}
	return n, span / float64(n)
}

// IntegrateFixed integrates y' = f(t, y) from t0 to tEnd using the explicit
// method described by bt with a fixed step of nominal size h. The step size is
// adjusted slightly so that a whole number of steps lands exactly on tEnd, and
// integration proceeds backward in time when tEnd < t0. Every step is recorded
// in the returned Solution.
func IntegrateFixed(f Field, bt *ButcherTableau, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	sol := newSolution(bt.Name, t0, y0)
	n, step := stepCount(t0, tEnd, h)
	y := Clone(y0)
	t := t0
	for i := 0; i < n; i++ {
		y = bt.Step(f, t, y, step)
		t = t0 + float64(i+1)*step
		sol.push(t, y)
	}
	return sol
}

// IntegrateFixedDense behaves like [IntegrateFixed] but also records the stage-0
// derivative at every sample so that the returned Solution supports cubic
// Hermite dense output.
func IntegrateFixedDense(f Field, bt *ButcherTableau, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	sol := &Solution{Method: bt.Name}
	n, step := stepCount(t0, tEnd, h)
	y := Clone(y0)
	t := t0
	sol.pushWithDeriv(t, y, f(t, y))
	for i := 0; i < n; i++ {
		y = bt.Step(f, t, y, step)
		t = t0 + float64(i+1)*step
		sol.pushWithDeriv(t, y, f(t, y))
	}
	return sol
}

// EulerStep advances one forward-Euler step and returns the new state.
func EulerStep(f Field, t float64, y []float64, h float64) []float64 {
	return AXPY(h, f(t, y), y)
}

// MidpointStep advances one explicit-midpoint step and returns the new state.
func MidpointStep(f Field, t float64, y []float64, h float64) []float64 {
	k1 := f(t, y)
	k2 := f(t+0.5*h, AXPY(0.5*h, k1, y))
	return AXPY(h, k2, y)
}

// HeunStep advances one step of Heun's method and returns the new state.
func HeunStep(f Field, t float64, y []float64, h float64) []float64 {
	k1 := f(t, y)
	k2 := f(t+h, AXPY(h, k1, y))
	out := Clone(y)
	for i := range out {
		out[i] += 0.5 * h * (k1[i] + k2[i])
	}
	return out
}

// RalstonStep advances one step of Ralston's method and returns the new state.
func RalstonStep(f Field, t float64, y []float64, h float64) []float64 {
	return RalstonTableau().Step(f, t, y, h)
}

// RK4Step advances one classical fourth-order Runge-Kutta step and returns the
// new state.
func RK4Step(f Field, t float64, y []float64, h float64) []float64 {
	k1 := f(t, y)
	k2 := f(t+0.5*h, AXPY(0.5*h, k1, y))
	k3 := f(t+0.5*h, AXPY(0.5*h, k2, y))
	k4 := f(t+h, AXPY(h, k3, y))
	out := Clone(y)
	for i := range out {
		out[i] += h / 6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
	}
	return out
}

// RK38Step advances one step of the fourth-order 3/8-rule and returns the new
// state.
func RK38Step(f Field, t float64, y []float64, h float64) []float64 {
	return RK38Tableau().Step(f, t, y, h)
}

// SolveEuler integrates with the forward Euler method and a fixed step h.
func SolveEuler(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	return IntegrateFixed(f, EulerTableau(), t0, y0, tEnd, h)
}

// SolveMidpoint integrates with the explicit midpoint method and a fixed step h.
func SolveMidpoint(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	return IntegrateFixed(f, MidpointTableau(), t0, y0, tEnd, h)
}

// SolveHeun integrates with Heun's method and a fixed step h.
func SolveHeun(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	return IntegrateFixed(f, HeunTableau(), t0, y0, tEnd, h)
}

// SolveRalston integrates with Ralston's method and a fixed step h.
func SolveRalston(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	return IntegrateFixed(f, RalstonTableau(), t0, y0, tEnd, h)
}

// SolveSSPRK3 integrates with the SSPRK3 method and a fixed step h.
func SolveSSPRK3(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	return IntegrateFixed(f, SSPRK3Tableau(), t0, y0, tEnd, h)
}

// SolveRK4 integrates with the classical fourth-order Runge-Kutta method and a
// fixed step h.
func SolveRK4(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	return IntegrateFixed(f, RK4Tableau(), t0, y0, tEnd, h)
}

// SolveRK38 integrates with the fourth-order 3/8-rule and a fixed step h.
func SolveRK38(f Field, t0 float64, y0 []float64, tEnd, h float64) *Solution {
	return IntegrateFixed(f, RK38Tableau(), t0, y0, tEnd, h)
}

// rk4FinalState integrates with RK4 and returns only the final state, used
// internally as a robust bootstrap for the multistep and shooting routines.
func rk4FinalState(f Field, t0 float64, y0 []float64, tEnd float64, steps int) []float64 {
	if steps < 1 {
		steps = 1
	}
	h := (tEnd - t0) / float64(steps)
	y := Clone(y0)
	t := t0
	for i := 0; i < steps; i++ {
		y = RK4Step(f, t, y, h)
		t = t0 + float64(i+1)*h
	}
	return y
}
