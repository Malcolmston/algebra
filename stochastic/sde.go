package stochastic

import (
	"math"
	"math/rand"
)

// DriftFunc is the drift coefficient a(t, x) of a scalar Ito SDE
// dX = a(t, X) dt + b(t, X) dW.
type DriftFunc func(t, x float64) float64

// DiffusionFunc is the diffusion coefficient b(t, x) of a scalar Ito SDE.
type DiffusionFunc func(t, x float64) float64

// SDE bundles the drift and diffusion coefficients of a scalar Ito stochastic
// differential equation dX = Drift(t, X) dt + Diffusion(t, X) dW.
type SDE struct {
	Drift     DriftFunc
	Diffusion DiffusionFunc
}

// EulerMaruyama integrates the SDE on [0, T] starting at x0 with the
// Euler-Maruyama scheme over nSteps steps and returns the sample path.
func EulerMaruyama(rng *rand.Rand, sde SDE, x0, T float64, nSteps int) Path {
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	sd := math.Sqrt(dt)
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	x := x0
	vals[0] = x0
	for i := 1; i <= nSteps; i++ {
		t := float64(i-1) * dt
		dW := sd * rng.NormFloat64()
		x += sde.Drift(t, x)*dt + sde.Diffusion(t, x)*dW
		times[i] = float64(i) * dt
		vals[i] = x
	}
	return Path{Times: times, Values: vals}
}

// EulerMaruyamaStep advances one Euler-Maruyama step from state x at time t over
// a step of size dt using the supplied Brownian increment dW.
func EulerMaruyamaStep(sde SDE, t, x, dt, dW float64) float64 {
	return x + sde.Drift(t, x)*dt + sde.Diffusion(t, x)*dW
}

// Milstein integrates the SDE with the Milstein scheme, which adds the
// first-order Ito correction 0.5*b*b'*(dW^2 - dt). The derivative of the
// diffusion coefficient with respect to x is supplied as dDiffusion.
func Milstein(rng *rand.Rand, sde SDE, dDiffusion DiffusionFunc, x0, T float64, nSteps int) Path {
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	sd := math.Sqrt(dt)
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	x := x0
	vals[0] = x0
	for i := 1; i <= nSteps; i++ {
		t := float64(i-1) * dt
		dW := sd * rng.NormFloat64()
		b := sde.Diffusion(t, x)
		x += sde.Drift(t, x)*dt + b*dW + 0.5*b*dDiffusion(t, x)*(dW*dW-dt)
		times[i] = float64(i) * dt
		vals[i] = x
	}
	return Path{Times: times, Values: vals}
}

// MilsteinNumeric is like Milstein but approximates the derivative of the
// diffusion coefficient by a central finite difference with step h.
func MilsteinNumeric(rng *rand.Rand, sde SDE, x0, T float64, nSteps int, h float64) Path {
	if h <= 0 {
		h = 1e-6
	}
	d := func(t, x float64) float64 {
		return (sde.Diffusion(t, x+h) - sde.Diffusion(t, x-h)) / (2 * h)
	}
	return Milstein(rng, sde, d, x0, T, nSteps)
}

// StratonovichHeun integrates the Stratonovich SDE dX = a dt + b o dW with the
// stochastic Heun (predictor-corrector) scheme, which converges to the
// Stratonovich solution.
func StratonovichHeun(rng *rand.Rand, sde SDE, x0, T float64, nSteps int) Path {
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	sd := math.Sqrt(dt)
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	x := x0
	vals[0] = x0
	for i := 1; i <= nSteps; i++ {
		t := float64(i-1) * dt
		tn := float64(i) * dt
		dW := sd * rng.NormFloat64()
		a0 := sde.Drift(t, x)
		b0 := sde.Diffusion(t, x)
		xPred := x + a0*dt + b0*dW
		a1 := sde.Drift(tn, xPred)
		b1 := sde.Diffusion(tn, xPred)
		x += 0.5*(a0+a1)*dt + 0.5*(b0+b1)*dW
		times[i] = tn
		vals[i] = x
	}
	return Path{Times: times, Values: vals}
}

// TamedEulerMaruyama integrates the SDE with the tamed Euler scheme, which
// bounds the drift term to remain stable for super-linearly growing drift.
func TamedEulerMaruyama(rng *rand.Rand, sde SDE, x0, T float64, nSteps int) Path {
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	sd := math.Sqrt(dt)
	times := make([]float64, nSteps+1)
	vals := make([]float64, nSteps+1)
	x := x0
	vals[0] = x0
	for i := 1; i <= nSteps; i++ {
		t := float64(i-1) * dt
		dW := sd * rng.NormFloat64()
		a := sde.Drift(t, x)
		tamed := a * dt / (1 + math.Abs(a)*dt)
		x += tamed + sde.Diffusion(t, x)*dW
		times[i] = float64(i) * dt
		vals[i] = x
	}
	return Path{Times: times, Values: vals}
}

// EulerMaruyamaSystem integrates a d-dimensional Ito SDE system with diagonal
// noise, dX_k = drift_k(t, X) dt + diffusion_k(t, X) dW_k, using independent
// Brownian motions per component. It returns one Path per component.
func EulerMaruyamaSystem(rng *rand.Rand, drift, diffusion []func(t float64, x []float64) float64, x0 []float64, T float64, nSteps int) []Path {
	d := len(x0)
	if nSteps < 1 {
		nSteps = 1
	}
	dt := T / float64(nSteps)
	sd := math.Sqrt(dt)
	paths := make([]Path, d)
	for k := 0; k < d; k++ {
		paths[k].Times = make([]float64, nSteps+1)
		paths[k].Values = make([]float64, nSteps+1)
		paths[k].Values[0] = x0[k]
	}
	x := append([]float64(nil), x0...)
	xn := make([]float64, d)
	for i := 1; i <= nSteps; i++ {
		t := float64(i-1) * dt
		for k := 0; k < d; k++ {
			dW := sd * rng.NormFloat64()
			xn[k] = x[k] + drift[k](t, x)*dt + diffusion[k](t, x)*dW
		}
		copy(x, xn)
		for k := 0; k < d; k++ {
			paths[k].Times[i] = float64(i) * dt
			paths[k].Values[i] = x[k]
		}
	}
	return paths
}

// SimulateEnsemble runs EulerMaruyama for the given SDE paths times, each with an
// independent stream derived from the base seed, and returns the paths.
func SimulateEnsemble(seed int64, sde SDE, x0, T float64, nSteps, paths int) []Path {
	out := make([]Path, 0, paths)
	for p := 0; p < paths; p++ {
		rng := NewRNG(seed + int64(p))
		out = append(out, EulerMaruyama(rng, sde, x0, T, nSteps))
	}
	return out
}

// MonteCarloExpectation estimates E[payoff(X)] where X is the terminal-value
// Path of the SDE, by averaging over the given number of independent paths. It
// returns the sample mean and the standard error of that mean.
func MonteCarloExpectation(seed int64, sde SDE, x0, T float64, nSteps, paths int, payoff func(Path) float64) (mean, stderr float64) {
	if paths <= 0 {
		return 0, 0
	}
	sum := 0.0
	sumSq := 0.0
	for p := 0; p < paths; p++ {
		rng := NewRNG(seed + int64(p))
		v := payoff(EulerMaruyama(rng, sde, x0, T, nSteps))
		sum += v
		sumSq += v * v
	}
	n := float64(paths)
	mean = sum / n
	if paths < 2 {
		return mean, 0
	}
	variance := (sumSq - n*mean*mean) / (n - 1)
	if variance < 0 {
		variance = 0
	}
	return mean, math.Sqrt(variance / n)
}
