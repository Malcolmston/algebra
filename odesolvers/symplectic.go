package odesolvers

import "math"

// Acceleration is the force field of a separable second-order system
// q” = a(t, q), returning the acceleration (force per unit mass) at
// configuration q and time t. The symplectic integrators integrate the
// first-order Hamiltonian system q' = v, v' = a(t, q).
type Acceleration func(t float64, q []float64) []float64

// SymplecticSolution stores the trajectory of a second-order system produced by
// a symplectic integrator: positions Q and velocities P sampled at times T.
type SymplecticSolution struct {
	T      []float64
	Q      [][]float64
	P      [][]float64
	Method string
}

// Len returns the number of stored samples.
func (s *SymplecticSolution) Len() int { return len(s.T) }

// Dim returns the configuration-space dimension.
func (s *SymplecticSolution) Dim() int {
	if len(s.Q) == 0 {
		return 0
	}
	return len(s.Q[0])
}

// FinalPosition returns a copy of the last stored position.
func (s *SymplecticSolution) FinalPosition() []float64 { return Clone(s.Q[len(s.Q)-1]) }

// FinalVelocity returns a copy of the last stored velocity.
func (s *SymplecticSolution) FinalVelocity() []float64 { return Clone(s.P[len(s.P)-1]) }

// PhaseState returns the combined phase-space vector [q, v] at sample i.
func (s *SymplecticSolution) PhaseState(i int) []float64 {
	return append(Clone(s.Q[i]), s.P[i]...)
}

func (s *SymplecticSolution) push(t float64, q, p []float64) {
	s.T = append(s.T, t)
	s.Q = append(s.Q, Clone(q))
	s.P = append(s.P, Clone(p))
}

// VelocityVerletStep advances one velocity-Verlet (kick-drift-kick) step of the
// system q” = a(t, q). Given the position q, velocity v and the precomputed
// acceleration aCur = a(t, q), it returns the updated position, velocity and
// acceleration a(t+h, qNew).
func VelocityVerletStep(a Acceleration, t float64, q, v, aCur []float64, h float64) (qNew, vNew, aNew []float64) {
	n := len(q)
	qNew = make([]float64, n)
	for i := 0; i < n; i++ {
		qNew[i] = q[i] + h*v[i] + 0.5*h*h*aCur[i]
	}
	aNew = a(t+h, qNew)
	vNew = make([]float64, n)
	for i := 0; i < n; i++ {
		vNew[i] = v[i] + 0.5*h*(aCur[i]+aNew[i])
	}
	return qNew, vNew, aNew
}

// PositionVerletStep advances one position-Verlet (drift-kick-drift) step of
// the system q” = a(t, q) and returns the updated position and velocity.
func PositionVerletStep(a Acceleration, t float64, q, v []float64, h float64) (qNew, vNew []float64) {
	n := len(q)
	qHalf := make([]float64, n)
	for i := 0; i < n; i++ {
		qHalf[i] = q[i] + 0.5*h*v[i]
	}
	acc := a(t+0.5*h, qHalf)
	vNew = make([]float64, n)
	for i := 0; i < n; i++ {
		vNew[i] = v[i] + h*acc[i]
	}
	qNew = make([]float64, n)
	for i := 0; i < n; i++ {
		qNew[i] = qHalf[i] + 0.5*h*vNew[i]
	}
	return qNew, vNew
}

// LeapfrogStep advances one leapfrog step. It is algebraically equivalent to
// velocity Verlet and shares its signature, recomputing the current
// acceleration internally for convenience.
func LeapfrogStep(a Acceleration, t float64, q, v []float64, h float64) (qNew, vNew []float64) {
	aCur := a(t, q)
	qNew, vNew, _ = VelocityVerletStep(a, t, q, v, aCur, h)
	return qNew, vNew
}

// SolveVelocityVerlet integrates q” = a(t, q) from t0 to tEnd with a fixed step
// h using the velocity-Verlet method, starting from position q0 and velocity
// v0.
func SolveVelocityVerlet(a Acceleration, t0 float64, q0, v0 []float64, tEnd, h float64) *SymplecticSolution {
	sol := &SymplecticSolution{Method: "Velocity Verlet"}
	nSteps, step := stepCount(t0, tEnd, h)
	q, v := Clone(q0), Clone(v0)
	t := t0
	sol.push(t, q, v)
	acc := a(t, q)
	for i := 0; i < nSteps; i++ {
		q, v, acc = VelocityVerletStep(a, t, q, v, acc, step)
		t = t0 + float64(i+1)*step
		sol.push(t, q, v)
	}
	return sol
}

// SolvePositionVerlet integrates q” = a(t, q) with the position-Verlet method
// and a fixed step h.
func SolvePositionVerlet(a Acceleration, t0 float64, q0, v0 []float64, tEnd, h float64) *SymplecticSolution {
	sol := &SymplecticSolution{Method: "Position Verlet"}
	nSteps, step := stepCount(t0, tEnd, h)
	q, v := Clone(q0), Clone(v0)
	t := t0
	sol.push(t, q, v)
	for i := 0; i < nSteps; i++ {
		q, v = PositionVerletStep(a, t, q, v, step)
		t = t0 + float64(i+1)*step
		sol.push(t, q, v)
	}
	return sol
}

// SolveLeapfrog integrates q” = a(t, q) with the leapfrog method and a fixed
// step h.
func SolveLeapfrog(a Acceleration, t0 float64, q0, v0 []float64, tEnd, h float64) *SymplecticSolution {
	sol := &SymplecticSolution{Method: "Leapfrog"}
	nSteps, step := stepCount(t0, tEnd, h)
	q, v := Clone(q0), Clone(v0)
	t := t0
	sol.push(t, q, v)
	acc := a(t, q)
	for i := 0; i < nSteps; i++ {
		q, v, acc = VelocityVerletStep(a, t, q, v, acc, step)
		t = t0 + float64(i+1)*step
		sol.push(t, q, v)
	}
	return sol
}

// YoshidaCoefficients returns the three composition weights (w1, w0, w1) of the
// fourth-order Yoshida method built from a symmetric second-order base step.
func YoshidaCoefficients() (w1, w0 float64) {
	cbrt2 := math.Cbrt(2)
	w1 = 1.0 / (2.0 - cbrt2)
	w0 = -cbrt2 / (2.0 - cbrt2)
	return w1, w0
}

// Yoshida4Step advances one fourth-order Yoshida step by composing three
// velocity-Verlet substeps of sizes w1*h, w0*h and w1*h. It returns the updated
// position and velocity.
func Yoshida4Step(a Acceleration, t float64, q, v []float64, h float64) (qNew, vNew []float64) {
	w1, w0 := YoshidaCoefficients()
	q1, v1 := Clone(q), Clone(v)
	// Substep 1.
	acc := a(t, q1)
	q1, v1, _ = VelocityVerletStep(a, t, q1, v1, acc, w1*h)
	t += w1 * h
	// Substep 2.
	acc = a(t, q1)
	q1, v1, _ = VelocityVerletStep(a, t, q1, v1, acc, w0*h)
	t += w0 * h
	// Substep 3.
	acc = a(t, q1)
	q1, v1, _ = VelocityVerletStep(a, t, q1, v1, acc, w1*h)
	return q1, v1
}

// SolveYoshida4 integrates q” = a(t, q) with the fourth-order Yoshida method
// and a fixed step h.
func SolveYoshida4(a Acceleration, t0 float64, q0, v0 []float64, tEnd, h float64) *SymplecticSolution {
	sol := &SymplecticSolution{Method: "Yoshida 4"}
	nSteps, step := stepCount(t0, tEnd, h)
	q, v := Clone(q0), Clone(v0)
	t := t0
	sol.push(t, q, v)
	for i := 0; i < nSteps; i++ {
		q, v = Yoshida4Step(a, t, q, v, step)
		t = t0 + float64(i+1)*step
		sol.push(t, q, v)
	}
	return sol
}

// HarmonicEnergy returns the total energy 0.5*|v|^2 + 0.5*omega^2*|q|^2 of a
// unit-mass isotropic harmonic oscillator, a convenient conserved quantity for
// checking symplectic integrators.
func HarmonicEnergy(q, v []float64, omega float64) float64 {
	return 0.5*Dot(v, v) + 0.5*omega*omega*Dot(q, q)
}
