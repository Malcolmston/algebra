package controltheory

// controltheoryDeriv computes x' = A x + B u for the state-space system.
func controltheoryDeriv(s StateSpace, x []float64, u float64) []float64 {
	dx := controltheoryMatVec(s.A, x)
	for i := range dx {
		dx[i] += s.B[i] * u
	}
	return dx
}

// controltheoryOutput computes y = C x + D u.
func controltheoryOutput(s StateSpace, x []float64, u float64) float64 {
	var y float64
	for i := range x {
		y += s.C[i] * x[i]
	}
	return y + s.D*u
}

// Simulate integrates the state-space system forward from initial state x0
// over the given monotonically increasing time points using the classical
// fourth-order Runge-Kutta method. The scalar input at time t is supplied by
// the function input. It returns the output y sampled at each time point.
// The returned slice has the same length as times.
func (s StateSpace) Simulate(times []float64, x0 []float64, input func(t float64) float64) []float64 {
	n := len(s.A)
	x := make([]float64, n)
	copy(x, x0)
	out := make([]float64, len(times))
	if len(times) == 0 {
		return out
	}
	out[0] = controltheoryOutput(s, x, input(times[0]))
	for k := 1; k < len(times); k++ {
		t0 := times[k-1]
		h := times[k] - t0
		x = controltheoryRK4Step(s, x, t0, h, input)
		out[k] = controltheoryOutput(s, x, input(times[k]))
	}
	return out
}

// controltheoryRK4Step advances the state by one RK4 step of size h from t0.
func controltheoryRK4Step(s StateSpace, x []float64, t0, h float64, input func(float64) float64) []float64 {
	n := len(x)
	add := func(a, b []float64, scale float64) []float64 {
		r := make([]float64, n)
		for i := range a {
			r[i] = a[i] + scale*b[i]
		}
		return r
	}
	k1 := controltheoryDeriv(s, x, input(t0))
	k2 := controltheoryDeriv(s, add(x, k1, h/2), input(t0+h/2))
	k3 := controltheoryDeriv(s, add(x, k2, h/2), input(t0+h/2))
	k4 := controltheoryDeriv(s, add(x, k3, h), input(t0+h))
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = x[i] + (h/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i])
	}
	return out
}

// StepResponse returns the unit-step response of the transfer function sampled
// at the given monotonically increasing time points. The system starts from
// rest (zero initial state) and the input is held at 1 for all t >= 0. The
// transfer function must be proper.
func (g TransferFunction) StepResponse(times []float64) []float64 {
	ss := TransferFunctionToStateSpace(g)
	x0 := make([]float64, ss.Order())
	return ss.Simulate(times, x0, func(float64) float64 { return 1 })
}

// StepResponse returns the unit-step response of the state-space system from
// rest, sampled at the given time points.
func (s StateSpace) StepResponse(times []float64) []float64 {
	x0 := make([]float64, s.Order())
	return s.Simulate(times, x0, func(float64) float64 { return 1 })
}

// ImpulseResponse returns the unit-impulse response of a strictly proper
// transfer function sampled at the given monotonically increasing time points.
// For a strictly proper realization the impulse response equals C·e^{At}·B, so
// it is computed by simulating the autonomous system with initial state B and
// zero input. The transfer function must be strictly proper.
func (g TransferFunction) ImpulseResponse(times []float64) []float64 {
	ss := TransferFunctionToStateSpace(g)
	return ss.ImpulseResponse(times)
}

// ImpulseResponse returns the unit-impulse response C·e^{At}·B of the
// state-space system sampled at the given time points, obtained by simulating
// the autonomous dynamics from initial state B. Any direct feedthrough D
// contributes only at t=0 (a Dirac impulse) and is not included in the samples.
func (s StateSpace) ImpulseResponse(times []float64) []float64 {
	x0 := append([]float64{}, s.B...)
	// Output without feedthrough on the input (u = 0).
	old := s.D
	s.D = 0
	res := s.Simulate(times, x0, func(float64) float64 { return 0 })
	s.D = old
	return res
}

// FinalValue returns the steady-state value of the step response predicted by
// the final-value theorem, lim_{t->inf} y(t) = G(0), when the closed system is
// stable. It equals the DC gain of the transfer function.
func (g TransferFunction) FinalValue() float64 {
	return g.DCGain()
}
