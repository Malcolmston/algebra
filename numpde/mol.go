package numpde

// ODEFunc is the right-hand side of an autonomous or non-autonomous system of
// ordinary differential equations, y' = f(t, y). Implementations must return a
// freshly allocated slice the same length as y and must not modify y.
type ODEFunc func(t float64, y []float64) []float64

// Stepper advances an ODE system by one step of size dt from state y at time t
// using the right-hand side f, returning the new state. EulerStep, RK2Step and
// RK4Step all satisfy this signature.
type Stepper func(f ODEFunc, t float64, y []float64, dt float64) []float64

// numpdeAxpy returns y + s*x for equal-length slices.
func numpdeAxpy(y []float64, s float64, x []float64) []float64 {
	out := make([]float64, len(y))
	for i := range y {
		out[i] = y[i] + s*x[i]
	}
	return out
}

// EulerStep advances the ODE system y' = f(t, y) by one explicit (forward)
// Euler step of size dt and returns the new state. It is first-order accurate.
func EulerStep(f ODEFunc, t float64, y []float64, dt float64) []float64 {
	k := f(t, y)
	return numpdeAxpy(y, dt, k)
}

// RK2Step advances the ODE system by one step of the explicit midpoint method
// (second-order Runge-Kutta) of size dt and returns the new state.
func RK2Step(f ODEFunc, t float64, y []float64, dt float64) []float64 {
	k1 := f(t, y)
	mid := numpdeAxpy(y, dt/2, k1)
	k2 := f(t+dt/2, mid)
	return numpdeAxpy(y, dt, k2)
}

// RK4Step advances the ODE system by one step of the classical fourth-order
// Runge-Kutta method of size dt and returns the new state. This is the default
// time integrator for the method-of-lines helpers.
func RK4Step(f ODEFunc, t float64, y []float64, dt float64) []float64 {
	n := len(y)
	k1 := f(t, y)
	k2 := f(t+dt/2, numpdeAxpy(y, dt/2, k1))
	k3 := f(t+dt/2, numpdeAxpy(y, dt/2, k2))
	k4 := f(t+dt, numpdeAxpy(y, dt, k3))
	out := make([]float64, n)
	for i := range out {
		out[i] = y[i] + dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i])
	}
	return out
}

// MOLIntegrate advances y0 through the given number of steps of size dt using
// the supplied stepper, starting at time t0, and returns the final state. It is
// the generic driver that turns a method-of-lines right-hand side into a full
// time integration. The input y0 is not modified.
func MOLIntegrate(f ODEFunc, y0 []float64, t0, dt float64, steps int, stepper Stepper) []float64 {
	y := Clone(y0)
	t := t0
	for s := 0; s < steps; s++ {
		y = stepper(f, t, y, dt)
		t += dt
	}
	return y
}

// MOLHeat1DRHS returns the semi-discrete right-hand side du/dt for the heat
// equation u_t = alpha*u_xx obtained by replacing u_xx with the three-point
// central difference on a grid of spacing dx. The endpoints are assigned a zero
// time derivative so that Dirichlet boundary values remain fixed under any of
// the time steppers. The returned slice is newly allocated.
func MOLHeat1DRHS(u []float64, alpha, dx float64) []float64 {
	n := len(u)
	out := make([]float64, n)
	inv := alpha / (dx * dx)
	for i := 1; i < n-1; i++ {
		out[i] = inv * (u[i-1] - 2*u[i] + u[i+1])
	}
	return out
}

// MOLHeat1DField returns an ODEFunc that evaluates MOLHeat1DRHS with the given
// diffusivity and grid spacing, suitable for passing to MOLIntegrate.
func MOLHeat1DField(alpha, dx float64) ODEFunc {
	return func(_ float64, u []float64) []float64 {
		return MOLHeat1DRHS(u, alpha, dx)
	}
}

// MOLAdvection1DRHS returns the semi-discrete right-hand side du/dt for the
// linear advection equation u_t + a*u_x = 0 using the second-order central
// difference for u_x on a periodic grid of spacing dx. The returned slice is
// newly allocated.
func MOLAdvection1DRHS(u []float64, a, dx float64) []float64 {
	n := len(u)
	out := make([]float64, n)
	c := a / (2 * dx)
	for i := 0; i < n; i++ {
		out[i] = -c * (u[numpdeWrap(i+1, n)] - u[numpdeWrap(i-1, n)])
	}
	return out
}

// MOLAdvection1DField returns an ODEFunc that evaluates MOLAdvection1DRHS with
// the given advection speed and grid spacing, suitable for MOLIntegrate.
func MOLAdvection1DField(a, dx float64) ODEFunc {
	return func(_ float64, u []float64) []float64 {
		return MOLAdvection1DRHS(u, a, dx)
	}
}
