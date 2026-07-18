package numpde

// Wave1DCourant returns the Courant number C = c*dt/dx for the 1D wave
// equation. The explicit leap-frog scheme is stable for C <= 1.
func Wave1DCourant(c, dx, dt float64) float64 { return c * dt / dx }

// Wave1DStableDt returns the largest stable explicit time step for the 1D wave
// equation, dt = dx/c (the Courant-Friedrichs-Lewy limit).
func Wave1DStableDt(c, dx float64) float64 { return dx / c }

// Wave1DStep advances one leap-frog step of the 1D wave equation given the
// solution at the previous two time levels. prev is u at time level n-1, curr
// is u at level n, and c2 is the squared Courant number (c*dt/dx)^2. Boundary
// values are copied from curr (fixed Dirichlet). A new slice with u at level
// n+1 is returned.
func Wave1DStep(prev, curr []float64, c2 float64) []float64 {
	n := len(curr)
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	out[0] = curr[0]
	out[n-1] = curr[n-1]
	for i := 1; i < n-1; i++ {
		out[i] = 2*curr[i] - prev[i] + c2*(curr[i-1]-2*curr[i]+curr[i+1])
	}
	return out
}

// Wave1DExplicit solves u_tt = c^2*u_xx with the explicit leap-frog (central
// time, central space) scheme. u0 is the initial displacement and v0 the
// initial velocity, both sampled on a grid of spacing dx. The solver performs
// the given number of time steps of size dt using the standard second-order
// accurate startup step and holds the boundary values of u0 fixed. It panics if
// the Courant number exceeds 1 (unstable). The inputs are not modified.
func Wave1DExplicit(u0, v0 []float64, c, dx, dt float64, steps int) []float64 {
	if len(u0) != len(v0) {
		panic("numpde: Wave1DExplicit u0 and v0 length mismatch")
	}
	C := Wave1DCourant(c, dx, dt)
	if C > 1+1e-12 {
		panic("numpde: Wave1DExplicit unstable, Courant number > 1")
	}
	n := len(u0)
	c2 := C * C
	prev := Clone(u0)
	// Second-order accurate first step using the initial velocity.
	curr := make([]float64, n)
	curr[0] = u0[0]
	curr[n-1] = u0[n-1]
	for i := 1; i < n-1; i++ {
		curr[i] = u0[i] + dt*v0[i] + 0.5*c2*(u0[i-1]-2*u0[i]+u0[i+1])
	}
	if steps <= 0 {
		return prev
	}
	if steps == 1 {
		return curr
	}
	for s := 1; s < steps; s++ {
		next := Wave1DStep(prev, curr, c2)
		prev, curr = curr, next
	}
	return curr
}

// Wave2DCourant returns the 2D Courant number c*dt*sqrt(1/dx^2 + 1/dy^2). The
// explicit scheme is stable when this quantity does not exceed 1.
func Wave2DCourant(c, dx, dy, dt float64) float64 {
	return c * dt * numpdeSqrt(1/(dx*dx)+1/(dy*dy))
}

// Wave2DExplicit solves u_tt = c^2*(u_xx + u_yy) with the explicit leap-frog
// scheme. u0 is the initial displacement and v0 the initial velocity, both
// stored as [i][j] on a grid with spacings dx and dy. The solver performs the
// given number of steps of size dt and holds all boundary nodes fixed. It
// panics if the 2D Courant number exceeds 1. The inputs are not modified.
func Wave2DExplicit(u0, v0 [][]float64, c, dx, dy, dt float64, steps int) [][]float64 {
	if Wave2DCourant(c, dx, dy, dt) > 1+1e-12 {
		panic("numpde: Wave2DExplicit unstable, Courant number > 1")
	}
	nx := len(u0)
	ny := len(u0[0])
	cx := c * c * dt * dt / (dx * dx)
	cy := c * c * dt * dt / (dy * dy)
	prev := Clone2D(u0)
	curr := Clone2D(u0)
	for i := 1; i < nx-1; i++ {
		for j := 1; j < ny-1; j++ {
			curr[i][j] = u0[i][j] + dt*v0[i][j] +
				0.5*cx*(u0[i-1][j]-2*u0[i][j]+u0[i+1][j]) +
				0.5*cy*(u0[i][j-1]-2*u0[i][j]+u0[i][j+1])
		}
	}
	if steps <= 0 {
		return prev
	}
	if steps == 1 {
		return curr
	}
	for s := 1; s < steps; s++ {
		next := Clone2D(curr)
		for i := 1; i < nx-1; i++ {
			for j := 1; j < ny-1; j++ {
				next[i][j] = 2*curr[i][j] - prev[i][j] +
					cx*(curr[i-1][j]-2*curr[i][j]+curr[i+1][j]) +
					cy*(curr[i][j-1]-2*curr[i][j]+curr[i][j+1])
			}
		}
		prev, curr = curr, next
	}
	return curr
}
