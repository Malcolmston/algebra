package numpde

// Heat1DMeshRatio returns the dimensionless mesh (Fourier) ratio
// r = alpha*dt/dx^2 that governs the stability and accuracy of the
// finite-difference heat-equation schemes.
func Heat1DMeshRatio(alpha, dx, dt float64) float64 {
	return alpha * dt / (dx * dx)
}

// Heat1DStableDt returns the largest time step for which the explicit FTCS
// scheme applied to u_t = alpha*u_xx is stable, namely dt = dx^2/(2*alpha).
// Steps larger than this cause the explicit solver to blow up.
func Heat1DStableDt(alpha, dx float64) float64 {
	return dx * dx / (2 * alpha)
}

// Heat1DExplicitStep advances the field u by one time step of the explicit
// forward-time centred-space (FTCS) scheme with mesh ratio r = alpha*dt/dx^2,
// holding the two boundary values fixed (homogeneous or inhomogeneous
// Dirichlet). It returns a new slice and leaves u unmodified. The scheme is
// stable only for r <= 1/2.
func Heat1DExplicitStep(u []float64, r float64) []float64 {
	n := len(u)
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	out[0] = u[0]
	out[n-1] = u[n-1]
	for i := 1; i < n-1; i++ {
		out[i] = u[i] + r*(u[i-1]-2*u[i]+u[i+1])
	}
	return out
}

// Heat1DExplicit solves u_t = alpha*u_xx on the grid implied by u0 (spacing dx)
// by taking the given number of explicit FTCS steps of size dt. The boundary
// values u0[0] and u0[last] are held fixed for the whole simulation. The input
// is not modified; the field at the final time is returned. It panics if the
// scheme is unstable (mesh ratio above 1/2) so that silent garbage is never
// produced.
func Heat1DExplicit(u0 []float64, alpha, dx, dt float64, steps int) []float64 {
	r := Heat1DMeshRatio(alpha, dx, dt)
	if r > 0.5+1e-12 {
		panic("numpde: Heat1DExplicit unstable, mesh ratio alpha*dt/dx^2 > 1/2")
	}
	u := Clone(u0)
	for s := 0; s < steps; s++ {
		u = Heat1DExplicitStep(u, r)
	}
	return u
}

// Heat1DImplicitStep advances u by one step of the fully implicit backward-time
// centred-space (BTCS) scheme with mesh ratio r, holding the boundary values
// fixed. The scheme is unconditionally stable for any r > 0. A new slice is
// returned and u is left unchanged.
func Heat1DImplicitStep(u []float64, r float64) []float64 {
	n := len(u)
	if n < 3 {
		return Clone(u)
	}
	a := make([]float64, n)
	b := make([]float64, n)
	c := make([]float64, n)
	d := make([]float64, n)
	// Dirichlet boundary rows: identity.
	b[0], d[0] = 1, u[0]
	b[n-1], d[n-1] = 1, u[n-1]
	for i := 1; i < n-1; i++ {
		a[i] = -r
		b[i] = 1 + 2*r
		c[i] = -r
		d[i] = u[i]
	}
	return ThomasSolve(a, b, c, d)
}

// Heat1DImplicit solves u_t = alpha*u_xx with the unconditionally stable BTCS
// scheme, taking the given number of steps of size dt on a grid of spacing dx.
// Boundary values are held fixed. The input is not modified.
func Heat1DImplicit(u0 []float64, alpha, dx, dt float64, steps int) []float64 {
	r := Heat1DMeshRatio(alpha, dx, dt)
	u := Clone(u0)
	for s := 0; s < steps; s++ {
		u = Heat1DImplicitStep(u, r)
	}
	return u
}

// Heat1DCrankNicolsonStep advances u by one step of the second-order accurate,
// unconditionally stable Crank-Nicolson scheme with mesh ratio r, holding the
// boundary values fixed. A new slice is returned and u is left unchanged.
func Heat1DCrankNicolsonStep(u []float64, r float64) []float64 {
	n := len(u)
	if n < 3 {
		return Clone(u)
	}
	a := make([]float64, n)
	b := make([]float64, n)
	c := make([]float64, n)
	d := make([]float64, n)
	half := r / 2
	b[0], d[0] = 1, u[0]
	b[n-1], d[n-1] = 1, u[n-1]
	for i := 1; i < n-1; i++ {
		a[i] = -half
		b[i] = 1 + r
		c[i] = -half
		d[i] = u[i] + half*(u[i-1]-2*u[i]+u[i+1])
	}
	return ThomasSolve(a, b, c, d)
}

// Heat1DCrankNicolson solves u_t = alpha*u_xx with the Crank-Nicolson scheme,
// taking the given number of steps of size dt on a grid of spacing dx. Boundary
// values are held fixed. The scheme is second-order accurate in both space and
// time and unconditionally stable. The input is not modified.
func Heat1DCrankNicolson(u0 []float64, alpha, dx, dt float64, steps int) []float64 {
	r := Heat1DMeshRatio(alpha, dx, dt)
	u := Clone(u0)
	for s := 0; s < steps; s++ {
		u = Heat1DCrankNicolsonStep(u, r)
	}
	return u
}
