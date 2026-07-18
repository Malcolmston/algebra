package numpde

// The advection solvers in this file integrate the constant-coefficient linear
// advection (transport) equation
//
//	u_t + a*u_x = 0
//
// on a uniform periodic grid. The field slice is treated as one full period of
// n independent cells with wrap-around indexing, so it should be sampled at
// x_i = x0 + i*dx for i in [0, n-1] with dx equal to the period length divided
// by n (the right endpoint is NOT duplicated). The exact solution simply
// translates the initial profile: u(x, t) = u0(x - a*t) taken periodically.

// numpdeWrap returns i reduced into the range [0, n) with periodic wrap-around.
func numpdeWrap(i, n int) int {
	i %= n
	if i < 0 {
		i += n
	}
	return i
}

// Advection1DCourant returns the (signed) Courant number nu = a*dt/dx. The
// explicit schemes are stable for |nu| <= 1.
func Advection1DCourant(a, dx, dt float64) float64 { return a * dt / dx }

// Advection1DStableDt returns the largest stable time step dx/|a| for the
// explicit advection schemes (the CFL limit). It panics if a is zero.
func Advection1DStableDt(a, dx float64) float64 {
	if a == 0 {
		panic("numpde: Advection1DStableDt undefined for zero advection speed")
	}
	return dx / numpdeAbs(a)
}

// Advection1DUpwindStep advances one step of the first-order upwind scheme on a
// periodic grid with Courant number nu = a*dt/dx. The upwind direction is
// chosen from the sign of nu. A new slice is returned; u is left unmodified.
func Advection1DUpwindStep(u []float64, nu float64) []float64 {
	n := len(u)
	out := make([]float64, n)
	if nu >= 0 {
		for i := 0; i < n; i++ {
			out[i] = u[i] - nu*(u[i]-u[numpdeWrap(i-1, n)])
		}
	} else {
		for i := 0; i < n; i++ {
			out[i] = u[i] - nu*(u[numpdeWrap(i+1, n)]-u[i])
		}
	}
	return out
}

// Advection1DUpwind solves u_t + a*u_x = 0 on a periodic grid of spacing dx with
// the first-order upwind scheme, taking the given number of steps of size dt.
// It panics if the Courant number exceeds 1 in magnitude (unstable). The input
// is not modified.
func Advection1DUpwind(u0 []float64, a, dx, dt float64, steps int) []float64 {
	nu := Advection1DCourant(a, dx, dt)
	if numpdeAbs(nu) > 1+1e-12 {
		panic("numpde: Advection1DUpwind unstable, |Courant| > 1")
	}
	u := Clone(u0)
	for s := 0; s < steps; s++ {
		u = Advection1DUpwindStep(u, nu)
	}
	return u
}

// Advection1DLaxFriedrichsStep advances one step of the Lax-Friedrichs scheme on
// a periodic grid with Courant number nu. The scheme is first-order accurate and
// quite diffusive but robust. A new slice is returned.
func Advection1DLaxFriedrichsStep(u []float64, nu float64) []float64 {
	n := len(u)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		l := u[numpdeWrap(i-1, n)]
		r := u[numpdeWrap(i+1, n)]
		out[i] = 0.5*(l+r) - 0.5*nu*(r-l)
	}
	return out
}

// Advection1DLaxFriedrichs solves u_t + a*u_x = 0 on a periodic grid with the
// Lax-Friedrichs scheme, taking the given number of steps. It panics if the
// Courant number exceeds 1 in magnitude. The input is not modified.
func Advection1DLaxFriedrichs(u0 []float64, a, dx, dt float64, steps int) []float64 {
	nu := Advection1DCourant(a, dx, dt)
	if numpdeAbs(nu) > 1+1e-12 {
		panic("numpde: Advection1DLaxFriedrichs unstable, |Courant| > 1")
	}
	u := Clone(u0)
	for s := 0; s < steps; s++ {
		u = Advection1DLaxFriedrichsStep(u, nu)
	}
	return u
}

// Advection1DLaxWendroffStep advances one step of the second-order accurate
// Lax-Wendroff scheme on a periodic grid with Courant number nu. A new slice is
// returned.
func Advection1DLaxWendroffStep(u []float64, nu float64) []float64 {
	n := len(u)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		l := u[numpdeWrap(i-1, n)]
		r := u[numpdeWrap(i+1, n)]
		out[i] = u[i] - 0.5*nu*(r-l) + 0.5*nu*nu*(l-2*u[i]+r)
	}
	return out
}

// Advection1DLaxWendroff solves u_t + a*u_x = 0 on a periodic grid with the
// second-order Lax-Wendroff scheme, taking the given number of steps. It panics
// if the Courant number exceeds 1 in magnitude. The input is not modified.
func Advection1DLaxWendroff(u0 []float64, a, dx, dt float64, steps int) []float64 {
	nu := Advection1DCourant(a, dx, dt)
	if numpdeAbs(nu) > 1+1e-12 {
		panic("numpde: Advection1DLaxWendroff unstable, |Courant| > 1")
	}
	u := Clone(u0)
	for s := 0; s < steps; s++ {
		u = Advection1DLaxWendroffStep(u, nu)
	}
	return u
}
