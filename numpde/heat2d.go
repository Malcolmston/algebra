package numpde

// Heat2DMeshRatio returns the combined explicit stability parameter
// alpha*dt*(1/dx^2 + 1/dy^2) for the 2D heat equation. The explicit scheme is
// stable when this quantity does not exceed 1/2.
func Heat2DMeshRatio(alpha, dx, dy, dt float64) float64 {
	return alpha * dt * (1/(dx*dx) + 1/(dy*dy))
}

// Heat2DStableDt returns the largest stable explicit time step for
// u_t = alpha*(u_xx + u_yy), namely dt = 1 / (2*alpha*(1/dx^2 + 1/dy^2)).
func Heat2DStableDt(alpha, dx, dy float64) float64 {
	return 1.0 / (2 * alpha * (1/(dx*dx) + 1/(dy*dy)))
}

// Heat2DExplicitStep advances the field u (stored as u[i][j]) by one explicit
// FTCS step for u_t = alpha*(u_xx + u_yy), holding every boundary node fixed.
// rx = alpha*dt/dx^2 and ry = alpha*dt/dy^2 are the directional mesh ratios. A
// new matrix is returned; u is left unmodified.
func Heat2DExplicitStep(u [][]float64, rx, ry float64) [][]float64 {
	nx := len(u)
	ny := len(u[0])
	out := Clone2D(u)
	for i := 1; i < nx-1; i++ {
		for j := 1; j < ny-1; j++ {
			out[i][j] = u[i][j] +
				rx*(u[i-1][j]-2*u[i][j]+u[i+1][j]) +
				ry*(u[i][j-1]-2*u[i][j]+u[i][j+1])
		}
	}
	return out
}

// Heat2DExplicit solves u_t = alpha*(u_xx + u_yy) on the grid implied by u0
// (spacings dx, dy) with the explicit FTCS scheme, taking the given number of
// steps of size dt. All boundary nodes are held fixed at their initial values.
// It panics if the combined mesh ratio exceeds 1/2 (unstable). The input is not
// modified.
func Heat2DExplicit(u0 [][]float64, alpha, dx, dy, dt float64, steps int) [][]float64 {
	if Heat2DMeshRatio(alpha, dx, dy, dt) > 0.5+1e-12 {
		panic("numpde: Heat2DExplicit unstable, combined mesh ratio > 1/2")
	}
	rx := alpha * dt / (dx * dx)
	ry := alpha * dt / (dy * dy)
	u := Clone2D(u0)
	for s := 0; s < steps; s++ {
		u = Heat2DExplicitStep(u, rx, ry)
	}
	return u
}

// Heat2DADIStep advances u by one Peaceman-Rachford alternating-direction
// implicit (ADI) step for u_t = alpha*(u_xx + u_yy). Each step consists of two
// half-steps: implicit in x then implicit in y, each solved with the Thomas
// algorithm. The scheme is second-order accurate in time and unconditionally
// stable. Boundary nodes are held fixed. A new matrix is returned.
func Heat2DADIStep(u [][]float64, alpha, dx, dy, dt float64) [][]float64 {
	nx := len(u)
	ny := len(u[0])
	rx := alpha * dt / (2 * dx * dx)
	ry := alpha * dt / (2 * dy * dy)

	// Half-step 1: implicit in x for each interior row j.
	star := Clone2D(u)
	for j := 1; j < ny-1; j++ {
		a := make([]float64, nx)
		b := make([]float64, nx)
		c := make([]float64, nx)
		d := make([]float64, nx)
		b[0], d[0] = 1, u[0][j]
		b[nx-1], d[nx-1] = 1, u[nx-1][j]
		for i := 1; i < nx-1; i++ {
			a[i] = -rx
			b[i] = 1 + 2*rx
			c[i] = -rx
			d[i] = u[i][j] + ry*(u[i][j-1]-2*u[i][j]+u[i][j+1])
		}
		col := ThomasSolve(a, b, c, d)
		for i := 0; i < nx; i++ {
			star[i][j] = col[i]
		}
	}

	// Half-step 2: implicit in y for each interior column i.
	next := Clone2D(star)
	for i := 1; i < nx-1; i++ {
		a := make([]float64, ny)
		b := make([]float64, ny)
		c := make([]float64, ny)
		d := make([]float64, ny)
		b[0], d[0] = 1, u[i][0]
		b[ny-1], d[ny-1] = 1, u[i][ny-1]
		for j := 1; j < ny-1; j++ {
			a[j] = -ry
			b[j] = 1 + 2*ry
			c[j] = -ry
			d[j] = star[i][j] + rx*(star[i-1][j]-2*star[i][j]+star[i+1][j])
		}
		row := ThomasSolve(a, b, c, d)
		for j := 0; j < ny; j++ {
			next[i][j] = row[j]
		}
	}
	return next
}

// Heat2DADI solves u_t = alpha*(u_xx + u_yy) with the unconditionally stable
// Peaceman-Rachford ADI scheme, taking the given number of steps of size dt on
// a grid with spacings dx and dy. All boundary nodes are held fixed at their
// initial values. The input is not modified.
func Heat2DADI(u0 [][]float64, alpha, dx, dy, dt float64, steps int) [][]float64 {
	u := Clone2D(u0)
	for s := 0; s < steps; s++ {
		u = Heat2DADIStep(u, alpha, dx, dy, dt)
	}
	return u
}
