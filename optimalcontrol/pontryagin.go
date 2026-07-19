package optimalcontrol

// TPBVPSolution holds the sampled solution of a two-point boundary-value problem
// arising from Pontryagin's minimum principle: the state, costate and optimal
// control trajectories on a time grid.
type TPBVPSolution struct {
	// Times holds the ascending sample times from 0 to T.
	Times []float64
	// X holds the state at each sample time.
	X [][]float64
	// P holds the costate (adjoint) at each sample time.
	P [][]float64
	// U holds the optimal control at each sample time.
	U [][]float64
}

// OptimalControlLQ returns the stationarity-condition control u = −R⁻¹ Bᵀ p that
// minimizes the linear-quadratic Hamiltonian for a given costate p.
func OptimalControlLQ(b, r *Matrix, p []float64) ([]float64, error) {
	rinv, err := Inverse(r)
	if err != nil {
		return nil, err
	}
	u := rinv.Mul(b.Transpose()).MulVec(p)
	for i := range u {
		u[i] = -u[i]
	}
	return u, nil
}

// CostateDynamicsLQ returns the costate derivative p' = −Q x − Aᵀ p for the
// linear-quadratic problem.
func CostateDynamicsLQ(a, q *Matrix, x, p []float64) []float64 {
	qx := q.MulVec(x)
	atp := a.Transpose().MulVec(p)
	out := make([]float64, len(x))
	for i := range out {
		out[i] = -qx[i] - atp[i]
	}
	return out
}

// PontryaginLQ solves the fixed-final-time linear-quadratic optimal-control
// problem
//
//	min  ½ x(T)ᵀ Qf x(T) + ½ ∫₀ᵀ (xᵀ Q x + uᵀ R u) dt
//	s.t. x' = A x + B u,  x(0) = x0,
//
// via Pontryagin's minimum principle. It forms the Hamiltonian two-point
// boundary-value problem, propagates it with the Hamiltonian state-transition
// matrix, solves the terminal transversality condition p(T) = Qf x(T) for the
// initial costate, and returns the sampled state, costate and control.
func PontryaginLQ(a, b, q, r, qf *Matrix, x0 []float64, tFinal float64, steps int) (*TPBVPSolution, error) {
	if steps < 1 {
		return nil, ErrDim
	}
	n := a.rows
	m, err := HamiltonianMatrix(a, b, q, r)
	if err != nil {
		return nil, err
	}
	phi := MatrixExp(m.Scale(tFinal))
	phi11 := phi.Submatrix(0, n, 0, n)
	phi12 := phi.Submatrix(0, n, n, 2*n)
	phi21 := phi.Submatrix(n, 2*n, 0, n)
	phi22 := phi.Submatrix(n, 2*n, n, 2*n)
	// (phi22 - Qf phi12) p0 = (Qf phi11 - phi21) x0.
	lhs := phi22.Minus(qf.Mul(phi12))
	rhsMat := qf.Mul(phi11).Minus(phi21)
	rhs := rhsMat.MulVec(x0)
	p0, err := Solve(lhs, rhs)
	if err != nil {
		return nil, err
	}
	// Initial augmented state [x0; p0].
	z0 := append(append([]float64{}, x0...), p0...)
	sol := &TPBVPSolution{
		Times: make([]float64, steps+1),
		X:     make([][]float64, steps+1),
		P:     make([][]float64, steps+1),
		U:     make([][]float64, steps+1),
	}
	dt := tFinal / float64(steps)
	for k := 0; k <= steps; k++ {
		t := float64(k) * dt
		phiT := MatrixExp(m.Scale(t))
		z := phiT.MulVec(z0)
		x := z[:n]
		p := z[n:]
		u, err := OptimalControlLQ(b, r, p)
		if err != nil {
			return nil, err
		}
		sol.Times[k] = t
		sol.X[k] = append([]float64{}, x...)
		sol.P[k] = append([]float64{}, p...)
		sol.U[k] = u
	}
	return sol, nil
}

// IndirectShooting solves a fixed-final-time two-point boundary-value problem by
// single shooting on the initial costate. The problem is specified by the
// augmented dynamics
//
//	x' = f(x, u),   p' = g(x, p, u),   u = uOpt(x, p),
//
// with x(0) = x0 and terminal condition p(T) = termGrad(x(T)). Newton's method
// with a finite-difference Jacobian is applied to the shooting residual. It
// returns the initial costate and the sampled trajectory.
func IndirectShooting(
	x0 []float64,
	f func(x, u []float64) []float64,
	g func(x, p, u []float64) []float64,
	uOpt func(x, p []float64) []float64,
	termGrad func(x []float64) []float64,
	tFinal float64, steps, maxIter int, tol float64,
) ([]float64, *TPBVPSolution, error) {
	n := len(x0)
	// residual(p0) = p(T) - termGrad(x(T)).
	integrate := func(p0 []float64) (xT, pT []float64, sol *TPBVPSolution) {
		x := append([]float64{}, x0...)
		p := append([]float64{}, p0...)
		dt := tFinal / float64(steps)
		sol = &TPBVPSolution{
			Times: make([]float64, steps+1),
			X:     make([][]float64, steps+1),
			P:     make([][]float64, steps+1),
			U:     make([][]float64, steps+1),
		}
		record := func(k int, x, p []float64) {
			u := uOpt(x, p)
			sol.Times[k] = float64(k) * dt
			sol.X[k] = append([]float64{}, x...)
			sol.P[k] = append([]float64{}, p...)
			sol.U[k] = u
		}
		record(0, x, p)
		for k := 0; k < steps; k++ {
			x, p = rk4Augmented(x, p, f, g, uOpt, dt)
			record(k+1, x, p)
		}
		return x, p, sol
	}
	residual := func(p0 []float64) []float64 {
		xT, pT, _ := integrate(p0)
		grad := termGrad(xT)
		res := make([]float64, n)
		for i := 0; i < n; i++ {
			res[i] = pT[i] - grad[i]
		}
		return res
	}
	p0 := make([]float64, n)
	var sol *TPBVPSolution
	for iter := 0; iter < maxIter; iter++ {
		r := residual(p0)
		var norm float64
		for _, v := range r {
			norm += v * v
		}
		if norm < tol*tol {
			_, _, sol = integrate(p0)
			return p0, sol, nil
		}
		// Finite-difference Jacobian J[i][j] = d r_i / d p0_j.
		jac := Zeros(n, n)
		h := 1e-6
		for j := 0; j < n; j++ {
			pp := append([]float64{}, p0...)
			pp[j] += h
			rp := residual(pp)
			for i := 0; i < n; i++ {
				jac.Set(i, j, (rp[i]-r[i])/h)
			}
		}
		delta, err := Solve(jac, r)
		if err != nil {
			return nil, nil, err
		}
		for j := 0; j < n; j++ {
			p0[j] -= delta[j]
		}
	}
	_, _, sol = integrate(p0)
	return p0, sol, ErrNotConverged
}

// rk4Augmented advances the augmented state/costate system by one RK4 step.
func rk4Augmented(
	x, p []float64,
	f func(x, u []float64) []float64,
	g func(x, p, u []float64) []float64,
	uOpt func(x, p []float64) []float64,
	dt float64,
) (nx, np []float64) {
	n := len(x)
	add := func(a, b []float64, s float64) []float64 {
		out := make([]float64, len(a))
		for i := range a {
			out[i] = a[i] + s*b[i]
		}
		return out
	}
	deriv := func(x, p []float64) (dx, dp []float64) {
		u := uOpt(x, p)
		return f(x, u), g(x, p, u)
	}
	k1x, k1p := deriv(x, p)
	k2x, k2p := deriv(add(x, k1x, dt/2), add(p, k1p, dt/2))
	k3x, k3p := deriv(add(x, k2x, dt/2), add(p, k2p, dt/2))
	k4x, k4p := deriv(add(x, k3x, dt), add(p, k3p, dt))
	nx = make([]float64, n)
	np = make([]float64, n)
	for i := 0; i < n; i++ {
		nx[i] = x[i] + dt/6*(k1x[i]+2*k2x[i]+2*k3x[i]+k4x[i])
		np[i] = p[i] + dt/6*(k1p[i]+2*k2p[i]+2*k3p[i]+k4p[i])
	}
	return nx, np
}
