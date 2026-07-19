package spectralpde

import (
	"math"
	"math/cmplx"
)

// BuildChebyshevLaplacian1D returns the second-derivative (Laplacian in 1-D)
// Chebyshev collocation operator on [a, b], i.e. D^2 scaled by (2/(b-a))^2.
func BuildChebyshevLaplacian1D(N int, a, b float64) [][]float64 {
	s := IntervalScale(a, b)
	return MatScale(ChebyshevDiffMatrix2(N), s*s)
}

// restrictInterior extracts the interior submatrix (rows/cols 1..N-1) of a
// square (N+1)-by-(N+1) matrix.
func restrictInterior(a [][]float64) [][]float64 {
	N := len(a) - 1
	m := N - 1
	out := Zeros(m, m)
	for i := 1; i < N; i++ {
		for j := 1; j < N; j++ {
			out[i-1][j-1] = a[i][j]
		}
	}
	return out
}

// solveDirichlet1D solves the collocation system L*u = rhs with the two
// boundary nodes (indices 0 and N, i.e. x=b and x=a) fixed to u0 and uN.
func solveDirichlet1D(L [][]float64, rhs []float64, u0, uN float64) ([]float64, error) {
	N := len(L) - 1
	m := N - 1
	A := restrictInterior(L)
	b := make([]float64, m)
	for i := 1; i < N; i++ {
		b[i-1] = rhs[i] - L[i][0]*u0 - L[i][N]*uN
	}
	sol, err := SolveLinearSystem(A, b)
	if err != nil {
		return nil, err
	}
	u := make([]float64, N+1)
	u[0] = u0
	u[N] = uN
	for i := 1; i < N; i++ {
		u[i] = sol[i-1]
	}
	return u, nil
}

// PoissonSolve1D solves the Dirichlet problem u” = f on [a, b] with u(a) = ua
// and u(b) = ub, using an (N+1)-node Chebyshev collocation method. It returns
// the physical nodes and the nodal solution values (ordered from x=b down to
// x=a, matching the Chebyshev-Gauss-Lobatto ordering).
func PoissonSolve1D(f func(float64) float64, N int, a, b, ua, ub float64) (nodes, u []float64, err error) {
	nodes = ChebyshevGaussLobattoNodesInterval(N, a, b)
	L := BuildChebyshevLaplacian1D(N, a, b)
	rhs := ApplyFunc(f, nodes)
	// nodes[0] = b, nodes[N] = a, so boundary values are ub and ua.
	u, err = solveDirichlet1D(L, rhs, ub, ua)
	if err != nil {
		return nil, nil, err
	}
	return nodes, u, nil
}

// HelmholtzSolve1D solves the Dirichlet problem u” + k2*u = f on [a, b] with
// u(a) = ua and u(b) = ub, using Chebyshev collocation. k2 may be negative
// (modified Helmholtz).
func HelmholtzSolve1D(f func(float64) float64, k2 float64, N int, a, b, ua, ub float64) (nodes, u []float64, err error) {
	nodes = ChebyshevGaussLobattoNodesInterval(N, a, b)
	L := BuildChebyshevLaplacian1D(N, a, b)
	for i := 0; i <= N; i++ {
		L[i][i] += k2
	}
	rhs := ApplyFunc(f, nodes)
	u, err = solveDirichlet1D(L, rhs, ub, ua)
	if err != nil {
		return nil, nil, err
	}
	return nodes, u, nil
}

// PoissonSolve1DFourier solves the periodic problem u” = f on [0, L) using a
// Fourier spectral method. The source f must have zero mean for a solution to
// exist; the returned solution also has zero mean. It returns the grid nodes
// and nodal solution values.
func PoissonSolve1DFourier(f func(float64) float64, N int, L float64) (nodes, u []float64) {
	nodes = FourierNodesInterval(N, 0, L)
	fv := ApplyFunc(f, nodes)
	fhat := FFT(RealToComplex(fv))
	uhat := make([]complex128, N)
	scale := 2 * math.Pi / L
	for k := 0; k < N; k++ {
		kk := float64(k)
		if k > N/2 {
			kk = float64(k - N)
		}
		kp := scale * kk
		if kk == 0 {
			uhat[k] = 0
			continue
		}
		uhat[k] = fhat[k] / complex(-(kp*kp), 0)
	}
	u = ComplexReal(IFFT(uhat))
	return nodes, u
}

// HeatSolveFourier evolves the periodic heat equation u_t = nu*u_xx on [0, L)
// exactly in Fourier space from the initial data u0 (sampled on the Fourier
// grid) to time t, returning the nodal solution.
func HeatSolveFourier(u0 []float64, nu, t, L float64) []float64 {
	N := len(u0)
	fhat := FFT(RealToComplex(u0))
	scale := 2 * math.Pi / L
	for k := 0; k < N; k++ {
		kk := float64(k)
		if k > N/2 {
			kk = float64(k - N)
		}
		kp := scale * kk
		fhat[k] *= complex(math.Exp(-nu*kp*kp*t), 0)
	}
	return ComplexReal(IFFT(fhat))
}

// AdvectionDiffusionSolveFourier evolves the periodic equation
// u_t + c*u_x = nu*u_xx on [0, L) exactly in Fourier space from u0 to time t.
func AdvectionDiffusionSolveFourier(u0 []float64, c, nu, t, L float64) []float64 {
	N := len(u0)
	fhat := FFT(RealToComplex(u0))
	scale := 2 * math.Pi / L
	for k := 0; k < N; k++ {
		kk := float64(k)
		if k > N/2 {
			kk = float64(k - N)
		}
		kp := scale * kk
		lambda := complex(-nu*kp*kp, -c*kp)
		fhat[k] *= cmplx.Exp(lambda * complex(t, 0))
	}
	return ComplexReal(IFFT(fhat))
}

// AdvectionSolveFourier evolves the periodic linear advection equation
// u_t + c*u_x = 0 on [0, L) exactly to time t.
func AdvectionSolveFourier(u0 []float64, c, t, L float64) []float64 {
	return AdvectionDiffusionSolveFourier(u0, c, 0, t, L)
}

// HeatStepChebyshevCN performs one Crank-Nicolson time step of size dt for the
// heat equation u_t = nu*u_xx on [a, b] with homogeneous Dirichlet boundary
// conditions, given the current interior nodal values (length N-1). It returns
// the updated interior values.
func HeatStepChebyshevCN(uInterior []float64, nu, dt float64, N int, a, b float64) ([]float64, error) {
	L := BuildChebyshevLaplacian1D(N, a, b)
	Aint := restrictInterior(L)
	m := N - 1
	// (I - dt/2 nu L) u^{n+1} = (I + dt/2 nu L) u^n.
	lhs := Zeros(m, m)
	rhsM := Zeros(m, m)
	for i := 0; i < m; i++ {
		for j := 0; j < m; j++ {
			lhs[i][j] = -0.5 * dt * nu * Aint[i][j]
			rhsM[i][j] = 0.5 * dt * nu * Aint[i][j]
			if i == j {
				lhs[i][j] += 1
				rhsM[i][j] += 1
			}
		}
	}
	rhs := MatVec(rhsM, uInterior)
	return SolveLinearSystem(lhs, rhs)
}

// HeatSolveChebyshev integrates u_t = nu*u_xx on [a, b] with homogeneous
// Dirichlet boundary conditions from the initial function u0 to time
// steps*dt, using Crank-Nicolson in time and Chebyshev collocation in space.
// It returns the physical nodes and the full nodal solution (boundaries
// included, set to zero).
func HeatSolveChebyshev(u0 func(float64) float64, nu, dt float64, steps, N int, a, b float64) (nodes, u []float64, err error) {
	nodes = ChebyshevGaussLobattoNodesInterval(N, a, b)
	full := ApplyFunc(u0, nodes)
	interior := make([]float64, N-1)
	for i := 1; i < N; i++ {
		interior[i-1] = full[i]
	}
	for s := 0; s < steps; s++ {
		interior, err = HeatStepChebyshevCN(interior, nu, dt, N, a, b)
		if err != nil {
			return nil, nil, err
		}
	}
	u = make([]float64, N+1)
	for i := 1; i < N; i++ {
		u[i] = interior[i-1]
	}
	return nodes, u, nil
}

// SpectralDerivativeChebyshev returns the nodal derivative values of a function
// sampled at the Chebyshev-Gauss-Lobatto nodes on [a, b].
func SpectralDerivativeChebyshev(values []float64, a, b float64) []float64 {
	N := len(values) - 1
	return MatVec(ChebyshevDiffMatrixInterval(N, a, b), values)
}
