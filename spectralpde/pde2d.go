package spectralpde

import "math"

// BuildLaplacian2D returns the 2-D Chebyshev collocation Laplacian on the
// tensor-product grid of [ax, bx] x [ay, by] with Nx+1 nodes in x and Ny+1
// nodes in y. The unknowns are ordered lexicographically with the x-index
// varying slowest: index(i, j) = i*(Ny+1) + j. The operator is
// D2x kron Iy + Ix kron D2y.
func BuildLaplacian2D(Nx, Ny int, ax, bx, ay, by float64) [][]float64 {
	d2x := BuildChebyshevLaplacian1D(Nx, ax, bx)
	d2y := BuildChebyshevLaplacian1D(Ny, ay, by)
	ix := Identity(Nx + 1)
	iy := Identity(Ny + 1)
	return MatAdd(Kron(d2x, iy), Kron(ix, d2y))
}

// index2D returns the flattened index for grid point (i, j).
func index2D(i, j, Ny int) int { return i*(Ny+1) + j }

// isBoundary2D reports whether grid point (i, j) lies on the boundary.
func isBoundary2D(i, j, Nx, Ny int) bool {
	return i == 0 || i == Nx || j == 0 || j == Ny
}

// Poisson2D solves the Dirichlet problem u_xx + u_yy = f on the rectangle
// [ax, bx] x [ay, by] with boundary values given by g(x, y). It uses Chebyshev
// collocation with Nx+1 by Ny+1 nodes and returns the x-nodes, y-nodes and the
// solution as a matrix U with U[i][j] = u(x_i, y_j).
func Poisson2D(f, g func(x, y float64) float64, Nx, Ny int, ax, bx, ay, by float64) (xn, yn []float64, U [][]float64, err error) {
	return helmholtz2D(f, g, 0, Nx, Ny, ax, bx, ay, by)
}

// Helmholtz2D solves u_xx + u_yy + k2*u = f on [ax, bx] x [ay, by] with
// Dirichlet data g, by Chebyshev collocation.
func Helmholtz2D(f, g func(x, y float64) float64, k2 float64, Nx, Ny int, ax, bx, ay, by float64) (xn, yn []float64, U [][]float64, err error) {
	return helmholtz2D(f, g, k2, Nx, Ny, ax, bx, ay, by)
}

func helmholtz2D(f, g func(x, y float64) float64, k2 float64, Nx, Ny int, ax, bx, ay, by float64) (xn, yn []float64, U [][]float64, err error) {
	xn = ChebyshevGaussLobattoNodesInterval(Nx, ax, bx)
	yn = ChebyshevGaussLobattoNodesInterval(Ny, ay, by)
	L := BuildLaplacian2D(Nx, Ny, ax, bx, ay, by)
	M := (Nx + 1) * (Ny + 1)
	if k2 != 0 {
		for i := 0; i < M; i++ {
			L[i][i] += k2
		}
	}
	// Full boundary/known values and source.
	uFull := make([]float64, M)
	fFull := make([]float64, M)
	interiorIdx := make([]int, 0, M)
	for i := 0; i <= Nx; i++ {
		for j := 0; j <= Ny; j++ {
			p := index2D(i, j, Ny)
			fFull[p] = f(xn[i], yn[j])
			if isBoundary2D(i, j, Nx, Ny) {
				uFull[p] = g(xn[i], yn[j])
			} else {
				interiorIdx = append(interiorIdx, p)
			}
		}
	}
	m := len(interiorIdx)
	A := Zeros(m, m)
	rhs := make([]float64, m)
	for a := 0; a < m; a++ {
		p := interiorIdx[a]
		s := fFull[p]
		// Subtract boundary contributions.
		for q := 0; q < M; q++ {
			if L[p][q] != 0 && isBoundaryIndex(q, Nx, Ny) {
				s -= L[p][q] * uFull[q]
			}
		}
		rhs[a] = s
		for b := 0; b < m; b++ {
			A[a][b] = L[p][interiorIdx[b]]
		}
	}
	sol, err := SolveLinearSystem(A, rhs)
	if err != nil {
		return nil, nil, nil, err
	}
	for a := 0; a < m; a++ {
		uFull[interiorIdx[a]] = sol[a]
	}
	U = Zeros(Nx+1, Ny+1)
	for i := 0; i <= Nx; i++ {
		for j := 0; j <= Ny; j++ {
			U[i][j] = uFull[index2D(i, j, Ny)]
		}
	}
	return xn, yn, U, nil
}

// isBoundaryIndex reports whether the flattened index p is a boundary node of
// the (Nx+1) x (Ny+1) grid.
func isBoundaryIndex(p, Nx, Ny int) bool {
	i := p / (Ny + 1)
	j := p % (Ny + 1)
	return isBoundary2D(i, j, Nx, Ny)
}

// Poisson2DFourier solves the doubly periodic problem u_xx + u_yy = f on
// [0, Lx) x [0, Ly) with a Fourier spectral method. The source f must have
// zero mean; the returned solution also has zero mean. It returns the x-nodes,
// y-nodes and the solution matrix U[i][j] = u(x_i, y_j).
func Poisson2DFourier(f func(x, y float64) float64, Nx, Ny int, Lx, Ly float64) (xn, yn []float64, U [][]float64) {
	xn = FourierNodesInterval(Nx, 0, Lx)
	yn = FourierNodesInterval(Ny, 0, Ly)
	// 2-D DFT via row/column 1-D transforms.
	grid := make([][]complex128, Nx)
	for i := 0; i < Nx; i++ {
		row := make([]complex128, Ny)
		for j := 0; j < Ny; j++ {
			row[j] = complex(f(xn[i], yn[j]), 0)
		}
		grid[i] = FFT(row)
	}
	// Transform along columns.
	hat := make([][]complex128, Nx)
	for i := range hat {
		hat[i] = make([]complex128, Ny)
	}
	col := make([]complex128, Nx)
	for j := 0; j < Ny; j++ {
		for i := 0; i < Nx; i++ {
			col[i] = grid[i][j]
		}
		colHat := FFT(col)
		for i := 0; i < Nx; i++ {
			hat[i][j] = colHat[i]
		}
	}
	sx := 2 * math.Pi / Lx
	sy := 2 * math.Pi / Ly
	for i := 0; i < Nx; i++ {
		kx := float64(i)
		if i > Nx/2 {
			kx = float64(i - Nx)
		}
		kx *= sx
		for j := 0; j < Ny; j++ {
			ky := float64(j)
			if j > Ny/2 {
				ky = float64(j - Ny)
			}
			ky *= sy
			denom := -(kx*kx + ky*ky)
			if denom == 0 {
				hat[i][j] = 0
			} else {
				hat[i][j] /= complex(denom, 0)
			}
		}
	}
	// Inverse 2-D transform.
	for j := 0; j < Ny; j++ {
		for i := 0; i < Nx; i++ {
			col[i] = hat[i][j]
		}
		colInv := IFFT(col)
		for i := 0; i < Nx; i++ {
			hat[i][j] = colInv[i]
		}
	}
	U = Zeros(Nx, Ny)
	for i := 0; i < Nx; i++ {
		rowInv := IFFT(hat[i])
		for j := 0; j < Ny; j++ {
			U[i][j] = real(rowInv[j])
		}
	}
	return xn, yn, U
}
