package spectralpde

import "math"

// ChebyshevDiffMatrix returns the (N+1)-by-(N+1) Chebyshev collocation
// differentiation matrix D on the Chebyshev-Gauss-Lobatto nodes. If u holds the
// nodal values of a function then D*u holds the nodal values of its derivative
// (Trefethen, Spectral Methods in MATLAB, program cheb).
func ChebyshevDiffMatrix(N int) [][]float64 {
	if N == 0 {
		return [][]float64{{0}}
	}
	x := ChebyshevGaussLobattoNodes(N)
	c := make([]float64, N+1)
	for j := 0; j <= N; j++ {
		c[j] = 1
		if j == 0 || j == N {
			c[j] = 2
		}
		if j%2 == 1 {
			c[j] = -c[j]
		}
	}
	D := Zeros(N+1, N+1)
	for i := 0; i <= N; i++ {
		for j := 0; j <= N; j++ {
			if i != j {
				D[i][j] = (c[i] / c[j]) / (x[i] - x[j])
			}
		}
	}
	// Diagonal entries from the negative row sums.
	for i := 0; i <= N; i++ {
		var s float64
		for j := 0; j <= N; j++ {
			if j != i {
				s += D[i][j]
			}
		}
		D[i][i] = -s
	}
	return D
}

// ChebyshevDiffMatrixOrder returns the m-th order Chebyshev differentiation
// matrix, computed as the m-th matrix power of the first-order matrix.
func ChebyshevDiffMatrixOrder(N, m int) [][]float64 {
	D := ChebyshevDiffMatrix(N)
	return MatPow(D, m)
}

// ChebyshevDiffMatrix2 returns the second-order Chebyshev differentiation
// matrix D^2.
func ChebyshevDiffMatrix2(N int) [][]float64 {
	return ChebyshevDiffMatrixOrder(N, 2)
}

// ChebyshevDiffMatrixInterval returns the first-order Chebyshev
// differentiation matrix scaled for the physical interval [a, b].
func ChebyshevDiffMatrixInterval(N int, a, b float64) [][]float64 {
	return MatScale(ChebyshevDiffMatrix(N), IntervalScale(a, b))
}

// FourierDiffMatrix returns the N-by-N first-order Fourier differentiation
// matrix on the periodic grid x_j = 2*pi*j/N. It is implemented for even N,
// following Trefethen's program 4.
func FourierDiffMatrix(N int) [][]float64 {
	D := Zeros(N, N)
	h := 2 * math.Pi / float64(N)
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			if i == j {
				D[i][j] = 0
				continue
			}
			diff := float64(i-j) * h / 2
			sign := 1.0
			if (i-j)%2 != 0 {
				sign = -1.0
			}
			if N%2 == 0 {
				D[i][j] = 0.5 * sign / math.Tan(diff)
			} else {
				D[i][j] = 0.5 * sign / math.Sin(diff)
			}
		}
	}
	return D
}

// FourierDiffMatrix2 returns the N-by-N second-order Fourier differentiation
// matrix on the periodic grid, for even N.
func FourierDiffMatrix2(N int) [][]float64 {
	D := Zeros(N, N)
	h := 2 * math.Pi / float64(N)
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			if i == j {
				if N%2 == 0 {
					D[i][j] = -math.Pi*math.Pi/(3*h*h) - 1.0/6.0
				} else {
					D[i][j] = -math.Pi*math.Pi/(3*h*h) + 1.0/12.0
				}
				continue
			}
			diff := float64(i-j) * h / 2
			sign := 1.0
			if (i-j)%2 != 0 {
				sign = -1.0
			}
			if N%2 == 0 {
				s := math.Sin(diff)
				D[i][j] = -sign / (2 * s * s)
			} else {
				s := math.Sin(diff)
				co := math.Cos(diff)
				D[i][j] = -sign * co / (2 * s * s)
			}
		}
	}
	return D
}

// FourierDiffMatrixOrder returns the m-th order Fourier differentiation matrix.
// For m == 1 and m == 2 the closed-form matrices are used; for higher orders
// the matrix power of the first-order matrix is returned.
func FourierDiffMatrixOrder(N, m int) [][]float64 {
	switch m {
	case 0:
		return Identity(N)
	case 1:
		return FourierDiffMatrix(N)
	case 2:
		return FourierDiffMatrix2(N)
	default:
		return MatPow(FourierDiffMatrix(N), m)
	}
}

// DifferentiateNodal applies a differentiation matrix D to nodal values u,
// returning D*u.
func DifferentiateNodal(D [][]float64, u []float64) []float64 {
	return MatVec(D, u)
}
