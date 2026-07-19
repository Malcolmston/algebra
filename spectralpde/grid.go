package spectralpde

import "math"

// ChebyshevGaussLobattoNodes returns the N+1 Chebyshev-Gauss-Lobatto nodes
// x_j = cos(pi*j/N) on [-1, 1], ordered from +1 (j=0) down to -1 (j=N). N must
// be >= 1.
func ChebyshevGaussLobattoNodes(N int) []float64 {
	x := make([]float64, N+1)
	for j := 0; j <= N; j++ {
		x[j] = math.Cos(math.Pi * float64(j) / float64(N))
	}
	return x
}

// ChebyshevGaussLobattoAngles returns the angles theta_j = pi*j/N whose
// cosines are the Chebyshev-Gauss-Lobatto nodes.
func ChebyshevGaussLobattoAngles(N int) []float64 {
	t := make([]float64, N+1)
	for j := 0; j <= N; j++ {
		t[j] = math.Pi * float64(j) / float64(N)
	}
	return t
}

// ChebyshevGaussLobattoWeights returns the quadrature weights associated with
// the Chebyshev-Gauss-Lobatto nodes for the weight function 1/sqrt(1-x^2).
// The weights are pi/N for interior nodes and pi/(2N) at the endpoints.
func ChebyshevGaussLobattoWeights(N int) []float64 {
	w := make([]float64, N+1)
	for j := 0; j <= N; j++ {
		if j == 0 || j == N {
			w[j] = math.Pi / (2 * float64(N))
		} else {
			w[j] = math.Pi / float64(N)
		}
	}
	return w
}

// ChebyshevGaussNodes returns the N Chebyshev-Gauss nodes, i.e. the roots of
// the Chebyshev polynomial T_N, x_k = cos(pi*(2k+1)/(2N)) for k = 0..N-1.
func ChebyshevGaussNodes(N int) []float64 {
	x := make([]float64, N)
	for k := 0; k < N; k++ {
		x[k] = math.Cos(math.Pi * float64(2*k+1) / float64(2*N))
	}
	return x
}

// ChebyshevGaussWeights returns the Chebyshev-Gauss quadrature weights, all
// equal to pi/N, for the weight function 1/sqrt(1-x^2).
func ChebyshevGaussWeights(N int) []float64 {
	w := make([]float64, N)
	for k := 0; k < N; k++ {
		w[k] = math.Pi / float64(N)
	}
	return w
}

// FourierNodes returns the N equispaced periodic nodes x_j = 2*pi*j/N on
// [0, 2*pi) for j = 0..N-1.
func FourierNodes(N int) []float64 {
	x := make([]float64, N)
	for j := 0; j < N; j++ {
		x[j] = 2 * math.Pi * float64(j) / float64(N)
	}
	return x
}

// FourierNodesInterval returns N equispaced periodic nodes on [a, b).
func FourierNodesInterval(N int, a, b float64) []float64 {
	x := make([]float64, N)
	h := (b - a) / float64(N)
	for j := 0; j < N; j++ {
		x[j] = a + float64(j)*h
	}
	return x
}

// FourierWeights returns the trapezoidal (spectrally accurate for periodic
// functions) quadrature weights 2*pi/N on the Fourier grid.
func FourierWeights(N int) []float64 {
	return VectorFill(N, 2*math.Pi/float64(N))
}

// FourierWavenumbers returns the signed wavenumbers associated with the
// length-N periodic FFT ordering: 0, 1, ..., N/2-1, -N/2, ..., -1 for even N.
func FourierWavenumbers(N int) []float64 {
	k := make([]float64, N)
	for i := 0; i < N; i++ {
		if i <= N/2 {
			k[i] = float64(i)
		} else {
			k[i] = float64(i - N)
		}
	}
	if N%2 == 0 {
		k[N/2] = float64(N / 2)
	}
	return k
}

// LinSpace returns n equally spaced points from a to b inclusive. For n == 1
// it returns {a}.
func LinSpace(a, b float64, n int) []float64 {
	if n <= 0 {
		return []float64{}
	}
	if n == 1 {
		return []float64{a}
	}
	x := make([]float64, n)
	h := (b - a) / float64(n-1)
	for i := 0; i < n; i++ {
		x[i] = a + float64(i)*h
	}
	x[n-1] = b
	return x
}

// UniformNodes returns n+1 equally spaced nodes on [a, b].
func UniformNodes(n int, a, b float64) []float64 {
	return LinSpace(a, b, n+1)
}

// AffineMap returns the affine transformation of x from the reference
// interval [-1, 1] to the interval [a, b].
func AffineMap(x, a, b float64) float64 {
	return 0.5*(a+b) + 0.5*(b-a)*x
}

// AffineMapInverse maps y from [a, b] back to the reference interval [-1, 1].
func AffineMapInverse(y, a, b float64) float64 {
	return (2*y - (a + b)) / (b - a)
}

// MapToInterval maps a slice of reference nodes on [-1, 1] to [a, b].
func MapToInterval(nodes []float64, a, b float64) []float64 {
	out := make([]float64, len(nodes))
	for i, x := range nodes {
		out[i] = AffineMap(x, a, b)
	}
	return out
}

// MapFromInterval maps a slice of nodes on [a, b] back to [-1, 1].
func MapFromInterval(nodes []float64, a, b float64) []float64 {
	out := make([]float64, len(nodes))
	for i, y := range nodes {
		out[i] = AffineMapInverse(y, a, b)
	}
	return out
}

// IntervalScale returns d/dx of the reference-to-physical map, i.e. the factor
// 2/(b-a) that multiplies reference derivatives to obtain physical ones.
func IntervalScale(a, b float64) float64 {
	return 2 / (b - a)
}

// ChebyshevGaussLobattoNodesInterval returns the CGL nodes mapped to [a, b].
func ChebyshevGaussLobattoNodesInterval(N int, a, b float64) []float64 {
	return MapToInterval(ChebyshevGaussLobattoNodes(N), a, b)
}
