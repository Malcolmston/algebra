package spectralpde

import "math"

// ChebyshevT evaluates the Chebyshev polynomial of the first kind T_n at x
// using the stable three-term recurrence.
func ChebyshevT(n int, x float64) float64 {
	if n == 0 {
		return 1
	}
	if n == 1 {
		return x
	}
	tm1, t := 1.0, x
	for k := 2; k <= n; k++ {
		tm1, t = t, 2*x*t-tm1
	}
	return t
}

// ChebyshevU evaluates the Chebyshev polynomial of the second kind U_n at x.
func ChebyshevU(n int, x float64) float64 {
	if n < 0 {
		return 0
	}
	if n == 0 {
		return 1
	}
	um1, u := 1.0, 2*x
	for k := 2; k <= n; k++ {
		um1, u = u, 2*x*u-um1
	}
	return u
}

// ChebyshevTDerivative evaluates the derivative of T_n at x, using the
// identity T_n'(x) = n*U_{n-1}(x).
func ChebyshevTDerivative(n int, x float64) float64 {
	if n == 0 {
		return 0
	}
	return float64(n) * ChebyshevU(n-1, x)
}

// ChebyshevTValues evaluates T_n at each point of xs.
func ChebyshevTValues(n int, xs []float64) []float64 {
	out := make([]float64, len(xs))
	for i, x := range xs {
		out[i] = ChebyshevT(n, x)
	}
	return out
}

// ChebyshevRoots returns the n roots of T_n (the Chebyshev-Gauss nodes).
func ChebyshevRoots(n int) []float64 {
	return ChebyshevGaussNodes(n)
}

// ChebyshevVandermonde returns the (len(nodes)) x (degree+1) matrix whose
// (i, j) entry is T_j(nodes[i]).
func ChebyshevVandermonde(nodes []float64, degree int) [][]float64 {
	m := len(nodes)
	v := Zeros(m, degree+1)
	for i, x := range nodes {
		v[i][0] = 1
		if degree >= 1 {
			v[i][1] = x
		}
		for j := 2; j <= degree; j++ {
			v[i][j] = 2*x*v[i][j-1] - v[i][j-2]
		}
	}
	return v
}

// ClenshavEval evaluates the Chebyshev series sum_k coeffs[k]*T_k(x) using
// Clenshaw's recurrence.
//
// Deprecated: use ClenshawEval; this misspelled alias is retained for
// backward compatibility.
func ClenshavEval(coeffs []float64, x float64) float64 {
	return ClenshawEval(coeffs, x)
}

// ClenshawEval evaluates the Chebyshev series sum_k coeffs[k]*T_k(x) using
// Clenshaw's backward recurrence, which is numerically stable on [-1, 1].
func ClenshawEval(coeffs []float64, x float64) float64 {
	n := len(coeffs)
	if n == 0 {
		return 0
	}
	var bk1, bk2 float64
	for k := n - 1; k >= 1; k-- {
		bk1, bk2 = 2*x*bk1-bk2+coeffs[k], bk1
	}
	return x*bk1 - bk2 + coeffs[0]
}

// ClenshawEvalDerivative evaluates the derivative of the Chebyshev series
// sum_k coeffs[k]*T_k at x.
func ClenshawEvalDerivative(coeffs []float64, x float64) float64 {
	d := ChebyshevDifferentiateCoeffs(coeffs)
	return ClenshawEval(d, x)
}

// ChebyshevCoefficients returns the discrete Chebyshev transform of the values
// f(x_j) sampled at the N+1 Chebyshev-Gauss-Lobatto nodes, where N =
// len(values)-1. The returned coefficients a_k satisfy f(x) = sum_k a_k T_k(x)
// as the polynomial interpolant through the nodes.
func ChebyshevCoefficients(values []float64) []float64 {
	N := len(values) - 1
	if N < 1 {
		out := make([]float64, len(values))
		copy(out, values)
		return out
	}
	a := make([]float64, N+1)
	for k := 0; k <= N; k++ {
		var s float64
		for j := 0; j <= N; j++ {
			cj := 1.0
			if j == 0 || j == N {
				cj = 2.0
			}
			s += values[j] / cj * math.Cos(math.Pi*float64(j)*float64(k)/float64(N))
		}
		ck := 1.0
		if k == 0 || k == N {
			ck = 2.0
		}
		a[k] = 2.0 / (float64(N) * ck) * s
	}
	return a
}

// ChebyshevValuesFromCoeffs reconstructs the nodal values at the N+1
// Chebyshev-Gauss-Lobatto nodes (N = len(coeffs)-1) from Chebyshev
// coefficients. It is the inverse of ChebyshevCoefficients.
func ChebyshevValuesFromCoeffs(coeffs []float64) []float64 {
	N := len(coeffs) - 1
	if N < 1 {
		out := make([]float64, len(coeffs))
		copy(out, coeffs)
		return out
	}
	x := ChebyshevGaussLobattoNodes(N)
	out := make([]float64, N+1)
	for j := 0; j <= N; j++ {
		out[j] = ClenshawEval(coeffs, x[j])
	}
	return out
}

// ChebyshevDifferentiateCoeffs returns the Chebyshev coefficients of the
// derivative of the series represented by coeffs.
func ChebyshevDifferentiateCoeffs(coeffs []float64) []float64 {
	N := len(coeffs) - 1
	if N < 1 {
		return []float64{0}
	}
	d := make([]float64, N+1)
	for k := N; k >= 1; k-- {
		var next float64
		if k+1 <= N {
			next = d[k+1]
		}
		d[k-1] = next + 2*float64(k)*coeffs[k]
	}
	d[0] /= 2
	return d[:N]
}

// ChebyshevIntegrateCoeffs returns the Chebyshev coefficients of an
// antiderivative of the series represented by coeffs. The integration constant
// (coefficient of T_0) is set to zero.
func ChebyshevIntegrateCoeffs(coeffs []float64) []float64 {
	N := len(coeffs) - 1
	b := make([]float64, N+2)
	for k := 1; k <= N+1; k++ {
		var akm1, akp1 float64
		if k-1 >= 0 && k-1 <= N {
			akm1 = coeffs[k-1]
		}
		if k+1 <= N {
			akp1 = coeffs[k+1]
		}
		b[k] = (akm1 - akp1) / (2 * float64(k))
	}
	return b
}

// ChebyshevInterpolate evaluates, at x, the polynomial interpolant of values
// sampled at the N+1 Chebyshev-Gauss-Lobatto nodes.
func ChebyshevInterpolate(values []float64, x float64) float64 {
	return ClenshawEval(ChebyshevCoefficients(values), x)
}

// ChebyshevFit samples f at the N+1 Chebyshev-Gauss-Lobatto nodes and returns
// the resulting Chebyshev coefficients of the interpolant.
func ChebyshevFit(f func(float64) float64, N int) []float64 {
	x := ChebyshevGaussLobattoNodes(N)
	v := ApplyFunc(f, x)
	return ChebyshevCoefficients(v)
}

// ChebyshevFitInterval samples f at the CGL nodes mapped to [a, b] and returns
// the Chebyshev coefficients of the interpolant in the reference variable.
func ChebyshevFitInterval(f func(float64) float64, N int, a, b float64) []float64 {
	ref := ChebyshevGaussLobattoNodes(N)
	v := make([]float64, N+1)
	for i, xr := range ref {
		v[i] = f(AffineMap(xr, a, b))
	}
	return ChebyshevCoefficients(v)
}

// ChebyshevEvalInterval evaluates a reference-variable Chebyshev series at the
// physical point y in [a, b].
func ChebyshevEvalInterval(coeffs []float64, y, a, b float64) float64 {
	return ClenshawEval(coeffs, AffineMapInverse(y, a, b))
}

// ChebyshevIntegral returns the definite integral over [-1, 1] of the
// Chebyshev series represented by coeffs. Only even-index terms contribute,
// with integral(T_2m) = -2/(4m^2-1) and integral(T_0) = 2.
func ChebyshevIntegral(coeffs []float64) float64 {
	var s float64
	for k := 0; k < len(coeffs); k += 2 {
		if k == 0 {
			s += coeffs[0] * 2
		} else {
			s += coeffs[k] * (-2.0 / float64(k*k-1))
		}
	}
	return s
}

// ChebyshevTruncationError estimates the interpolation error of a Chebyshev
// series by the sum of absolute values of the tail coefficients beyond index
// keep-1.
func ChebyshevTruncationError(coeffs []float64, keep int) float64 {
	var s float64
	for k := keep; k < len(coeffs); k++ {
		s += math.Abs(coeffs[k])
	}
	return s
}
