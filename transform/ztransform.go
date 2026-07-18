package transform

import (
	"math"
	"math/cmplx"
)

// ZTransform evaluates the one-sided Z-transform of the real sequence x at the
// complex point z, returning sum_{n=0}^{N-1} x[n] z^{-n}.
func ZTransform(x []float64, z complex128) complex128 {
	c := make([]complex128, len(x))
	for i, v := range x {
		c[i] = complex(v, 0)
	}
	return ZTransformComplex(c, z)
}

// ZTransformComplex evaluates the one-sided Z-transform of the complex
// sequence x at the point z, returning sum_{n=0}^{N-1} x[n] z^{-n}. The sum is
// evaluated by Horner's method for numerical stability.
func ZTransformComplex(x []complex128, z complex128) complex128 {
	n := len(x)
	if n == 0 {
		return 0
	}
	// Horner: (((x[N-1]/z + x[N-2])/z + ...)/z + x[0]).
	var acc complex128
	for i := n - 1; i >= 1; i-- {
		acc = (acc + x[i]) / z
	}
	return acc + x[0]
}

// InverseZTransform recovers the first n samples of a causal sequence from its
// Z-transform X(z), which is supplied as a callable. The samples are obtained
// by numerically evaluating the inversion contour integral on a circle of the
// given radius using m equally spaced points:
//
//	x[j] = (radius^j / m) * sum_{k=0}^{m-1} X(radius*e^{i*theta_k}) * e^{i*j*theta_k}
//
// The radius must lie outside every pole of X(z); m controls accuracy and
// should be noticeably larger than n.
func InverseZTransform(X func(complex128) complex128, n int, radius float64, m int) []complex128 {
	vals := make([]complex128, m)
	for k := 0; k < m; k++ {
		theta := 2 * math.Pi * float64(k) / float64(m)
		z := complex(radius, 0) * cmplx.Rect(1, theta)
		vals[k] = X(z)
	}
	out := make([]complex128, n)
	for j := 0; j < n; j++ {
		var sum complex128
		for k := 0; k < m; k++ {
			theta := 2 * math.Pi * float64(k) / float64(m)
			sum += vals[k] * cmplx.Rect(1, float64(j)*theta)
		}
		out[j] = sum * complex(math.Pow(radius, float64(j))/float64(m), 0)
	}
	return out
}

// transformChirpPow returns w raised to the (possibly fractional) power e,
// computed as exp(e*log(w)).
func transformChirpPow(w complex128, e float64) complex128 {
	return cmplx.Exp(complex(e, 0) * cmplx.Log(w))
}

// ChirpZTransform computes the chirp-z transform of x, sampling the
// Z-transform along a spiral contour. It returns m points
//
//	X[k] = sum_{j=0}^{N-1} x[j] * a^{-j} * w^{j*k},   k = 0 .. m-1,
//
// where a is the complex starting point and w is the ratio between successive
// spiral points. Setting a = 1 and w = exp(-2*pi*i/N) with m = N reproduces
// the DFT. The transform is evaluated with Bluestein's convolution method in
// O((N+m) log(N+m)) time.
func ChirpZTransform(x []complex128, m int, w, a complex128) []complex128 {
	n := len(x)
	if n == 0 || m == 0 {
		return make([]complex128, m)
	}
	L := NextPow2(n + m - 1)
	g := make([]complex128, L)
	for j := 0; j < n; j++ {
		g[j] = x[j] * transformChirpPow(a, -float64(j)) * transformChirpPow(w, float64(j*j)/2)
	}
	h := make([]complex128, L)
	for i := 0; i < m; i++ {
		h[i] = transformChirpPow(w, -float64(i*i)/2)
	}
	for i := 1; i < n; i++ {
		h[L-i] = transformChirpPow(w, -float64(i*i)/2)
	}
	G := FFT(g)
	H := FFT(h)
	for i := range G {
		G[i] *= H[i]
	}
	r := IFFT(G)
	out := make([]complex128, m)
	for k := 0; k < m; k++ {
		out[k] = r[k] * transformChirpPow(w, float64(k*k)/2)
	}
	return out
}

// Bluestein computes the discrete Fourier transform of x for an arbitrary
// length by expressing it as a convolution (the chirp-z transform with a = 1,
// w = exp(-2*pi*i/N) and m = N). It runs in O(n log n) time and returns the
// same result as [DFT].
func Bluestein(x []complex128) []complex128 {
	n := len(x)
	if n == 0 {
		return []complex128{}
	}
	if n == 1 {
		return []complex128{x[0]}
	}
	w := cmplx.Rect(1, -2*math.Pi/float64(n))
	return ChirpZTransform(x, n, w, complex(1, 0))
}
