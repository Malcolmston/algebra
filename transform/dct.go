package transform

import "math"

// DCT computes the unnormalized type-II discrete cosine transform of x:
//
//	X[k] = sum_{n=0}^{N-1} x[n] * cos(pi/N * (n + 1/2) * k),   k = 0 .. N-1.
//
// This is the DCT most commonly meant by "the DCT" and used in signal and
// image compression. It is inverted by [IDCT].
func DCT(x []float64) []float64 {
	n := len(x)
	out := make([]float64, n)
	for k := 0; k < n; k++ {
		var sum float64
		for m := 0; m < n; m++ {
			sum += x[m] * math.Cos(math.Pi/float64(n)*(float64(m)+0.5)*float64(k))
		}
		out[k] = sum
	}
	return out
}

// IDCT computes the inverse of the type-II transform produced by [DCT] (it is
// the scaled type-III transform):
//
//	x[n] = (2/N) * ( X[0]/2 + sum_{k=1}^{N-1} X[k] * cos(pi/N*(n+1/2)*k) ).
//
// Applying IDCT to the output of [DCT] recovers the original signal.
func IDCT(X []float64) []float64 {
	n := len(X)
	out := make([]float64, n)
	for m := 0; m < n; m++ {
		sum := X[0] / 2
		for k := 1; k < n; k++ {
			sum += X[k] * math.Cos(math.Pi/float64(n)*(float64(m)+0.5)*float64(k))
		}
		out[m] = sum * 2 / float64(n)
	}
	return out
}

// DCT1 computes the unnormalized type-I discrete cosine transform of x, which
// requires at least two samples:
//
//	X[k] = (x[0] + (-1)^k x[N-1])/2 + sum_{n=1}^{N-2} x[n] cos(pi*n*k/(N-1)).
//
// The type-I transform is symmetric: applying DCT1 twice and scaling by
// 2/(N-1) recovers the input.
func DCT1(x []float64) []float64 {
	n := len(x)
	out := make([]float64, n)
	if n < 2 {
		copy(out, x)
		return out
	}
	for k := 0; k < n; k++ {
		sum := 0.5 * (x[0] + transformSignPow(k)*x[n-1])
		for m := 1; m < n-1; m++ {
			sum += x[m] * math.Cos(math.Pi*float64(m)*float64(k)/float64(n-1))
		}
		out[k] = sum
	}
	return out
}

// DCT4 computes the unnormalized type-IV discrete cosine transform of x:
//
//	X[k] = sum_{n=0}^{N-1} x[n] cos(pi/N * (n+1/2) * (k+1/2)).
//
// The type-IV transform is its own inverse up to the scale factor 2/N and is
// the building block of the modified DCT used in audio coding.
func DCT4(x []float64) []float64 {
	n := len(x)
	out := make([]float64, n)
	for k := 0; k < n; k++ {
		var sum float64
		for m := 0; m < n; m++ {
			sum += x[m] * math.Cos(math.Pi/float64(n)*(float64(m)+0.5)*(float64(k)+0.5))
		}
		out[k] = sum
	}
	return out
}

// DST computes the unnormalized type-I discrete sine transform of x:
//
//	X[k] = sum_{n=0}^{N-1} x[n] sin(pi*(n+1)*(k+1)/(N+1)),   k = 0 .. N-1.
//
// It is inverted by [IDST].
func DST(x []float64) []float64 {
	n := len(x)
	out := make([]float64, n)
	for k := 0; k < n; k++ {
		var sum float64
		for m := 0; m < n; m++ {
			sum += x[m] * math.Sin(math.Pi*float64(m+1)*float64(k+1)/float64(n+1))
		}
		out[k] = sum
	}
	return out
}

// IDST computes the inverse of the type-I sine transform produced by [DST].
// The type-I DST is orthogonal, so the inverse is the same transform scaled by
// 2/(N+1).
func IDST(X []float64) []float64 {
	n := len(X)
	out := DST(X)
	scale := 2 / float64(n+1)
	for i := range out {
		out[i] *= scale
	}
	return out
}

// transformSignPow returns (-1)^k.
func transformSignPow(k int) float64 {
	if k&1 == 0 {
		return 1
	}
	return -1
}
