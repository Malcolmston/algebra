package spectralpde

import (
	"math"
	"math/cmplx"
)

// isPowerOfTwo reports whether n is a positive power of two.
func isPowerOfTwo(n int) bool { return n > 0 && n&(n-1) == 0 }

// RealToComplex returns a complex slice whose real parts are x and imaginary
// parts are zero.
func RealToComplex(x []float64) []complex128 {
	out := make([]complex128, len(x))
	for i, v := range x {
		out[i] = complex(v, 0)
	}
	return out
}

// ComplexReal returns the real parts of x.
func ComplexReal(x []complex128) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = real(v)
	}
	return out
}

// ComplexImag returns the imaginary parts of x.
func ComplexImag(x []complex128) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = imag(v)
	}
	return out
}

// ComplexAbs returns the magnitudes of x.
func ComplexAbs(x []complex128) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = cmplx.Abs(v)
	}
	return out
}

// DFT computes the (unnormalized) discrete Fourier transform
// X_k = sum_j x_j exp(-2*pi*i*j*k/N) directly in O(N^2).
func DFT(x []complex128) []complex128 {
	n := len(x)
	out := make([]complex128, n)
	for k := 0; k < n; k++ {
		var s complex128
		for j := 0; j < n; j++ {
			angle := -2 * math.Pi * float64(j) * float64(k) / float64(n)
			s += x[j] * cmplx.Exp(complex(0, angle))
		}
		out[k] = s
	}
	return out
}

// IDFT computes the inverse discrete Fourier transform
// x_j = (1/N) sum_k X_k exp(2*pi*i*j*k/N) directly in O(N^2).
func IDFT(x []complex128) []complex128 {
	n := len(x)
	out := make([]complex128, n)
	for j := 0; j < n; j++ {
		var s complex128
		for k := 0; k < n; k++ {
			angle := 2 * math.Pi * float64(j) * float64(k) / float64(n)
			s += x[k] * cmplx.Exp(complex(0, angle))
		}
		out[j] = s / complex(float64(n), 0)
	}
	return out
}

// FFT computes the discrete Fourier transform. When the length is a power of
// two it uses a radix-2 Cooley-Tukey algorithm; otherwise it falls back to the
// direct DFT. The normalization matches DFT.
func FFT(x []complex128) []complex128 {
	n := len(x)
	if !isPowerOfTwo(n) {
		return DFT(x)
	}
	out := make([]complex128, n)
	copy(out, x)
	fftRadix2(out, false)
	return out
}

// IFFT computes the inverse discrete Fourier transform, matching IDFT.
func IFFT(x []complex128) []complex128 {
	n := len(x)
	if !isPowerOfTwo(n) {
		return IDFT(x)
	}
	out := make([]complex128, n)
	copy(out, x)
	fftRadix2(out, true)
	inv := complex(1/float64(n), 0)
	for i := range out {
		out[i] *= inv
	}
	return out
}

// fftRadix2 performs an in-place iterative radix-2 FFT. If inverse is true the
// sign of the exponent is flipped (normalization is applied by the caller).
func fftRadix2(a []complex128, inverse bool) {
	n := len(a)
	// Bit-reversal permutation.
	for i, j := 1, 0; i < n; i++ {
		bit := n >> 1
		for ; j&bit != 0; bit >>= 1 {
			j ^= bit
		}
		j ^= bit
		if i < j {
			a[i], a[j] = a[j], a[i]
		}
	}
	for length := 2; length <= n; length <<= 1 {
		ang := 2 * math.Pi / float64(length)
		if !inverse {
			ang = -ang
		}
		wlen := cmplx.Exp(complex(0, ang))
		for i := 0; i < n; i += length {
			w := complex(1, 0)
			half := length >> 1
			for k := 0; k < half; k++ {
				u := a[i+k]
				v := a[i+k+half] * w
				a[i+k] = u + v
				a[i+k+half] = u - v
				w *= wlen
			}
		}
	}
}

// FFTReal computes the DFT of a real-valued signal.
func FFTReal(x []float64) []complex128 {
	return FFT(RealToComplex(x))
}

// FourierCoefficients returns the complex Fourier coefficients
// c_k = (1/N) sum_j f_j exp(-i*k*x_j) of the values sampled on the periodic
// Fourier grid, in DFT ordering (k = 0..N-1).
func FourierCoefficients(values []float64) []complex128 {
	n := len(values)
	f := FFT(RealToComplex(values))
	inv := complex(1/float64(n), 0)
	for i := range f {
		f[i] *= inv
	}
	return f
}

// FourierInterpolate evaluates, at x, the band-limited trigonometric
// interpolant of values sampled at the N periodic nodes x_j = 2*pi*j/N.
func FourierInterpolate(values []float64, x float64) float64 {
	n := len(values)
	f := FFT(RealToComplex(values))
	inv := 1 / float64(n)
	var sum complex128
	for k := 0; k < n; k++ {
		ck := f[k] * complex(inv, 0)
		if n%2 == 0 && k == n/2 {
			sum += ck * complex(math.Cos(float64(n/2)*x), 0)
			continue
		}
		kk := float64(k)
		if k > n/2 {
			kk = float64(k - n)
		}
		sum += ck * cmplx.Exp(complex(0, kk*x))
	}
	return real(sum)
}

// FourierDifferentiateOrder returns the m-th derivative of the trigonometric
// interpolant, evaluated at the periodic grid nodes.
func FourierDifferentiateOrder(values []float64, m int) []float64 {
	n := len(values)
	f := FFT(RealToComplex(values))
	for k := 0; k < n; k++ {
		kk := float64(k)
		if k > n/2 {
			kk = float64(k - n)
		}
		if n%2 == 0 && k == n/2 {
			kk = float64(n / 2)
			if m%2 == 1 {
				f[k] = 0
				continue
			}
		}
		factor := complex(1, 0)
		base := complex(0, kk)
		for p := 0; p < m; p++ {
			factor *= base
		}
		f[k] *= factor
	}
	inv := IFFT(f)
	return ComplexReal(inv)
}

// FourierDifferentiate returns the first derivative of the trigonometric
// interpolant at the grid nodes.
func FourierDifferentiate(values []float64) []float64 {
	return FourierDifferentiateOrder(values, 1)
}

// FourierDifferentiate2 returns the second derivative of the trigonometric
// interpolant at the grid nodes.
func FourierDifferentiate2(values []float64) []float64 {
	return FourierDifferentiateOrder(values, 2)
}

// DCT1 computes the type-I discrete cosine transform of the N+1 input samples
// (N = len(x)-1): X_k = 0.5*(x_0 + (-1)^k x_N) + sum_{n=1}^{N-1} x_n
// cos(pi*n*k/N), for k = 0..N.
func DCT1(x []float64) []float64 {
	N := len(x) - 1
	out := make([]float64, N+1)
	for k := 0; k <= N; k++ {
		s := 0.5 * (x[0] + math.Pow(-1, float64(k))*x[N])
		for n := 1; n < N; n++ {
			s += x[n] * math.Cos(math.Pi*float64(n)*float64(k)/float64(N))
		}
		out[k] = s
	}
	return out
}

// IDCT1 computes the inverse of DCT1, satisfying IDCT1(DCT1(x)) = x.
func IDCT1(x []float64) []float64 {
	N := len(x) - 1
	tmp := DCT1(x)
	scale := 2.0 / float64(N)
	for i := range tmp {
		tmp[i] *= scale
	}
	return tmp
}

// DST1 computes the type-I discrete sine transform of the N-1 interior samples
// (indexing n, k = 1..N-1): X_k = sum_{n=1}^{N-1} x_{n-1} sin(pi*n*k/N).
func DST1(x []float64) []float64 {
	M := len(x)
	N := M + 1
	out := make([]float64, M)
	for k := 1; k <= M; k++ {
		var s float64
		for n := 1; n <= M; n++ {
			s += x[n-1] * math.Sin(math.Pi*float64(n)*float64(k)/float64(N))
		}
		out[k-1] = s
	}
	return out
}

// IDST1 computes the inverse of DST1, satisfying IDST1(DST1(x)) = x.
func IDST1(x []float64) []float64 {
	M := len(x)
	N := M + 1
	tmp := DST1(x)
	scale := 2.0 / float64(N)
	for i := range tmp {
		tmp[i] *= scale
	}
	return tmp
}
