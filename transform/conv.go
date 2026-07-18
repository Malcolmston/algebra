package transform

// transformReverse returns a reversed copy of s.
func transformReverse(s []float64) []float64 {
	n := len(s)
	out := make([]float64, n)
	for i := range s {
		out[n-1-i] = s[i]
	}
	return out
}

// Convolve returns the full linear convolution of the real sequences a and b,
// computed directly. The result has length len(a)+len(b)-1 with
//
//	c[k] = sum_i a[i] * b[k-i].
//
// For long inputs [ConvolveFFT] is asymptotically faster.
func Convolve(a, b []float64) []float64 {
	if len(a) == 0 || len(b) == 0 {
		return []float64{}
	}
	out := make([]float64, len(a)+len(b)-1)
	for i, av := range a {
		for j, bv := range b {
			out[i+j] += av * bv
		}
	}
	return out
}

// ConvolveComplex returns the full linear convolution of the complex sequences
// a and b, computed directly. The result has length len(a)+len(b)-1.
func ConvolveComplex(a, b []complex128) []complex128 {
	if len(a) == 0 || len(b) == 0 {
		return []complex128{}
	}
	out := make([]complex128, len(a)+len(b)-1)
	for i, av := range a {
		for j, bv := range b {
			out[i+j] += av * bv
		}
	}
	return out
}

// ConvolveFFT returns the full linear convolution of the real sequences a and
// b using the FFT (multiplication in the frequency domain). The result has
// length len(a)+len(b)-1 and matches [Convolve] to within floating-point
// rounding; it is much faster for large inputs.
func ConvolveFFT(a, b []float64) []float64 {
	la, lb := len(a), len(b)
	if la == 0 || lb == 0 {
		return []float64{}
	}
	n := la + lb - 1
	L := NextPow2(n)
	ca := make([]complex128, L)
	cb := make([]complex128, L)
	for i, v := range a {
		ca[i] = complex(v, 0)
	}
	for i, v := range b {
		cb[i] = complex(v, 0)
	}
	fa := FFT(ca)
	fb := FFT(cb)
	for i := range fa {
		fa[i] *= fb[i]
	}
	inv := IFFT(fa)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = real(inv[i])
	}
	return out
}

// CircularConvolve returns the circular (cyclic) convolution of a and b. The
// shorter input is zero-extended to the length N of the longer one and the
// result has length N with
//
//	c[k] = sum_{n=0}^{N-1} a[n] * b[(k-n) mod N].
func CircularConvolve(a, b []float64) []float64 {
	n := len(a)
	if len(b) > n {
		n = len(b)
	}
	if n == 0 {
		return []float64{}
	}
	aa := ZeroPadReal(a, n)
	bb := ZeroPadReal(b, n)
	out := make([]float64, n)
	for k := 0; k < n; k++ {
		var sum float64
		for m := 0; m < n; m++ {
			sum += aa[m] * bb[((k-m)%n+n)%n]
		}
		out[k] = sum
	}
	return out
}

// Correlate returns the full cross-correlation of a and b, defined as the
// convolution of a with the time-reversed b. The result has length
// len(a)+len(b)-1; the peak indicates the lag at which the two signals are
// best aligned.
func Correlate(a, b []float64) []float64 {
	return Convolve(a, transformReverse(b))
}

// CrossCorrelateFFT returns the full cross-correlation of a and b computed via
// the FFT. It matches [Correlate] to within floating-point rounding.
func CrossCorrelateFFT(a, b []float64) []float64 {
	return ConvolveFFT(a, transformReverse(b))
}

// AutoCorrelate returns the full autocorrelation of x, i.e. the
// cross-correlation of x with itself. The result has length 2*len(x)-1 and is
// symmetric about its center, whose value equals the signal energy.
func AutoCorrelate(x []float64) []float64 {
	return Correlate(x, x)
}
