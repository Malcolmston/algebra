package transform

import "math"

// Hann returns the n-point Hann (raised-cosine) window,
// w[i] = 0.5*(1 - cos(2*pi*i/(n-1))). The Hann window offers a good balance
// between main-lobe width and side-lobe suppression for spectral analysis.
func Hann(n int) []float64 {
	w := make([]float64, n)
	if n == 1 {
		w[0] = 1
		return w
	}
	for i := 0; i < n; i++ {
		w[i] = 0.5 * (1 - math.Cos(2*math.Pi*float64(i)/float64(n-1)))
	}
	return w
}

// Hamming returns the n-point Hamming window,
// w[i] = 0.54 - 0.46*cos(2*pi*i/(n-1)). It minimizes the nearest side lobe at
// the cost of slower far-off side-lobe roll-off compared with the Hann window.
func Hamming(n int) []float64 {
	w := make([]float64, n)
	if n == 1 {
		w[0] = 1
		return w
	}
	for i := 0; i < n; i++ {
		w[i] = 0.54 - 0.46*math.Cos(2*math.Pi*float64(i)/float64(n-1))
	}
	return w
}

// Blackman returns the n-point Blackman window,
// w[i] = 0.42 - 0.5*cos(2*pi*i/(n-1)) + 0.08*cos(4*pi*i/(n-1)). It provides
// strong side-lobe suppression at the cost of a wider main lobe.
func Blackman(n int) []float64 {
	w := make([]float64, n)
	if n == 1 {
		w[0] = 1
		return w
	}
	for i := 0; i < n; i++ {
		t := 2 * math.Pi * float64(i) / float64(n-1)
		w[i] = 0.42 - 0.5*math.Cos(t) + 0.08*math.Cos(2*t)
	}
	return w
}

// Bartlett returns the n-point Bartlett (triangular) window, which tapers
// linearly to zero at both ends.
func Bartlett(n int) []float64 {
	w := make([]float64, n)
	if n == 1 {
		w[0] = 1
		return w
	}
	half := float64(n-1) / 2
	for i := 0; i < n; i++ {
		w[i] = 1 - math.Abs((float64(i)-half)/half)
	}
	return w
}

// Welch returns the n-point Welch (parabolic) window,
// w[i] = 1 - ((i-(n-1)/2)/((n-1)/2))^2.
func Welch(n int) []float64 {
	w := make([]float64, n)
	if n == 1 {
		w[0] = 1
		return w
	}
	half := float64(n-1) / 2
	for i := 0; i < n; i++ {
		t := (float64(i) - half) / half
		w[i] = 1 - t*t
	}
	return w
}

// ApplyWindow returns the element-wise product of the signal x and the window
// w. It panics if the lengths differ.
func ApplyWindow(x, w []float64) []float64 {
	if len(x) != len(w) {
		panic("transform: ApplyWindow length mismatch")
	}
	out := make([]float64, len(x))
	for i := range x {
		out[i] = x[i] * w[i]
	}
	return out
}

// Periodogram returns the periodogram power-spectral-density estimate of the
// real signal x, defined as (1/N) |DFT(x)[k]|^2 for k = 0 .. N-1.
func Periodogram(x []float64) []float64 {
	n := len(x)
	if n == 0 {
		return []float64{}
	}
	X := FFTReal(x)
	out := make([]float64, n)
	for k := 0; k < n; k++ {
		m := realAbs(X[k])
		out[k] = m * m / float64(n)
	}
	return out
}

// realAbs returns the magnitude of a complex value without importing cmplx.
func realAbs(c complex128) float64 {
	return math.Hypot(real(c), imag(c))
}
