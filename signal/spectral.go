package signal

import "math"

// DFT returns the discrete Fourier transform of the real signal x, computed
// directly as X[k] = Σ x[n]·exp(-j2πkn/N) for k = 0…N-1. The result has the
// same length as x. The algorithm is O(N²) and intended for correctness and
// modest sizes rather than speed. An empty input yields an empty result.
func DFT(x []float64) []complex128 {
	n := len(x)
	out := make([]complex128, n)
	if n == 0 {
		return out
	}
	for k := 0; k < n; k++ {
		var re, im float64
		for t := 0; t < n; t++ {
			ang := -2 * math.Pi * float64(k) * float64(t) / float64(n)
			re += x[t] * math.Cos(ang)
			im += x[t] * math.Sin(ang)
		}
		out[k] = complex(re, im)
	}
	return out
}

// IDFT returns the inverse discrete Fourier transform of the spectrum X,
// x[n] = (1/N)·Σ X[k]·exp(+j2πkn/N), as a complex slice of the same length. For
// a spectrum produced by [DFT] of a real signal the imaginary parts of the
// result are zero to within rounding. An empty input yields an empty result.
func IDFT(X []complex128) []complex128 {
	n := len(X)
	out := make([]complex128, n)
	if n == 0 {
		return out
	}
	for t := 0; t < n; t++ {
		var re, im float64
		for k := 0; k < n; k++ {
			ang := 2 * math.Pi * float64(k) * float64(t) / float64(n)
			c := math.Cos(ang)
			s := math.Sin(ang)
			xr := real(X[k])
			xi := imag(X[k])
			re += xr*c - xi*s
			im += xr*s + xi*c
		}
		out[t] = complex(re/float64(n), im/float64(n))
	}
	return out
}

// Magnitude returns the element-wise magnitude |X[k]| of a complex spectrum X.
// The result has the same length as X.
func Magnitude(X []complex128) []float64 {
	out := make([]float64, len(X))
	for i, c := range X {
		out[i] = math.Hypot(real(c), imag(c))
	}
	return out
}

// Phase returns the element-wise phase angle arg(X[k]), in radians in the range
// (-π, π], of a complex spectrum X. The result has the same length as X.
func Phase(X []complex128) []float64 {
	out := make([]float64, len(X))
	for i, c := range X {
		out[i] = math.Atan2(imag(c), real(c))
	}
	return out
}

// PowerSpectrum returns the element-wise squared magnitude |X[k]|² of a complex
// spectrum X. The result has the same length as X.
func PowerSpectrum(X []complex128) []float64 {
	out := make([]float64, len(X))
	for i, c := range X {
		out[i] = real(c)*real(c) + imag(c)*imag(c)
	}
	return out
}

// FrequencyBins returns the n bin-centre frequencies of an n-point DFT for a
// signal sampled at fs hertz, i.e. k·fs/n for k = 0…n-1. Bins above n/2
// correspond to negative frequencies fs·(k-n)/n if reinterpreted as a two-sided
// spectrum. An n <= 0 yields an empty slice.
func FrequencyBins(n int, fs float64) []float64 {
	if n <= 0 {
		return []float64{}
	}
	out := make([]float64, n)
	for k := 0; k < n; k++ {
		out[k] = float64(k) * fs / float64(n)
	}
	return out
}

// signalOneSidedPSD converts a full complex spectrum of a length-N real signal
// into a one-sided power spectral density of length N/2+1, using the window
// power normalization u and sampling rate fs. Interior bins are doubled to
// account for the discarded negative-frequency half.
func signalOneSidedPSD(X []complex128, u, fs float64) []float64 {
	n := len(X)
	half := n/2 + 1
	psd := make([]float64, half)
	scale := 1.0 / (fs * u)
	for k := 0; k < half; k++ {
		p := (real(X[k])*real(X[k]) + imag(X[k])*imag(X[k])) * scale
		if k != 0 && !(n%2 == 0 && k == n/2) {
			p *= 2
		}
		psd[k] = p
	}
	return psd
}

// Periodogram estimates the one-sided power spectral density of the real signal
// x sampled at fs hertz. If window is non-nil it is applied element-wise to x
// first and must have the same length; if it is nil a rectangular window is
// used. It returns the bin frequencies (length len(x)/2+1, from 0 to fs/2) and
// the matching PSD values, normalized by the window power so that the estimate
// is a density in units of power per hertz. An empty x yields two empty slices.
func Periodogram(x []float64, window []float64, fs float64) (freqs, psd []float64) {
	n := len(x)
	if n == 0 {
		return []float64{}, []float64{}
	}
	seg := x
	u := float64(n) // Σ w² for a rectangular window is N.
	if window != nil {
		if len(window) != n {
			panic("signal: Periodogram window length must equal len(x)")
		}
		seg = ApplyWindow(x, window)
		u = 0
		for _, w := range window {
			u += w * w
		}
	}
	X := DFT(seg)
	psd = signalOneSidedPSD(X, u, fs)
	half := n/2 + 1
	freqs = make([]float64, half)
	for k := 0; k < half; k++ {
		freqs[k] = float64(k) * fs / float64(n)
	}
	return freqs, psd
}

// WelchPSD estimates the one-sided power spectral density of x sampled at fs
// hertz using Welch's method: x is split into overlapping segments of length
// segLen with the given overlap (number of shared samples, 0 <= overlap <
// segLen), each segment is windowed with a Hann window, its periodogram is
// computed, and the periodograms are averaged. Averaging reduces the variance
// of the estimate at the cost of frequency resolution. It returns the bin
// frequencies (length segLen/2+1) and the averaged PSD. It panics if segLen < 1
// or overlap is out of range; if x is shorter than one segment a single
// zero-padded segment is used.
func WelchPSD(x []float64, segLen, overlap int, fs float64) (freqs, psd []float64) {
	if segLen < 1 {
		panic("signal: WelchPSD segLen must be >= 1")
	}
	if overlap < 0 || overlap >= segLen {
		panic("signal: WelchPSD overlap must satisfy 0 <= overlap < segLen")
	}
	win := Hann(segLen)
	half := segLen/2 + 1
	psd = make([]float64, half)
	step := segLen - overlap
	count := 0
	for start := 0; start+segLen <= len(x); start += step {
		seg := x[start : start+segLen]
		_, p := Periodogram(seg, win, fs)
		for k := range psd {
			psd[k] += p[k]
		}
		count++
	}
	if count == 0 {
		// Signal shorter than one segment: zero-pad a single segment.
		seg := make([]float64, segLen)
		copy(seg, x)
		_, p := Periodogram(seg, win, fs)
		copy(psd, p)
		count = 1
	} else {
		for k := range psd {
			psd[k] /= float64(count)
		}
	}
	freqs = make([]float64, half)
	for k := 0; k < half; k++ {
		freqs[k] = float64(k) * fs / float64(segLen)
	}
	return freqs, psd
}
