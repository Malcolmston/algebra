package transform

import (
	"math"
	"math/cmplx"
)

// Hilbert returns the analytic signal of the real sequence x, whose real part
// is x and whose imaginary part is the Hilbert transform of x. It is computed
// with the standard FFT method: the negative-frequency components are zeroed
// and the positive-frequency components doubled before inverting. The result
// has the same length as x.
func Hilbert(x []float64) []complex128 {
	n := len(x)
	if n == 0 {
		return []complex128{}
	}
	c := make([]complex128, n)
	for i, v := range x {
		c[i] = complex(v, 0)
	}
	X := FFTAny(c)
	h := make([]float64, n)
	if n%2 == 0 {
		h[0] = 1
		h[n/2] = 1
		for i := 1; i < n/2; i++ {
			h[i] = 2
		}
	} else {
		h[0] = 1
		for i := 1; i < (n+1)/2; i++ {
			h[i] = 2
		}
	}
	for i := range X {
		X[i] *= complex(h[i], 0)
	}
	return IFFTAny(X)
}

// HilbertTransform returns the Hilbert transform of the real sequence x, i.e.
// the imaginary part of its analytic signal (see [Hilbert]). It applies a
// 90-degree phase shift to every frequency component.
func HilbertTransform(x []float64) []float64 {
	a := Hilbert(x)
	out := make([]float64, len(a))
	for i, v := range a {
		out[i] = imag(v)
	}
	return out
}

// Envelope returns the amplitude envelope of the real sequence x, i.e. the
// magnitude of its analytic signal (see [Hilbert]).
func Envelope(x []float64) []float64 {
	a := Hilbert(x)
	out := make([]float64, len(a))
	for i, v := range a {
		out[i] = cmplx.Abs(v)
	}
	return out
}

// InstantaneousPhase returns the unwrapped instantaneous phase, in radians, of
// the real sequence x, obtained as the argument of its analytic signal (see
// [Hilbert]).
func InstantaneousPhase(x []float64) []float64 {
	a := Hilbert(x)
	ph := make([]float64, len(a))
	for i, v := range a {
		ph[i] = cmplx.Phase(v)
	}
	return PhaseUnwrap(ph)
}

// InstantaneousFrequency returns the instantaneous frequency, in the same
// units as sampleRate (typically Hz), of the real sequence x. It is the time
// derivative of the unwrapped [InstantaneousPhase] scaled by
// sampleRate/(2*pi), computed with central differences.
func InstantaneousFrequency(x []float64, sampleRate float64) []float64 {
	ph := InstantaneousPhase(x)
	n := len(ph)
	out := make([]float64, n)
	if n < 2 {
		return out
	}
	scale := sampleRate / (2 * math.Pi)
	out[0] = (ph[1] - ph[0]) * scale
	out[n-1] = (ph[n-1] - ph[n-2]) * scale
	for i := 1; i < n-1; i++ {
		out[i] = (ph[i+1] - ph[i-1]) / 2 * scale
	}
	return out
}

// PhaseUnwrap returns a copy of the phase sequence with discontinuities larger
// than pi removed by adding integer multiples of 2*pi, producing a continuous
// phase curve.
func PhaseUnwrap(phase []float64) []float64 {
	out := make([]float64, len(phase))
	if len(phase) == 0 {
		return out
	}
	out[0] = phase[0]
	offset := 0.0
	for i := 1; i < len(phase); i++ {
		d := phase[i] - phase[i-1]
		for d > math.Pi {
			offset -= 2 * math.Pi
			d -= 2 * math.Pi
		}
		for d < -math.Pi {
			offset += 2 * math.Pi
			d += 2 * math.Pi
		}
		out[i] = phase[i] + offset
	}
	return out
}

// DTFT evaluates the discrete-time Fourier transform of the real sequence x at
// the angular frequency omega (radians/sample), returning
// sum_{n=0}^{N-1} x[n] e^{-i omega n}.
func DTFT(x []float64, omega float64) complex128 {
	var sum complex128
	for n, v := range x {
		sum += complex(v, 0) * cmplx.Rect(1, -omega*float64(n))
	}
	return sum
}

// DTFTComplex evaluates the discrete-time Fourier transform of the complex
// sequence x at the angular frequency omega (radians/sample).
func DTFTComplex(x []complex128, omega float64) complex128 {
	var sum complex128
	for n, v := range x {
		sum += v * cmplx.Rect(1, -omega*float64(n))
	}
	return sum
}

// DTFTSample evaluates the discrete-time Fourier transform of the real
// sequence x at each of the supplied angular frequencies, returning one
// complex value per frequency.
func DTFTSample(x []float64, omegas []float64) []complex128 {
	out := make([]complex128, len(omegas))
	for i, w := range omegas {
		out[i] = DTFT(x, w)
	}
	return out
}

// SampleDTFT samples the discrete-time Fourier transform of the real sequence
// x at m frequencies equally spaced over [0, 2*pi), returning
//
//	X[k] = sum_{n} x[n] e^{-2*pi*i*k*n/m},   k = 0 .. m-1.
//
// When m equals len(x) this coincides with the DFT; larger m interpolates the
// spectrum and smaller m produces the aliased (folded) samples.
func SampleDTFT(x []float64, m int) []complex128 {
	out := make([]complex128, m)
	for k := 0; k < m; k++ {
		omega := 2 * math.Pi * float64(k) / float64(m)
		out[k] = DTFT(x, omega)
	}
	return out
}
