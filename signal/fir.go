package signal

import "math"

// Sinc returns the normalized cardinal sine, sinc(x) = sin(πx)/(πx), with the
// removable singularity sinc(0) = 1 handled exactly. It is the ideal
// continuous-time impulse response of a brick-wall low-pass filter and the
// building block of every windowed-sinc design in this package.
func Sinc(x float64) float64 {
	if x == 0 {
		return 1
	}
	px := math.Pi * x
	return math.Sin(px) / px
}

// signalLowpassKernel builds a length-numtaps windowed-sinc low-pass impulse
// response with the given normalized cutoff (in fractions of Nyquist) and
// window w, then scales it to unit DC gain. It is the shared core of the
// exported low-pass designers.
func signalLowpassKernel(numtaps int, cutoff float64, w []float64) []float64 {
	h := make([]float64, numtaps)
	center := float64(numtaps-1) / 2
	var sum float64
	for n := 0; n < numtaps; n++ {
		x := float64(n) - center
		h[n] = cutoff * Sinc(cutoff*x) * w[n]
		sum += h[n]
	}
	if sum != 0 {
		for n := range h {
			h[n] /= sum
		}
	}
	return h
}

// signalCheckTaps panics if numtaps is not a positive integer, guarding the
// exported designers against degenerate lengths.
func signalCheckTaps(numtaps int) {
	if numtaps < 1 {
		panic("signal: numtaps must be >= 1")
	}
}

// FIRLowpassWin designs a length-numtaps linear-phase low-pass FIR filter with
// the given normalized cutoff (in fractions of Nyquist, 0 < cutoff < 1) using
// the supplied window w, whose length must equal numtaps. The returned taps
// are scaled to unit gain at DC. It panics on a length mismatch or a
// non-positive numtaps.
func FIRLowpassWin(numtaps int, cutoff float64, w []float64) []float64 {
	signalCheckTaps(numtaps)
	if len(w) != numtaps {
		panic("signal: window length must equal numtaps")
	}
	return signalLowpassKernel(numtaps, cutoff, w)
}

// FIRHighpassWin designs a length-numtaps linear-phase high-pass FIR filter
// with the given normalized cutoff using the supplied window w. The design
// uses spectral inversion of a low-pass prototype and therefore requires an
// odd numtaps so that a single centre tap exists; an even numtaps panics, as
// does a window whose length differs from numtaps. The taps have unit gain at
// Nyquist and zero gain at DC.
func FIRHighpassWin(numtaps int, cutoff float64, w []float64) []float64 {
	signalCheckTaps(numtaps)
	if numtaps%2 == 0 {
		panic("signal: high-pass FIR requires odd numtaps")
	}
	if len(w) != numtaps {
		panic("signal: window length must equal numtaps")
	}
	h := signalLowpassKernel(numtaps, cutoff, w)
	for i := range h {
		h[i] = -h[i]
	}
	h[(numtaps-1)/2] += 1
	return h
}

// FIRBandpassWin designs a length-numtaps linear-phase band-pass FIR filter
// that passes normalized frequencies between low and high (both in fractions
// of Nyquist, 0 < low < high < 1) using the supplied window w. It is formed as
// the difference of two unit-DC-gain low-pass prototypes and therefore has
// zero gain at DC and at Nyquist. It panics on a length mismatch or a
// non-positive numtaps.
func FIRBandpassWin(numtaps int, low, high float64, w []float64) []float64 {
	signalCheckTaps(numtaps)
	if len(w) != numtaps {
		panic("signal: window length must equal numtaps")
	}
	lo := signalLowpassKernel(numtaps, low, w)
	hi := signalLowpassKernel(numtaps, high, w)
	h := make([]float64, numtaps)
	for i := range h {
		h[i] = hi[i] - lo[i]
	}
	return h
}

// FIRBandstopWin designs a length-numtaps linear-phase band-stop (notch) FIR
// filter that rejects normalized frequencies between low and high using the
// supplied window w. It is the spectral complement of the corresponding
// band-pass design and requires an odd numtaps so that a single centre tap
// exists; an even numtaps panics, as does a window length mismatch. The taps
// have unit gain at DC and at Nyquist.
func FIRBandstopWin(numtaps int, low, high float64, w []float64) []float64 {
	signalCheckTaps(numtaps)
	if numtaps%2 == 0 {
		panic("signal: band-stop FIR requires odd numtaps")
	}
	if len(w) != numtaps {
		panic("signal: window length must equal numtaps")
	}
	lo := signalLowpassKernel(numtaps, low, w)
	hi := signalLowpassKernel(numtaps, high, w)
	h := make([]float64, numtaps)
	for i := range h {
		h[i] = lo[i] - hi[i]
	}
	h[(numtaps-1)/2] += 1
	return h
}

// FIRLowpass designs a length-numtaps linear-phase low-pass FIR filter with
// the given normalized cutoff (in fractions of Nyquist) using a Hamming
// window. It is the windowed convenience wrapper over [FIRLowpassWin].
func FIRLowpass(numtaps int, cutoff float64) []float64 {
	signalCheckTaps(numtaps)
	return signalLowpassKernel(numtaps, cutoff, Hamming(numtaps))
}

// FIRHighpass designs a length-numtaps linear-phase high-pass FIR filter with
// the given normalized cutoff using a Hamming window. numtaps must be odd; see
// [FIRHighpassWin].
func FIRHighpass(numtaps int, cutoff float64) []float64 {
	return FIRHighpassWin(numtaps, cutoff, Hamming(numtaps))
}

// FIRBandpass designs a length-numtaps linear-phase band-pass FIR filter
// passing normalized frequencies between low and high using a Hamming window.
// See [FIRBandpassWin].
func FIRBandpass(numtaps int, low, high float64) []float64 {
	return FIRBandpassWin(numtaps, low, high, Hamming(numtaps))
}

// FIRBandstop designs a length-numtaps linear-phase band-stop FIR filter
// rejecting normalized frequencies between low and high using a Hamming
// window. numtaps must be odd; see [FIRBandstopWin].
func FIRBandstop(numtaps int, low, high float64) []float64 {
	return FIRBandstopWin(numtaps, low, high, Hamming(numtaps))
}

// ApplyFIR filters the signal x through the FIR filter defined by taps and
// returns the full linear convolution, of length len(x)+len(taps)-1. Because
// the filter is linear-phase, the useful output is delayed by (len(taps)-1)/2
// samples relative to the input. The inputs are not modified.
func ApplyFIR(taps, x []float64) []float64 {
	return Convolve(x, taps)
}

// FIRFrequencyResponse evaluates the complex frequency response H(ω) of the
// FIR filter taps at each normalized frequency in freqs, where a value of 1.0
// corresponds to Nyquist (ω = π rad/sample). H(ω) = Σ h[n]·exp(-jωn). The
// returned slice has the same length as freqs.
func FIRFrequencyResponse(taps, freqs []float64) []complex128 {
	out := make([]complex128, len(freqs))
	for i, f := range freqs {
		w := math.Pi * f
		var re, im float64
		for n, h := range taps {
			re += h * math.Cos(w*float64(n))
			im -= h * math.Sin(w*float64(n))
		}
		out[i] = complex(re, im)
	}
	return out
}

// FIRGroupDelay returns the group delay, in samples, of the FIR filter taps at
// each normalized frequency in freqs (1.0 corresponds to Nyquist). It is
// computed as τ(ω) = Re(Σ n·h[n]·e^{-jωn} / Σ h[n]·e^{-jωn}). For an exactly
// linear-phase filter the result is the constant (len(taps)-1)/2 wherever the
// response does not vanish. Where the response magnitude is negligible the
// group delay is reported as 0.
func FIRGroupDelay(taps, freqs []float64) []float64 {
	out := make([]float64, len(freqs))
	for i, f := range freqs {
		w := math.Pi * f
		var re, im, nre, nim float64
		for n, h := range taps {
			c := math.Cos(w * float64(n))
			s := math.Sin(w * float64(n))
			re += h * c
			im -= h * s
			nre += float64(n) * h * c
			nim -= float64(n) * h * s
		}
		denom := re*re + im*im
		if denom < 1e-30 {
			out[i] = 0
			continue
		}
		// Re( (nre+jnim)/(re+jim) ) = (nre·re + nim·im)/denom.
		out[i] = (nre*re + nim*im) / denom
	}
	return out
}

// FIRFilter is a reusable streaming finite-impulse-response filter. It holds a
// copy of the tap coefficients and an internal delay line so that a long
// signal can be processed block by block with correct state carried across
// calls. The zero value is not usable; construct one with [NewFIRFilter].
type FIRFilter struct {
	taps  []float64
	state []float64 // circular delay line
	pos   int
}

// NewFIRFilter returns a streaming [FIRFilter] initialized with a copy of taps
// and a zeroed delay line. The supplied slice is not retained.
func NewFIRFilter(taps []float64) *FIRFilter {
	t := make([]float64, len(taps))
	copy(t, taps)
	return &FIRFilter{
		taps:  t,
		state: make([]float64, len(taps)),
	}
}

// Reset clears the filter's delay line, returning it to the state of a freshly
// constructed filter without reallocating.
func (f *FIRFilter) Reset() {
	for i := range f.state {
		f.state[i] = 0
	}
	f.pos = 0
}

// ProcessSample feeds a single input sample through the filter and returns the
// corresponding output sample, updating the internal delay line. Successive
// calls implement the causal convolution y[k] = Σ h[n]·x[k-n].
func (f *FIRFilter) ProcessSample(x float64) float64 {
	n := len(f.taps)
	if n == 0 {
		return 0
	}
	f.state[f.pos] = x
	var acc float64
	idx := f.pos
	for i := 0; i < n; i++ {
		acc += f.taps[i] * f.state[idx]
		idx--
		if idx < 0 {
			idx = n - 1
		}
	}
	f.pos++
	if f.pos >= n {
		f.pos = 0
	}
	return acc
}

// Process filters an entire block x through the streaming filter, returning a
// new slice of the same length. State is carried across calls, so feeding the
// concatenation of two blocks yields the same result as feeding them one after
// another. The input is not modified.
func (f *FIRFilter) Process(x []float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = f.ProcessSample(v)
	}
	return out
}
