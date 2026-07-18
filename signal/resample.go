package signal

import "math"

// Upsample increases the sampling rate of x by an integer factor by inserting
// factor-1 zeros between each input sample (zero stuffing). The result has
// length len(x)·factor, or is the input unchanged when factor <= 1. No
// anti-imaging filter is applied; follow with a low-pass filter or use
// [Interpolate] for a smoothed result. It panics if factor < 1.
func Upsample(x []float64, factor int) []float64 {
	if factor < 1 {
		panic("signal: Upsample factor must be >= 1")
	}
	if factor == 1 {
		out := make([]float64, len(x))
		copy(out, x)
		return out
	}
	out := make([]float64, len(x)*factor)
	for i, v := range x {
		out[i*factor] = v
	}
	return out
}

// Downsample decreases the sampling rate of x by an integer factor by keeping
// every factor-th sample starting from index 0 (decimation without
// filtering). The result has length ceil(len(x)/factor). No anti-aliasing
// filter is applied; use [Decimate] to suppress aliasing first. It panics if
// factor < 1.
func Downsample(x []float64, factor int) []float64 {
	if factor < 1 {
		panic("signal: Downsample factor must be >= 1")
	}
	if factor == 1 {
		out := make([]float64, len(x))
		copy(out, x)
		return out
	}
	n := (len(x) + factor - 1) / factor
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = x[i*factor]
	}
	return out
}

// Decimate reduces the sampling rate of x by an integer factor after applying
// a linear-phase FIR low-pass anti-aliasing filter with cutoff 1/factor of
// Nyquist, then compensating for the filter's group delay so the output stays
// time-aligned with the input. The result has length ceil(len(x)/factor). It
// panics if factor < 1. For factor == 1 the input is returned unchanged.
func Decimate(x []float64, factor int) []float64 {
	if factor < 1 {
		panic("signal: Decimate factor must be >= 1")
	}
	if factor == 1 {
		out := make([]float64, len(x))
		copy(out, x)
		return out
	}
	numtaps := 8*factor + 1 // odd, integer group delay
	taps := FIRLowpass(numtaps, 1.0/float64(factor))
	full := Convolve(x, taps)
	delay := (numtaps - 1) / 2
	// Aligned, same-length filtered signal.
	filt := make([]float64, len(x))
	copy(filt, full[delay:delay+len(x)])
	return Downsample(filt, factor)
}

// Interpolate increases the sampling rate of x by an integer factor,
// zero-stuffing and then applying a linear-phase FIR low-pass anti-imaging
// filter of cutoff 1/factor scaled by factor to preserve amplitude. The result
// has length len(x)·factor and is time-aligned with the input. It panics if
// factor < 1; for factor == 1 the input is returned unchanged.
func Interpolate(x []float64, factor int) []float64 {
	if factor < 1 {
		panic("signal: Interpolate factor must be >= 1")
	}
	if factor == 1 {
		out := make([]float64, len(x))
		copy(out, x)
		return out
	}
	up := Upsample(x, factor)
	numtaps := 8*factor + 1
	taps := FIRLowpass(numtaps, 1.0/float64(factor))
	// Scale so the interpolated samples keep the original amplitude.
	for i := range taps {
		taps[i] *= float64(factor)
	}
	full := Convolve(up, taps)
	delay := (numtaps - 1) / 2
	out := make([]float64, len(up))
	copy(out, full[delay:delay+len(up)])
	return out
}

// ResampleLinear resamples x to exactly numOut samples using linear
// interpolation over the original sample positions. It is a lightweight,
// arbitrary-ratio resampler suitable when a fractional rate change is needed
// and full band-limited interpolation is not required. It returns an empty
// slice when numOut <= 0 or x is empty, and a length-numOut constant slice when
// x has a single sample. The input is not modified.
func ResampleLinear(x []float64, numOut int) []float64 {
	if numOut <= 0 || len(x) == 0 {
		return []float64{}
	}
	out := make([]float64, numOut)
	if len(x) == 1 {
		for i := range out {
			out[i] = x[0]
		}
		return out
	}
	if numOut == 1 {
		out[0] = x[0]
		return out
	}
	step := float64(len(x)-1) / float64(numOut-1)
	for i := 0; i < numOut; i++ {
		pos := float64(i) * step
		lo := int(math.Floor(pos))
		if lo >= len(x)-1 {
			out[i] = x[len(x)-1]
			continue
		}
		frac := pos - float64(lo)
		out[i] = x[lo]*(1-frac) + x[lo+1]*frac
	}
	return out
}
