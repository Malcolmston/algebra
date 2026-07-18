package signal

import "math"

// Biquad is a second-order recursive (IIR) filter section realising the
// transfer function
//
//	H(z) = (B0 + B1·z⁻¹ + B2·z⁻²) / (1 + A1·z⁻¹ + A2·z⁻²).
//
// The leading denominator coefficient is normalized to 1. The Z1 and Z2 fields
// hold the two state variables of a transposed direct-form II realisation and
// are updated by [Biquad.Process]; call [Biquad.Reset] to clear them. The
// design constructors in this package (for example [BiquadLowpass] and
// [BiquadPeaking]) return fully-populated, ready-to-run sections.
type Biquad struct {
	B0, B1, B2 float64 // numerator (feed-forward) coefficients
	A1, A2     float64 // denominator (feedback) coefficients, A0 == 1
	Z1, Z2     float64 // transposed direct-form II state
}

// Reset clears the section's internal state variables so that subsequent
// filtering starts from rest.
func (bq *Biquad) Reset() {
	bq.Z1, bq.Z2 = 0, 0
}

// Process filters one input sample x through the section and returns the output
// sample, advancing the internal state. It implements the transposed
// direct-form II recurrence, which is numerically well behaved for
// fixed-frequency sections.
func (bq *Biquad) Process(x float64) float64 {
	y := bq.B0*x + bq.Z1
	bq.Z1 = bq.B1*x - bq.A1*y + bq.Z2
	bq.Z2 = bq.B2*x - bq.A2*y
	return y
}

// ProcessBlock filters an entire block x through the section, returning a new
// slice of the same length, with state carried across calls. The input is not
// modified.
func (bq *Biquad) ProcessBlock(x []float64) []float64 {
	out := make([]float64, len(x))
	for i, v := range x {
		out[i] = bq.Process(v)
	}
	return out
}

// FrequencyResponse evaluates the section's complex frequency response H(ω) at
// each normalized frequency in freqs, where 1.0 corresponds to Nyquist
// (ω = π rad/sample). The returned slice has the same length as freqs and does
// not depend on the state variables.
func (bq *Biquad) FrequencyResponse(freqs []float64) []complex128 {
	out := make([]complex128, len(freqs))
	for i, f := range freqs {
		w := math.Pi * f
		e1 := complex(math.Cos(-w), math.Sin(-w))
		e2 := complex(math.Cos(-2*w), math.Sin(-2*w))
		num := complex(bq.B0, 0) + complex(bq.B1, 0)*e1 + complex(bq.B2, 0)*e2
		den := complex(1, 0) + complex(bq.A1, 0)*e1 + complex(bq.A2, 0)*e2
		out[i] = num / den
	}
	return out
}

// signalRBJ normalizes the raw analog-cookbook coefficients (before dividing
// by a0) into a ready-to-use Biquad with unit A0.
func signalRBJ(b0, b1, b2, a0, a1, a2 float64) Biquad {
	return Biquad{
		B0: b0 / a0, B1: b1 / a0, B2: b2 / a0,
		A1: a1 / a0, A2: a2 / a0,
	}
}

// signalOmega returns the pre-warped angular frequency w0 = 2π·f0/fs together
// with its cosine and sine, the common preamble of the RBJ cookbook designs.
func signalOmega(f0, fs float64) (w0, cw, sw float64) {
	w0 = 2 * math.Pi * f0 / fs
	return w0, math.Cos(w0), math.Sin(w0)
}

// BiquadLowpass returns a second-order low-pass [Biquad] with corner frequency
// f0 hertz, sampling rate fs hertz and quality factor q, following the RBJ
// audio-EQ cookbook. A q of 1/√2 ≈ 0.7071 gives a maximally-flat (Butterworth)
// response. The DC gain is unity.
func BiquadLowpass(f0, fs, q float64) Biquad {
	_, cw, sw := signalOmega(f0, fs)
	alpha := sw / (2 * q)
	b0 := (1 - cw) / 2
	b1 := 1 - cw
	b2 := (1 - cw) / 2
	a0 := 1 + alpha
	a1 := -2 * cw
	a2 := 1 - alpha
	return signalRBJ(b0, b1, b2, a0, a1, a2)
}

// BiquadHighpass returns a second-order high-pass [Biquad] with corner
// frequency f0 hertz, sampling rate fs hertz and quality factor q, following
// the RBJ cookbook. The gain at Nyquist is unity and the gain at DC is zero.
func BiquadHighpass(f0, fs, q float64) Biquad {
	_, cw, sw := signalOmega(f0, fs)
	alpha := sw / (2 * q)
	b0 := (1 + cw) / 2
	b1 := -(1 + cw)
	b2 := (1 + cw) / 2
	a0 := 1 + alpha
	a1 := -2 * cw
	a2 := 1 - alpha
	return signalRBJ(b0, b1, b2, a0, a1, a2)
}

// BiquadBandpass returns a second-order band-pass [Biquad] centred on f0 hertz
// with sampling rate fs hertz and quality factor q, normalized to unit
// (0 dB) peak gain at the centre frequency, following the RBJ cookbook. Larger
// q gives a narrower pass band.
func BiquadBandpass(f0, fs, q float64) Biquad {
	_, cw, sw := signalOmega(f0, fs)
	alpha := sw / (2 * q)
	b0 := alpha
	b1 := 0.0
	b2 := -alpha
	a0 := 1 + alpha
	a1 := -2 * cw
	a2 := 1 - alpha
	return signalRBJ(b0, b1, b2, a0, a1, a2)
}

// BiquadNotch returns a second-order band-reject (notch) [Biquad] centred on
// f0 hertz with sampling rate fs hertz and quality factor q, following the RBJ
// cookbook. The gain is zero at f0 and unity far from it; larger q gives a
// narrower notch.
func BiquadNotch(f0, fs, q float64) Biquad {
	_, cw, sw := signalOmega(f0, fs)
	alpha := sw / (2 * q)
	b0 := 1.0
	b1 := -2 * cw
	b2 := 1.0
	a0 := 1 + alpha
	a1 := -2 * cw
	a2 := 1 - alpha
	return signalRBJ(b0, b1, b2, a0, a1, a2)
}

// BiquadAllpass returns a second-order all-pass [Biquad] centred on f0 hertz
// with sampling rate fs hertz and quality factor q, following the RBJ
// cookbook. Its magnitude response is unity at every frequency while its phase
// varies, which is useful for phase equalisation and fractional delay.
func BiquadAllpass(f0, fs, q float64) Biquad {
	_, cw, sw := signalOmega(f0, fs)
	alpha := sw / (2 * q)
	b0 := 1 - alpha
	b1 := -2 * cw
	b2 := 1 + alpha
	a0 := 1 + alpha
	a1 := -2 * cw
	a2 := 1 - alpha
	return signalRBJ(b0, b1, b2, a0, a1, a2)
}

// BiquadPeaking returns a second-order peaking-EQ [Biquad] centred on f0 hertz
// with sampling rate fs hertz, quality factor q and a boost or cut of gainDB
// decibels (positive boosts, negative cuts), following the RBJ cookbook. The
// gain far from f0 is unity.
func BiquadPeaking(f0, fs, q, gainDB float64) Biquad {
	_, cw, sw := signalOmega(f0, fs)
	alpha := sw / (2 * q)
	A := math.Pow(10, gainDB/40)
	b0 := 1 + alpha*A
	b1 := -2 * cw
	b2 := 1 - alpha*A
	a0 := 1 + alpha/A
	a1 := -2 * cw
	a2 := 1 - alpha/A
	return signalRBJ(b0, b1, b2, a0, a1, a2)
}

// BiquadLowShelf returns a second-order low-shelf [Biquad] with shelf midpoint
// f0 hertz, sampling rate fs hertz, quality factor q and shelf gain gainDB
// decibels, following the RBJ cookbook. Frequencies below f0 are scaled by the
// shelf gain while high frequencies pass at unity.
func BiquadLowShelf(f0, fs, q, gainDB float64) Biquad {
	_, cw, sw := signalOmega(f0, fs)
	alpha := sw / (2 * q)
	A := math.Pow(10, gainDB/40)
	tsa := 2 * math.Sqrt(A) * alpha
	b0 := A * ((A + 1) - (A-1)*cw + tsa)
	b1 := 2 * A * ((A - 1) - (A+1)*cw)
	b2 := A * ((A + 1) - (A-1)*cw - tsa)
	a0 := (A + 1) + (A-1)*cw + tsa
	a1 := -2 * ((A - 1) + (A+1)*cw)
	a2 := (A + 1) + (A-1)*cw - tsa
	return signalRBJ(b0, b1, b2, a0, a1, a2)
}

// BiquadHighShelf returns a second-order high-shelf [Biquad] with shelf
// midpoint f0 hertz, sampling rate fs hertz, quality factor q and shelf gain
// gainDB decibels, following the RBJ cookbook. Frequencies above f0 are scaled
// by the shelf gain while low frequencies pass at unity.
func BiquadHighShelf(f0, fs, q, gainDB float64) Biquad {
	_, cw, sw := signalOmega(f0, fs)
	alpha := sw / (2 * q)
	A := math.Pow(10, gainDB/40)
	tsa := 2 * math.Sqrt(A) * alpha
	b0 := A * ((A + 1) + (A-1)*cw + tsa)
	b1 := -2 * A * ((A - 1) + (A+1)*cw)
	b2 := A * ((A + 1) + (A-1)*cw - tsa)
	a0 := (A + 1) - (A-1)*cw + tsa
	a1 := 2 * ((A - 1) - (A+1)*cw)
	a2 := (A + 1) - (A-1)*cw - tsa
	return signalRBJ(b0, b1, b2, a0, a1, a2)
}

// signalButterQ returns the quality factors of the second-order sections of an
// order-N Butterworth filter. For each complex-conjugate pole pair the section
// Q is 1/(2·|cos θ|) with θ = π/2 + (2k+1)π/(2N). Odd orders contribute one
// additional real pole that is handled separately by the designers.
func signalButterQ(order int) []float64 {
	qs := make([]float64, 0, order/2)
	for k := 0; k < order/2; k++ {
		theta := math.Pi/2 + float64(2*k+1)*math.Pi/(2*float64(order))
		qs = append(qs, 1/(2*math.Abs(math.Cos(theta))))
	}
	return qs
}

// signalFirstOrderLowpass returns the bilinear-transformed first-order
// low-pass section (as a Biquad with B2 = A2 = 0) with corner f0 and rate fs,
// used for the real pole of odd-order Butterworth low-pass filters.
func signalFirstOrderLowpass(f0, fs float64) Biquad {
	k := math.Tan(math.Pi * f0 / fs)
	norm := 1 / (1 + k)
	return Biquad{B0: k * norm, B1: k * norm, B2: 0, A1: (k - 1) * norm, A2: 0}
}

// signalFirstOrderHighpass returns the bilinear-transformed first-order
// high-pass section with corner f0 and rate fs, used for the real pole of
// odd-order Butterworth high-pass filters.
func signalFirstOrderHighpass(f0, fs float64) Biquad {
	k := math.Tan(math.Pi * f0 / fs)
	norm := 1 / (1 + k)
	return Biquad{B0: norm, B1: -norm, B2: 0, A1: (k - 1) * norm, A2: 0}
}

// ButterworthLowpass designs an order-order Butterworth low-pass filter with
// corner frequency f0 hertz and sampling rate fs hertz, returned as a cascade
// of second-order [Biquad] sections (plus one first-order section when order is
// odd). Apply the cascade with [FilterSOS]. It panics if order < 1. The
// response is maximally flat in the pass band and is 3 dB down at f0.
func ButterworthLowpass(order int, f0, fs float64) []Biquad {
	if order < 1 {
		panic("signal: Butterworth order must be >= 1")
	}
	var sos []Biquad
	for _, q := range signalButterQ(order) {
		sos = append(sos, BiquadLowpass(f0, fs, q))
	}
	if order%2 == 1 {
		sos = append(sos, signalFirstOrderLowpass(f0, fs))
	}
	return sos
}

// ButterworthHighpass designs an order-order Butterworth high-pass filter with
// corner frequency f0 hertz and sampling rate fs hertz, returned as a cascade
// of second-order [Biquad] sections (plus one first-order section when order is
// odd). Apply the cascade with [FilterSOS]. It panics if order < 1. The
// response is maximally flat and 3 dB down at f0.
func ButterworthHighpass(order int, f0, fs float64) []Biquad {
	if order < 1 {
		panic("signal: Butterworth order must be >= 1")
	}
	var sos []Biquad
	for _, q := range signalButterQ(order) {
		sos = append(sos, BiquadHighpass(f0, fs, q))
	}
	if order%2 == 1 {
		sos = append(sos, signalFirstOrderHighpass(f0, fs))
	}
	return sos
}

// FilterSOS filters the signal x through a cascade of second-order sections
// sos and returns a new slice of the same length. Each section is applied in
// turn using a fresh copy so that the caller's sections and their state are
// left unmodified. The input is not modified.
func FilterSOS(sos []Biquad, x []float64) []float64 {
	out := make([]float64, len(x))
	copy(out, x)
	for _, s := range sos {
		section := s
		section.Reset()
		out = section.ProcessBlock(out)
	}
	return out
}
