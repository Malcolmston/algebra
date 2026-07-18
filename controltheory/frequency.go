package controltheory

import (
	"math"
	"math/cmplx"
)

// BodePoint is a single sample of a Bode plot at one angular frequency.
type BodePoint struct {
	// Omega is the angular frequency in radians per second.
	Omega float64
	// MagnitudeDB is the magnitude of G(jω) expressed in decibels.
	MagnitudeDB float64
	// PhaseDeg is the phase of G(jω) in degrees.
	PhaseDeg float64
}

// NyquistPoint is a single sample of a Nyquist plot at one angular frequency.
type NyquistPoint struct {
	// Omega is the angular frequency in radians per second.
	Omega float64
	// Real is the real part of G(jω).
	Real float64
	// Imag is the imaginary part of G(jω).
	Imag float64
}

// MagnitudeDB returns 20·log10(|value|) in decibels.
func MagnitudeDB(value complex128) float64 {
	return 20 * math.Log10(cmplx.Abs(value))
}

// PhaseDeg returns the phase angle of value in degrees.
func PhaseDeg(value complex128) float64 {
	return cmplx.Phase(value) * 180 / math.Pi
}

// LogSpace returns n logarithmically spaced values between 10^startExp and
// 10^endExp inclusive. It is convenient for generating frequency grids for
// Bode and Nyquist sampling. It returns a single-element slice when n is 1 and
// an empty slice when n is less than 1.
func LogSpace(startExp, endExp float64, n int) []float64 {
	if n < 1 {
		return nil
	}
	if n == 1 {
		return []float64{math.Pow(10, startExp)}
	}
	out := make([]float64, n)
	step := (endExp - startExp) / float64(n-1)
	for i := 0; i < n; i++ {
		out[i] = math.Pow(10, startExp+step*float64(i))
	}
	return out
}

// Bode returns the Bode-plot samples (magnitude in dB and phase in degrees) of
// the transfer function at each supplied angular frequency. Phase is computed
// with a continuous unwrap so successive samples do not jump by 2π.
func (g TransferFunction) Bode(omegas []float64) []BodePoint {
	out := make([]BodePoint, len(omegas))
	var prev float64
	first := true
	for i, w := range omegas {
		resp := g.FrequencyResponse(w)
		ph := cmplx.Phase(resp)
		if !first {
			ph = controltheoryUnwrap(prev, ph)
		}
		prev = ph
		first = false
		out[i] = BodePoint{
			Omega:       w,
			MagnitudeDB: 20 * math.Log10(cmplx.Abs(resp)),
			PhaseDeg:    ph * 180 / math.Pi,
		}
	}
	return out
}

// controltheoryUnwrap adjusts angle by multiples of 2π so it is within π of the
// previous (already unwrapped) angle.
func controltheoryUnwrap(prev, angle float64) float64 {
	for angle-prev > math.Pi {
		angle -= 2 * math.Pi
	}
	for angle-prev < -math.Pi {
		angle += 2 * math.Pi
	}
	return angle
}

// Nyquist returns the Nyquist-plot samples (real and imaginary parts of G(jω))
// of the transfer function at each supplied angular frequency.
func (g TransferFunction) Nyquist(omegas []float64) []NyquistPoint {
	out := make([]NyquistPoint, len(omegas))
	for i, w := range omegas {
		resp := g.FrequencyResponse(w)
		out[i] = NyquistPoint{Omega: w, Real: real(resp), Imag: imag(resp)}
	}
	return out
}

// GainCrossoverFrequency returns the lowest angular frequency in the supplied
// grid at which the open-loop magnitude |G(jω)| crosses unity (0 dB), refined
// by linear interpolation in the log-magnitude versus log-frequency plane. The
// boolean result reports whether such a crossing was found within the grid.
func (g TransferFunction) GainCrossoverFrequency(omegas []float64) (float64, bool) {
	prevMag := 0.0
	prevW := 0.0
	have := false
	for _, w := range omegas {
		mag := math.Log10(cmplx.Abs(g.FrequencyResponse(w)))
		if have {
			if (prevMag >= 0 && mag <= 0) || (prevMag <= 0 && mag >= 0) {
				if prevMag == mag {
					return w, true
				}
				frac := prevMag / (prevMag - mag)
				lw := math.Log10(prevW) + frac*(math.Log10(w)-math.Log10(prevW))
				return math.Pow(10, lw), true
			}
		}
		prevMag = mag
		prevW = w
		have = true
	}
	return 0, false
}

// PhaseCrossoverFrequency returns the lowest angular frequency in the supplied
// grid at which the open-loop phase crosses -180 degrees, refined by linear
// interpolation. The boolean result reports whether such a crossing was found.
func (g TransferFunction) PhaseCrossoverFrequency(omegas []float64) (float64, bool) {
	var prevPh, prevW float64
	have := false
	var prevUnwrapped float64
	firstPhase := true
	for _, w := range omegas {
		ph := cmplx.Phase(g.FrequencyResponse(w))
		if !firstPhase {
			ph = controltheoryUnwrap(prevUnwrapped, ph)
		}
		prevUnwrapped = ph
		firstPhase = false
		deg := ph * 180 / math.Pi
		if have {
			if (prevPh+180)*(deg+180) <= 0 && prevPh != deg {
				frac := (-180 - prevPh) / (deg - prevPh)
				lw := math.Log10(prevW) + frac*(math.Log10(w)-math.Log10(prevW))
				return math.Pow(10, lw), true
			}
		}
		prevPh = deg
		prevW = w
		have = true
	}
	return 0, false
}

// PhaseMargin returns the phase margin in degrees, defined as 180 plus the
// open-loop phase at the gain crossover frequency, evaluated over the supplied
// frequency grid. The boolean result reports whether a gain crossover was
// found. A positive phase margin indicates a stable closed loop under unity
// negative feedback.
func (g TransferFunction) PhaseMargin(omegas []float64) (float64, bool) {
	wc, ok := g.GainCrossoverFrequency(omegas)
	if !ok {
		return 0, false
	}
	ph := cmplx.Phase(g.FrequencyResponse(wc)) * 180 / math.Pi
	// Normalize phase to (-360, 0] region typical of open-loop responses.
	for ph > 0 {
		ph -= 360
	}
	for ph <= -360 {
		ph += 360
	}
	return 180 + ph, true
}

// GainMargin returns the gain margin in decibels, defined as -20·log10|G(jω)|
// at the phase crossover frequency (where the phase is -180 degrees), evaluated
// over the supplied frequency grid. The boolean result reports whether a phase
// crossover was found. A positive gain margin indicates a stable closed loop
// under unity negative feedback.
func (g TransferFunction) GainMargin(omegas []float64) (float64, bool) {
	wc, ok := g.PhaseCrossoverFrequency(omegas)
	if !ok {
		return 0, false
	}
	mag := cmplx.Abs(g.FrequencyResponse(wc))
	return -20 * math.Log10(mag), true
}
