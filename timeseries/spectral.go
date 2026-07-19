package timeseries

import (
	"math"
	"math/cmplx"
)

// DFT returns the discrete Fourier transform of the real series x, a complex
// spectrum of the same length computed by direct summation:
// X_k = Σ_t x_t·exp(−2πi·k·t/n).
func DFT(x []float64) []complex128 {
	n := len(x)
	out := make([]complex128, n)
	for k := 0; k < n; k++ {
		var s complex128
		for t := 0; t < n; t++ {
			ang := -2 * math.Pi * float64(k) * float64(t) / float64(n)
			s += complex(x[t], 0) * cmplx.Exp(complex(0, ang))
		}
		out[k] = s
	}
	return out
}

// InverseDFT returns the inverse discrete Fourier transform of the complex
// spectrum X, a complex series of the same length. For a spectrum produced by
// [DFT] the imaginary parts are numerically negligible.
func InverseDFT(X []complex128) []complex128 {
	n := len(X)
	out := make([]complex128, n)
	for t := 0; t < n; t++ {
		var s complex128
		for k := 0; k < n; k++ {
			ang := 2 * math.Pi * float64(k) * float64(t) / float64(n)
			s += X[k] * cmplx.Exp(complex(0, ang))
		}
		out[t] = s / complex(float64(n), 0)
	}
	return out
}

// FourierFrequencies returns the n Fourier frequencies (in cycles per sample
// unit) associated with a length-n DFT: k/(n·d) for k = 0,…,n−1, where d is the
// sample spacing.
func FourierFrequencies(n int, d float64) []float64 {
	out := make([]float64, n)
	if n == 0 || d == 0 {
		return out
	}
	for k := 0; k < n; k++ {
		out[k] = float64(k) / (float64(n) * d)
	}
	return out
}

// Periodogram returns the one-sided raw periodogram of the series, an estimate
// of the power at each Fourier frequency. The returned slices freqs and power
// have length ⌊n/2⌋+1; power[j] = |Σ x_t·e^{−2πi j t/n}|² / n. The mean is
// removed before the transform.
func Periodogram(x []float64) (freqs, power []float64) {
	n := len(x)
	if n == 0 {
		return []float64{}, []float64{}
	}
	m := mean(x)
	w := make([]float64, n)
	for i, v := range x {
		w[i] = v - m
	}
	half := n/2 + 1
	freqs = make([]float64, half)
	power = make([]float64, half)
	for k := 0; k < half; k++ {
		var re, im float64
		for t := 0; t < n; t++ {
			ang := -2 * math.Pi * float64(k) * float64(t) / float64(n)
			re += w[t] * math.Cos(ang)
			im += w[t] * math.Sin(ang)
		}
		freqs[k] = float64(k) / float64(n)
		power[k] = (re*re + im*im) / float64(n)
	}
	return freqs, power
}

// CumulativePeriodogram returns the normalized cumulative periodogram of the
// series over the non-zero Fourier frequencies, a monotone sequence rising from
// near 0 to 1 used to inspect departures from white noise.
func CumulativePeriodogram(x []float64) []float64 {
	_, power := Periodogram(x)
	if len(power) <= 1 {
		return []float64{}
	}
	p := power[1:] // drop the zero frequency
	total := sumf(p)
	out := make([]float64, len(p))
	if total == 0 {
		return out
	}
	var c float64
	for i, v := range p {
		c += v
		out[i] = c / total
	}
	return out
}

// DominantFrequency returns the non-zero Fourier frequency with the largest
// periodogram power, i.e. the strongest cyclical component of the series. It
// returns NaN if the series is too short.
func DominantFrequency(x []float64) float64 {
	freqs, power := Periodogram(x)
	if len(power) <= 1 {
		return math.NaN()
	}
	best := 1
	for k := 2; k < len(power); k++ {
		if power[k] > power[best] {
			best = k
		}
	}
	return freqs[best]
}

// DominantPeriod returns the reciprocal of [DominantFrequency], i.e. the length
// (in samples) of the strongest cycle in the series. It returns NaN if no
// non-zero frequency dominates.
func DominantPeriod(x []float64) float64 {
	f := DominantFrequency(x)
	if math.IsNaN(f) || f == 0 {
		return math.NaN()
	}
	return 1 / f
}

// ARSpectralDensity evaluates the AR power spectral density implied by the
// model at nf equally spaced frequencies in [0, 0.5]. The density at frequency f
// is Sigma2 / |1 − Σ Phi_k·e^{−2πi k f}|². It returns the frequencies and the
// spectral density values.
func ARSpectralDensity(m *ARModel, nf int) (freqs, density []float64) {
	if m == nil || nf < 1 {
		return nil, nil
	}
	freqs = make([]float64, nf)
	density = make([]float64, nf)
	for i := 0; i < nf; i++ {
		f := 0.5 * float64(i) / float64(nf-1)
		if nf == 1 {
			f = 0
		}
		var re, im float64
		re = 1
		for k := 1; k <= m.Order; k++ {
			ang := -2 * math.Pi * f * float64(k)
			re -= m.Phi[k-1] * math.Cos(ang)
			im -= m.Phi[k-1] * math.Sin(ang)
		}
		denom := re*re + im*im
		freqs[i] = f
		if denom == 0 {
			density[i] = math.Inf(1)
		} else {
			density[i] = m.Sigma2 / denom
		}
	}
	return freqs, density
}

// SpectralEntropy returns the normalized spectral (Shannon) entropy of the
// series computed from its periodogram, a value in [0,1] that is near 1 for
// white noise and near 0 for a pure sinusoid. It returns NaN if the series is
// too short.
func SpectralEntropy(x []float64) float64 {
	_, power := Periodogram(x)
	if len(power) <= 1 {
		return math.NaN()
	}
	p := power[1:]
	total := sumf(p)
	if total == 0 {
		return math.NaN()
	}
	var h float64
	for _, v := range p {
		pr := v / total
		if pr > 0 {
			h -= pr * math.Log(pr)
		}
	}
	return h / math.Log(float64(len(p)))
}

// EstimateSeasonalPeriod returns the lag (≥ 2) of the largest peak in the
// autocorrelation function up to maxLag, a data-driven estimate of the dominant
// seasonal period. It returns 0 if no interior peak is found.
func EstimateSeasonalPeriod(x []float64, maxLag int) int {
	if maxLag < 2 || maxLag >= len(x) {
		if len(x) > 2 {
			maxLag = len(x) - 1
		} else {
			return 0
		}
	}
	acf := AutoCorrelation(x, maxLag)
	bestLag := 0
	bestVal := math.Inf(-1)
	for k := 2; k < len(acf)-1; k++ {
		if acf[k] > acf[k-1] && acf[k] >= acf[k+1] && acf[k] > bestVal {
			bestVal = acf[k]
			bestLag = k
		}
	}
	return bestLag
}
