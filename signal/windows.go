package signal

import "math"

// signalCosineSum evaluates a generalized cosine-sum window of length n at
// every index and returns the resulting slice. The coefficients a are applied
// as w[k] = a0 - a1·cos(θ) + a2·cos(2θ) - a3·cos(3θ) + …, with
// θ = 2π·k/(n-1). This is the shared engine behind the Hann, Hamming,
// Blackman, Blackman-Harris, Nuttall and flat-top windows.
func signalCosineSum(n int, a ...float64) []float64 {
	if n <= 0 {
		return []float64{}
	}
	if n == 1 {
		return []float64{1}
	}
	w := make([]float64, n)
	denom := float64(n - 1)
	for k := 0; k < n; k++ {
		theta := 2 * math.Pi * float64(k) / denom
		var sum float64
		sign := 1.0
		for j, c := range a {
			sum += sign * c * math.Cos(float64(j)*theta)
			sign = -sign
		}
		w[k] = sum
	}
	return w
}

// Rectangular returns the length-n rectangular (boxcar) window, whose samples
// are all 1. It applies no tapering and is equivalent to not windowing at all.
// For n <= 0 an empty slice is returned.
func Rectangular(n int) []float64 {
	if n <= 0 {
		return []float64{}
	}
	w := make([]float64, n)
	for i := range w {
		w[i] = 1
	}
	return w
}

// Hann returns the length-n Hann (raised-cosine) window,
// w[k] = 0.5·(1 - cos(2π·k/(n-1))). It has zero-valued endpoints, a −31 dB
// peak side lobe and a side-lobe roll-off of 18 dB/octave, making it a good
// general-purpose choice for spectral analysis. For n <= 0 an empty slice is
// returned; for n == 1 the single sample is 1.
func Hann(n int) []float64 {
	return signalCosineSum(n, 0.5, 0.5)
}

// Hamming returns the length-n Hamming window,
// w[k] = 0.54 - 0.46·cos(2π·k/(n-1)). Its coefficients are optimized to
// cancel the nearest side lobe, giving a −43 dB peak side lobe at the cost of
// non-zero endpoints. For n <= 0 an empty slice is returned.
func Hamming(n int) []float64 {
	return signalCosineSum(n, 0.54, 0.46)
}

// Blackman returns the length-n Blackman window,
// w[k] = 0.42 - 0.5·cos(2π·k/(n-1)) + 0.08·cos(4π·k/(n-1)). The three-term
// design pushes the peak side lobe down to about −58 dB, trading a wider main
// lobe for far lower spectral leakage. For n <= 0 an empty slice is returned.
func Blackman(n int) []float64 {
	return signalCosineSum(n, 0.42, 0.5, 0.08)
}

// BlackmanHarris returns the length-n four-term Blackman-Harris window using
// the coefficients 0.35875, 0.48829, 0.14128 and 0.01168. It reaches a peak
// side lobe near −92 dB, one of the lowest of any fixed window, and is well
// suited to detecting weak tones near strong ones. For n <= 0 an empty slice
// is returned.
func BlackmanHarris(n int) []float64 {
	return signalCosineSum(n, 0.35875, 0.48829, 0.14128, 0.01168)
}

// Nuttall returns the length-n four-term Nuttall window using the coefficients
// 0.355768, 0.487396, 0.144232 and 0.012604. It has a continuous first
// derivative at the endpoints and a peak side lobe near −93 dB. For n <= 0 an
// empty slice is returned.
func Nuttall(n int) []float64 {
	return signalCosineSum(n, 0.355768, 0.487396, 0.144232, 0.012604)
}

// FlatTop returns the length-n flat-top window using the SR785 five-term
// coefficients. Its main lobe is deliberately broadened and flattened so that
// the amplitude of a tone is recovered accurately regardless of where it falls
// between DFT bins, which makes it the standard choice for calibrated
// amplitude measurements. For n <= 0 an empty slice is returned.
func FlatTop(n int) []float64 {
	return signalCosineSum(n,
		0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368)
}

// Bartlett returns the length-n Bartlett (triangular) window with zero-valued
// endpoints, w[k] = 1 - |(k - (n-1)/2)/((n-1)/2)|. It is the simplest tapered
// window and equals the normalized autoconvolution of two rectangular
// windows. For n <= 0 an empty slice is returned; for n == 1 the sample is 1.
func Bartlett(n int) []float64 {
	if n <= 0 {
		return []float64{}
	}
	if n == 1 {
		return []float64{1}
	}
	w := make([]float64, n)
	half := float64(n-1) / 2
	for k := 0; k < n; k++ {
		w[k] = 1 - math.Abs((float64(k)-half)/half)
	}
	return w
}

// Welch returns the length-n Welch (parabolic) window,
// w[k] = 1 - ((k - (n-1)/2)/((n-1)/2))². Its samples trace a downward
// parabola that vanishes at the endpoints. For n <= 0 an empty slice is
// returned; for n == 1 the sample is 1.
func Welch(n int) []float64 {
	if n <= 0 {
		return []float64{}
	}
	if n == 1 {
		return []float64{1}
	}
	w := make([]float64, n)
	half := float64(n-1) / 2
	for k := 0; k < n; k++ {
		t := (float64(k) - half) / half
		w[k] = 1 - t*t
	}
	return w
}

// Cosine returns the length-n cosine (sine) window,
// w[k] = sin(π·k/(n-1)). It is a single half-cycle of a sine and provides a
// gentle taper with a −23 dB peak side lobe. For n <= 0 an empty slice is
// returned; for n == 1 the sample is 1.
func Cosine(n int) []float64 {
	if n <= 0 {
		return []float64{}
	}
	if n == 1 {
		return []float64{1}
	}
	w := make([]float64, n)
	denom := float64(n - 1)
	for k := 0; k < n; k++ {
		w[k] = math.Sin(math.Pi * float64(k) / denom)
	}
	return w
}

// Gaussian returns the length-n Gaussian window with normalized standard
// deviation sigma (0 < sigma <= 0.5 is typical),
// w[k] = exp(-½·((k - (n-1)/2)/(sigma·(n-1)/2))²). Smaller sigma yields a
// narrower window and greater tapering. For n <= 0 an empty slice is returned;
// for n == 1 the sample is 1.
func Gaussian(n int, sigma float64) []float64 {
	if n <= 0 {
		return []float64{}
	}
	if n == 1 {
		return []float64{1}
	}
	w := make([]float64, n)
	half := float64(n-1) / 2
	for k := 0; k < n; k++ {
		t := (float64(k) - half) / (sigma * half)
		w[k] = math.Exp(-0.5 * t * t)
	}
	return w
}

// BesselI0 returns the modified Bessel function of the first kind of order
// zero, I0(x), evaluated by its always-convergent power series
// Σ ((x/2)^(2m)/(m!)²). It underlies the Kaiser window and is accurate to
// full double precision for the arguments encountered in filter design.
func BesselI0(x float64) float64 {
	sum := 1.0
	term := 1.0
	y := x * x / 4
	for m := 1; m <= 100; m++ {
		term *= y / (float64(m) * float64(m))
		sum += term
		if term < 1e-18*sum {
			break
		}
	}
	return sum
}

// Kaiser returns the length-n Kaiser window with shape parameter beta,
// w[k] = I0(beta·√(1 - (2k/(n-1) - 1)²)) / I0(beta), where I0 is [BesselI0].
// Larger beta trades a wider main lobe for lower side lobes; beta = 0 reduces
// to the rectangular window. For n <= 0 an empty slice is returned; for n == 1
// the sample is 1.
func Kaiser(n int, beta float64) []float64 {
	if n <= 0 {
		return []float64{}
	}
	if n == 1 {
		return []float64{1}
	}
	w := make([]float64, n)
	denom := float64(n - 1)
	i0beta := BesselI0(beta)
	for k := 0; k < n; k++ {
		r := 2*float64(k)/denom - 1
		w[k] = BesselI0(beta*math.Sqrt(1-r*r)) / i0beta
	}
	return w
}

// KaiserBeta returns the Kaiser shape parameter beta that achieves a stop-band
// attenuation of atten decibels, using Kaiser's empirical formula. For
// atten > 50 it returns 0.1102·(atten-8.7); for 21 <= atten <= 50 it returns
// 0.5842·(atten-21)^0.4 + 0.07886·(atten-21); and for atten < 21 it returns 0.
func KaiserBeta(atten float64) float64 {
	switch {
	case atten > 50:
		return 0.1102 * (atten - 8.7)
	case atten >= 21:
		a := atten - 21
		return 0.5842*math.Pow(a, 0.4) + 0.07886*a
	default:
		return 0
	}
}

// Tukey returns the length-n Tukey (tapered-cosine) window with taper fraction
// alpha in [0, 1]. A fraction alpha of the window is a cosine taper split
// between the two ends and the remaining centre is flat. alpha = 0 gives the
// rectangular window and alpha = 1 gives the Hann window. Values of alpha
// outside [0, 1] are clamped. For n <= 0 an empty slice is returned; for
// n == 1 the sample is 1.
func Tukey(n int, alpha float64) []float64 {
	if n <= 0 {
		return []float64{}
	}
	if n == 1 {
		return []float64{1}
	}
	if alpha <= 0 {
		return Rectangular(n)
	}
	if alpha > 1 {
		alpha = 1
	}
	w := make([]float64, n)
	denom := float64(n - 1)
	for k := 0; k < n; k++ {
		x := float64(k) / denom
		switch {
		case x < alpha/2:
			w[k] = 0.5 * (1 + math.Cos(math.Pi*(2*x/alpha-1)))
		case x <= 1-alpha/2:
			w[k] = 1
		default:
			w[k] = 0.5 * (1 + math.Cos(math.Pi*(2*x/alpha-2/alpha+1)))
		}
	}
	return w
}

// ApplyWindow returns the element-wise product of the signal x and the window
// w. The two slices must have the same length; if they do not, ApplyWindow
// panics. The input x is not modified.
func ApplyWindow(x, w []float64) []float64 {
	if len(x) != len(w) {
		panic("signal: ApplyWindow length mismatch")
	}
	out := make([]float64, len(x))
	for i := range x {
		out[i] = x[i] * w[i]
	}
	return out
}
