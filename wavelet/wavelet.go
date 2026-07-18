package wavelet

import (
	"errors"
	"math"
)

// Wavelet is an orthonormal discrete wavelet, described by its analysis
// (decomposition) and synthesis (reconstruction) filter banks. The zero value
// is not usable; construct one with [Haar], [Daubechies], [DB2] or [DB4].
//
// For the orthogonal families provided here the synthesis filters are equal to
// the analysis filters, because the periodic analysis operator is orthogonal
// and its inverse is its transpose.
type Wavelet struct {
	family string
	vm     int
	lo     []float64 // analysis low-pass (scaling) filter
	hi     []float64 // analysis high-pass (wavelet) filter
}

// waveletQMF returns the high-pass quadrature-mirror filter derived from a
// low-pass filter lo by the relation hi[n] = (-1)^n * lo[L-1-n].
func waveletQMF(lo []float64) []float64 {
	l := len(lo)
	hi := make([]float64, l)
	for n := 0; n < l; n++ {
		s := 1.0
		if n%2 == 1 {
			s = -1.0
		}
		hi[n] = s * lo[l-1-n]
	}
	return hi
}

// waveletFromLowPass builds a Wavelet from a low-pass filter, deriving the
// high-pass filter through the quadrature-mirror relation.
func waveletFromLowPass(family string, vm int, lo []float64) Wavelet {
	cp := make([]float64, len(lo))
	copy(cp, lo)
	return Wavelet{family: family, vm: vm, lo: cp, hi: waveletQMF(cp)}
}

// Haar returns the Haar wavelet, the orthonormal Daubechies wavelet with one
// vanishing moment (equivalently db1). Its scaling filter is
// [1/sqrt(2), 1/sqrt(2)].
func Haar() Wavelet {
	s := 1.0 / math.Sqrt2
	return waveletFromLowPass("haar", 1, []float64{s, s})
}

// daubechiesLowPass holds the low-pass scaling coefficients of the Daubechies
// wavelets, indexed by the number of vanishing moments p. Each filter has 2p
// taps and sums to sqrt(2).
var daubechiesLowPass = map[int][]float64{
	2: {
		(1 + math.Sqrt(3)) / (4 * math.Sqrt2),
		(3 + math.Sqrt(3)) / (4 * math.Sqrt2),
		(3 - math.Sqrt(3)) / (4 * math.Sqrt2),
		(1 - math.Sqrt(3)) / (4 * math.Sqrt2),
	},
	3: {
		0.035226291882100656,
		-0.08544127388224149,
		-0.13501102001039084,
		0.4598775021193313,
		0.8068915093133388,
		0.3326705529509569,
	},
	4: {
		0.23037781330885523,
		0.7148465705525415,
		0.6308807679295904,
		-0.02798376941698385,
		-0.18703481171888114,
		0.030841381835986965,
		0.032883011666982945,
		-0.010597401784997278,
	},
}

// Daubechies returns the Daubechies wavelet with p vanishing moments (a filter
// of 2p taps). Supported orders are p = 1 (which is [Haar]), 2, 3 and 4. It
// returns an error for any other order.
func Daubechies(p int) (Wavelet, error) {
	if p == 1 {
		return Haar(), nil
	}
	lo, ok := daubechiesLowPass[p]
	if !ok {
		return Wavelet{}, errors.New("wavelet: unsupported Daubechies order (supported: 1, 2, 3, 4)")
	}
	name := "db" + string(rune('0'+p))
	return waveletFromLowPass(name, p, lo), nil
}

// DB2 returns the Daubechies wavelet with two vanishing moments (four taps).
func DB2() Wavelet {
	w, _ := Daubechies(2)
	return w
}

// DB4 returns the Daubechies wavelet with four vanishing moments (eight taps).
func DB4() Wavelet {
	w, _ := Daubechies(4)
	return w
}

// Name returns the short name of the wavelet family, such as "haar" or "db4".
func (w Wavelet) Name() string { return w.family }

// Taps returns the number of filter coefficients (the filter length).
func (w Wavelet) Taps() int { return len(w.lo) }

// VanishingMoments returns the number of vanishing moments of the wavelet.
// Polynomial signals of degree below this number produce zero detail
// coefficients away from the boundary.
func (w Wavelet) VanishingMoments() int { return w.vm }

// DecLo returns a copy of the analysis (decomposition) low-pass filter.
func (w Wavelet) DecLo() []float64 { return append([]float64(nil), w.lo...) }

// DecHi returns a copy of the analysis (decomposition) high-pass filter.
func (w Wavelet) DecHi() []float64 { return append([]float64(nil), w.hi...) }

// RecLo returns a copy of the synthesis (reconstruction) low-pass filter. For
// the orthogonal wavelets in this package it equals the analysis low-pass
// filter.
func (w Wavelet) RecLo() []float64 { return append([]float64(nil), w.lo...) }

// RecHi returns a copy of the synthesis (reconstruction) high-pass filter. For
// the orthogonal wavelets in this package it equals the analysis high-pass
// filter.
func (w Wavelet) RecHi() []float64 { return append([]float64(nil), w.hi...) }

// IsOrthogonal reports whether the scaling filter satisfies the orthonormality
// conditions of an orthogonal wavelet to within tol: unit norm and vanishing
// autocorrelation at every non-zero even lag.
func (w Wavelet) IsOrthogonal(tol float64) bool {
	l := len(w.lo)
	// Unit norm.
	var norm float64
	for _, c := range w.lo {
		norm += c * c
	}
	if math.Abs(norm-1) > tol {
		return false
	}
	// Zero autocorrelation at even lags 2, 4, ...
	for m := 2; m < l; m += 2 {
		var s float64
		for n := 0; n+m < l; n++ {
			s += w.lo[n] * w.lo[n+m]
		}
		if math.Abs(s) > tol {
			return false
		}
	}
	return true
}
