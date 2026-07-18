package signal

// Convolve returns the full discrete linear convolution of a and b,
// (a∗b)[k] = Σ a[i]·b[k-i], as a new slice of length len(a)+len(b)-1. If
// either input is empty the result is empty. The inputs are not modified.
func Convolve(a, b []float64) []float64 {
	na, nb := len(a), len(b)
	if na == 0 || nb == 0 {
		return []float64{}
	}
	out := make([]float64, na+nb-1)
	for i := 0; i < na; i++ {
		ai := a[i]
		if ai == 0 {
			continue
		}
		for j := 0; j < nb; j++ {
			out[i+j] += ai * b[j]
		}
	}
	return out
}

// ConvolveSame returns the central part of the full convolution of a and b
// that is the same length as a. It is the full [Convolve] result trimmed by
// (len(b)-1)/2 samples at the front, matching the "same" mode of common
// numerical libraries. The inputs are not modified.
func ConvolveSame(a, b []float64) []float64 {
	na, nb := len(a), len(b)
	if na == 0 || nb == 0 {
		return []float64{}
	}
	full := Convolve(a, b)
	start := (nb - 1) / 2
	out := make([]float64, na)
	copy(out, full[start:start+na])
	return out
}

// ConvolveValid returns only those elements of the convolution of a and b that
// are computed without any zero-padding of the shorter sequence, i.e. where
// the two sequences overlap completely. The result has length
// max(len(a),len(b)) - min(len(a),len(b)) + 1, or is empty if either input is
// empty. The inputs are not modified.
func ConvolveValid(a, b []float64) []float64 {
	na, nb := len(a), len(b)
	if na == 0 || nb == 0 {
		return []float64{}
	}
	full := Convolve(a, b)
	small := nb
	if na < nb {
		small = na
	}
	start := small - 1
	n := na + nb - 1 - 2*(small-1)
	out := make([]float64, n)
	copy(out, full[start:start+n])
	return out
}

// CrossCorrelate returns the full cross-correlation of a and b as a slice of
// length len(a)+len(b)-1. Element index i corresponds to lag ℓ = i-(len(b)-1),
// so out[len(b)-1] is the zero-lag correlation Σ a[n]·b[n]. Positive lags shift
// b to the right relative to a. The inputs are not modified.
func CrossCorrelate(a, b []float64) []float64 {
	na, nb := len(a), len(b)
	if na == 0 || nb == 0 {
		return []float64{}
	}
	out := make([]float64, na+nb-1)
	for i := range out {
		lag := i - (nb - 1)
		var s float64
		for n := 0; n < na; n++ {
			m := n - lag
			if m >= 0 && m < nb {
				s += a[n] * b[m]
			}
		}
		out[i] = s
	}
	return out
}

// AutoCorrelate returns the full autocorrelation of a, equal to
// CrossCorrelate(a, a). The result has length 2·len(a)-1 and is symmetric about
// its centre element, which holds the signal energy Σ a[n]². The input is not
// modified.
func AutoCorrelate(a []float64) []float64 {
	return CrossCorrelate(a, a)
}
