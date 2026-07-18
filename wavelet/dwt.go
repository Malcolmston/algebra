package wavelet

// wavelet1DForward applies one level of the periodic discrete wavelet transform
// to x using the analysis filters lo and hi. It returns the approximation and
// detail coefficients, each of length len(x)/2. len(x) must be even.
func wavelet1DForward(x, lo, hi []float64) (approx, detail []float64) {
	n := len(x)
	half := n / 2
	l := len(lo)
	approx = make([]float64, half)
	detail = make([]float64, half)
	for k := 0; k < half; k++ {
		var sa, sd float64
		base := 2 * k
		for j := 0; j < l; j++ {
			idx := (base + j) % n
			v := x[idx]
			sa += lo[j] * v
			sd += hi[j] * v
		}
		approx[k] = sa
		detail[k] = sd
	}
	return approx, detail
}

// wavelet1DInverse reconstructs a signal of length 2*len(approx) from the
// approximation and detail coefficients using the synthesis filters lo and hi.
// It is the exact transpose of wavelet1DForward. len(approx) must equal
// len(detail).
func wavelet1DInverse(approx, detail, lo, hi []float64) []float64 {
	half := len(approx)
	n := 2 * half
	l := len(lo)
	out := make([]float64, n)
	for k := 0; k < half; k++ {
		a := approx[k]
		d := detail[k]
		base := 2 * k
		for j := 0; j < l; j++ {
			idx := (base + j) % n
			out[idx] += lo[j]*a + hi[j]*d
		}
	}
	return out
}

// DWT computes one level of the forward discrete wavelet transform of signal
// using the wavelet w with periodic boundary handling. It returns the
// approximation (low-pass) and detail (high-pass) coefficients, each of length
// len(signal)/2. It panics if len(signal) is odd or smaller than two.
func DWT(signal []float64, w Wavelet) (approx, detail []float64) {
	if len(signal) < 2 || len(signal)%2 != 0 {
		panic("wavelet: DWT requires an even signal length of at least 2")
	}
	return wavelet1DForward(signal, w.lo, w.hi)
}

// IDWT reconstructs a signal from its approximation and detail coefficients
// using the wavelet w, inverting a single level of [DWT]. The returned signal
// has length 2*len(approx). It panics if the two coefficient slices differ in
// length.
func IDWT(approx, detail []float64, w Wavelet) []float64 {
	if len(approx) != len(detail) {
		panic("wavelet: IDWT requires approx and detail of equal length")
	}
	return wavelet1DInverse(approx, detail, w.lo, w.hi)
}

// Coefficients holds the approximation and detail coefficients produced by a
// single level of the discrete wavelet transform.
type Coefficients struct {
	// Approx holds the low-pass approximation coefficients.
	Approx []float64
	// Detail holds the high-pass detail coefficients.
	Detail []float64
}

// Transform computes one level of the forward transform of signal with wavelet
// w and returns the result as a [Coefficients] value. It panics under the same
// conditions as [DWT].
func Transform(signal []float64, w Wavelet) Coefficients {
	a, d := DWT(signal, w)
	return Coefficients{Approx: a, Detail: d}
}

// Inverse reconstructs the signal from the coefficients using wavelet w,
// inverting [Transform].
func (c Coefficients) Inverse(w Wavelet) []float64 {
	return IDWT(c.Approx, c.Detail, w)
}

// Length returns the reconstructed signal length, 2*len(c.Approx).
func (c Coefficients) Length() int { return 2 * len(c.Approx) }

// Energy returns the total energy of the coefficients, the sum of squares of
// the approximation and detail coefficients. For an orthogonal wavelet this
// equals the energy of the original signal.
func (c Coefficients) Energy() float64 {
	return Energy(c.Approx) + Energy(c.Detail)
}
