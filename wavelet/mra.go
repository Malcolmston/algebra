package wavelet

import "errors"

// Decomposition is the result of a multi-level (multiresolution) discrete
// wavelet transform. It stores the coarsest approximation together with the
// detail coefficients captured at each level.
type Decomposition struct {
	// Wavelet is the wavelet used to produce and reconstruct the coefficients.
	Wavelet Wavelet
	// Approx holds the coarsest-scale approximation coefficients.
	Approx []float64
	// Details holds the detail coefficients per level, with Details[0] the
	// finest (first) level and Details[len-1] the coarsest level.
	Details [][]float64
}

// MaxLevel returns the maximum number of dyadic decomposition levels that a
// signal of length n admits under periodic boundary handling, which is the
// 2-adic valuation of n (the number of times n can be halved while remaining
// even). It returns 0 for non-positive or odd n.
func MaxLevel(n int) int {
	return wavelet2Valuation(n)
}

// WaveDec computes a multi-level forward wavelet transform of signal with
// wavelet w, cascading the low-pass branch for the requested number of levels.
// It returns an error if levels is not positive, if it exceeds [MaxLevel] for
// the signal length, or if any intermediate length is not even.
func WaveDec(signal []float64, w Wavelet, levels int) (*Decomposition, error) {
	if levels < 1 {
		return nil, errors.New("wavelet: WaveDec requires levels >= 1")
	}
	if levels > MaxLevel(len(signal)) {
		return nil, errors.New("wavelet: levels exceeds the maximum for this signal length")
	}
	approx := append([]float64(nil), signal...)
	details := make([][]float64, levels)
	for i := 0; i < levels; i++ {
		if len(approx)%2 != 0 {
			return nil, errors.New("wavelet: intermediate length is not even")
		}
		a, d := wavelet1DForward(approx, w.lo, w.hi)
		details[i] = d
		approx = a
	}
	return &Decomposition{Wavelet: w, Approx: approx, Details: details}, nil
}

// WaveRec reconstructs the original signal from a [Decomposition], inverting
// [WaveDec] to within floating-point rounding.
func WaveRec(d *Decomposition) []float64 {
	return d.Reconstruct()
}

// Levels returns the number of decomposition levels stored in d.
func (d *Decomposition) Levels() int { return len(d.Details) }

// DetailAt returns the detail coefficients captured at the given level, where
// level 1 is the finest scale and d.Levels() is the coarsest. It panics if
// level is out of range.
func (d *Decomposition) DetailAt(level int) []float64 {
	if level < 1 || level > len(d.Details) {
		panic("wavelet: DetailAt level out of range")
	}
	return d.Details[level-1]
}

// Reconstruct rebuilds the original signal from the stored coefficients by
// cascading the inverse transform from the coarsest to the finest level.
func (d *Decomposition) Reconstruct() []float64 {
	approx := append([]float64(nil), d.Approx...)
	for i := len(d.Details) - 1; i >= 0; i-- {
		approx = wavelet1DInverse(approx, d.Details[i], d.Wavelet.lo, d.Wavelet.hi)
	}
	return approx
}

// Flatten returns all coefficients concatenated into a single slice in the
// order [approx, detail_coarsest, ..., detail_finest], matching the common
// wavelet coefficient layout. The total length equals the original signal
// length.
func (d *Decomposition) Flatten() []float64 {
	out := append([]float64(nil), d.Approx...)
	for i := len(d.Details) - 1; i >= 0; i-- {
		out = append(out, d.Details[i]...)
	}
	return out
}

// Energy returns the total energy (sum of squares) of every stored coefficient.
// For an orthogonal wavelet this equals the energy of the original signal.
func (d *Decomposition) Energy() float64 {
	e := Energy(d.Approx)
	for _, det := range d.Details {
		e += Energy(det)
	}
	return e
}

// ApproxComponent projects the coarsest approximation back to the original
// signal resolution, yielding the smooth multiresolution component of the
// signal. The sum of ApproxComponent and every slice of [DetailComponents]
// reproduces the original signal.
func (d *Decomposition) ApproxComponent() []float64 {
	approx := append([]float64(nil), d.Approx...)
	for i := len(d.Details) - 1; i >= 0; i-- {
		zeros := make([]float64, len(d.Details[i]))
		approx = wavelet1DInverse(approx, zeros, d.Wavelet.lo, d.Wavelet.hi)
	}
	return approx
}

// DetailComponents projects the detail coefficients of each level back to the
// original signal resolution. The returned slice has one full-length component
// per level, with index 0 the finest scale. Together with [ApproxComponent]
// these components form an additive multiresolution decomposition of the signal.
func (d *Decomposition) DetailComponents() [][]float64 {
	levels := len(d.Details)
	out := make([][]float64, levels)
	for lvl := 0; lvl < levels; lvl++ {
		// Reconstruct with a zero approximation and only the level-lvl detail
		// retained; every other band is replaced by zeros of the same length.
		approx := make([]float64, len(d.Approx))
		for i := levels - 1; i >= 0; i-- {
			det := make([]float64, len(d.Details[i]))
			if i == lvl {
				copy(det, d.Details[i])
			}
			approx = wavelet1DInverse(approx, det, d.Wavelet.lo, d.Wavelet.hi)
		}
		out[lvl] = approx
	}
	return out
}

// MRAComponents returns the additive multiresolution components of the signal:
// the full-length smooth approximation and one full-length detail component per
// level (index 0 finest). Their pointwise sum reproduces the original signal to
// within floating-point rounding.
func (d *Decomposition) MRAComponents() (approx []float64, details [][]float64) {
	return d.ApproxComponent(), d.DetailComponents()
}
