package probability

// Transform returns the distribution of the random variable Y = g(X), applying
// g to every outcome. Outcomes that g maps to the same value have their
// probabilities merged, so the result keeps the canonical sorted-unique form.
func (d Distribution) Transform(g func(float64) float64) Distribution {
	outs := make([]float64, len(d.Outcomes))
	probs := make([]float64, len(d.Probs))
	for i, o := range d.Outcomes {
		outs[i] = g(o)
		probs[i] = d.Probs[i]
	}
	mo, mp := probabilityMerge(outs, probs)
	return Distribution{Outcomes: mo, Probs: mp}
}

// Scale returns the distribution of aX, multiplying every outcome by a. The
// probabilities are unchanged. When a is zero the result is a point mass at zero.
func (d Distribution) Scale(a float64) Distribution {
	return d.Transform(func(x float64) float64 { return a * x })
}

// Shift returns the distribution of X + b, translating every outcome by b. The
// probabilities are unchanged.
func (d Distribution) Shift(b float64) Distribution {
	return d.Transform(func(x float64) float64 { return x + b })
}

// Affine returns the distribution of aX + b, applying a scale followed by a
// shift in a single pass.
func (d Distribution) Affine(a, b float64) Distribution {
	return d.Transform(func(x float64) float64 { return a*x + b })
}

// Standardize returns the distribution of (X - E[X]) / σ, the standardized
// random variable with mean zero and unit variance. It returns an error when the
// standard deviation is zero.
func (d Distribution) Standardize() (Distribution, error) {
	sd := d.StdDev()
	if sd == 0 {
		return Distribution{}, probabilityErrorf("Standardize: zero standard deviation")
	}
	mean := d.Mean()
	return d.Affine(1/sd, -mean/sd), nil
}

// Convolve returns the distribution of the sum X + Y of two independent random
// variables X (the receiver) and Y (the argument): every pair of outcomes is
// added and the joint probabilities Probs[i]·other.Probs[j] accumulated on the
// resulting sums.
func (d Distribution) Convolve(other Distribution) Distribution {
	n := len(d.Outcomes) * len(other.Outcomes)
	outs := make([]float64, 0, n)
	probs := make([]float64, 0, n)
	for i, ox := range d.Outcomes {
		for j, oy := range other.Outcomes {
			outs = append(outs, ox+oy)
			probs = append(probs, d.Probs[i]*other.Probs[j])
		}
	}
	mo, mp := probabilityMerge(outs, probs)
	return Distribution{Outcomes: mo, Probs: mp}
}

// ConvolvePower returns the distribution of the sum of n independent copies of
// X (the receiver): X_1 + X_2 + … + X_n. ConvolvePower(0) is a point mass at
// zero (the empty sum) and ConvolvePower(1) is the receiver itself. It uses
// exponentiation by squaring, so it performs O(log n) convolutions. It returns
// an error for negative n.
func (d Distribution) ConvolvePower(n int) (Distribution, error) {
	if n < 0 {
		return Distribution{}, probabilityErrorf("ConvolvePower: negative n=%d", n)
	}
	result := Distribution{Outcomes: []float64{0}, Probs: []float64{1}} // point mass at 0
	base := d
	for n > 0 {
		if n&1 == 1 {
			result = result.Convolve(base)
		}
		n >>= 1
		if n > 0 {
			base = base.Convolve(base)
		}
	}
	return result, nil
}

// Mixture returns the finite mixture Σ_k weights[k]·components[k], a single
// Distribution whose outcomes are the union of the components' supports with
// probabilities blended according to weights. The weights must be non-negative
// and sum to one within [probabilityTol]. It returns an error if the slice
// lengths differ, either slice is empty, or the weights are invalid.
func Mixture(weights []float64, components []Distribution) (Distribution, error) {
	if len(weights) != len(components) {
		return Distribution{}, probabilityErrorf("Mixture: weights/components length mismatch %d != %d", len(weights), len(components))
	}
	if len(weights) == 0 {
		return Distribution{}, probabilityErrorf("Mixture: empty mixture")
	}
	for k, w := range weights {
		if w < 0 {
			return Distribution{}, probabilityErrorf("Mixture: negative weight %g at index %d", w, k)
		}
	}
	if s := probabilitySum(weights); probabilityAbs(s-1) > probabilityTol {
		return Distribution{}, probabilityErrorf("Mixture: weights sum to %g, not 1", s)
	}
	var outs, probs []float64
	for k, comp := range components {
		for i, o := range comp.Outcomes {
			outs = append(outs, o)
			probs = append(probs, weights[k]*comp.Probs[i])
		}
	}
	mo, mp := probabilityMerge(outs, probs)
	return Distribution{Outcomes: mo, Probs: mp}, nil
}
