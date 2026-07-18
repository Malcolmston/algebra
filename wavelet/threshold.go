package wavelet

import "math"

// ThresholdMode selects the shrinkage rule applied to wavelet coefficients.
type ThresholdMode int

const (
	// Soft selects soft thresholding, which shrinks each coefficient toward
	// zero by the threshold and zeroes those below it. It yields a continuous
	// estimator.
	Soft ThresholdMode = iota
	// Hard selects hard thresholding, which keeps coefficients whose magnitude
	// exceeds the threshold unchanged and zeroes the rest.
	Hard
)

// SoftThreshold applies the soft-thresholding (shrinkage) rule to x with
// threshold lambda: sign(x) * max(|x| - lambda, 0). A negative lambda is
// treated as zero.
func SoftThreshold(x, lambda float64) float64 {
	if lambda < 0 {
		lambda = 0
	}
	a := math.Abs(x) - lambda
	if a <= 0 {
		return 0
	}
	return math.Copysign(a, x)
}

// HardThreshold applies the hard-thresholding rule to x with threshold lambda:
// x is kept when |x| > lambda and set to zero otherwise. A negative lambda is
// treated as zero.
func HardThreshold(x, lambda float64) float64 {
	if lambda < 0 {
		lambda = 0
	}
	if math.Abs(x) > lambda {
		return x
	}
	return 0
}

// Threshold applies the thresholding rule selected by mode to x with threshold
// lambda. It panics for an unknown mode.
func Threshold(x, lambda float64, mode ThresholdMode) float64 {
	switch mode {
	case Soft:
		return SoftThreshold(x, lambda)
	case Hard:
		return HardThreshold(x, lambda)
	default:
		panic("wavelet: unknown ThresholdMode")
	}
}

// ThresholdSlice returns a new slice with the thresholding rule selected by
// mode applied elementwise to xs with threshold lambda. The input is not
// modified.
func ThresholdSlice(xs []float64, lambda float64, mode ThresholdMode) []float64 {
	out := make([]float64, len(xs))
	for i, v := range xs {
		out[i] = Threshold(v, lambda, mode)
	}
	return out
}

// UniversalThreshold returns the Donoho-Johnstone universal (VisuShrink)
// threshold sigma * sqrt(2 * ln(n)) for a signal of length n and noise standard
// deviation sigma. It returns 0 for n <= 1.
func UniversalThreshold(n int, sigma float64) float64 {
	if n <= 1 {
		return 0
	}
	return sigma * math.Sqrt(2*math.Log(float64(n)))
}

// EstimateNoiseSigma estimates the standard deviation of white Gaussian noise
// from a set of finest-scale detail coefficients using the robust median
// estimator sigma = median(|detail|) / 0.6745, where 0.6745 is the standard
// normal 75th percentile. It returns 0 for an empty slice.
func EstimateNoiseSigma(detail []float64) float64 {
	if len(detail) == 0 {
		return 0
	}
	abs := make([]float64, len(detail))
	for i, v := range detail {
		abs[i] = math.Abs(v)
	}
	return Median(abs) / 0.6745
}

// Denoise removes noise from signal by wavelet shrinkage (the VisuShrink
// procedure). It performs a levels-deep decomposition with wavelet w, estimates
// the noise level from the finest detail band via [EstimateNoiseSigma],
// computes the [UniversalThreshold], applies the selected thresholding mode to
// every detail band while leaving the approximation untouched, and reconstructs
// the signal. It panics if levels is not valid for the signal length; on a
// noise-free constant signal it returns the input unchanged.
func Denoise(signal []float64, w Wavelet, levels int, mode ThresholdMode) []float64 {
	dec, err := WaveDec(signal, w, levels)
	if err != nil {
		panic("wavelet: Denoise: " + err.Error())
	}
	sigma := EstimateNoiseSigma(dec.Details[0])
	lambda := UniversalThreshold(len(signal), sigma)
	for i := range dec.Details {
		dec.Details[i] = ThresholdSlice(dec.Details[i], lambda, mode)
	}
	return dec.Reconstruct()
}
