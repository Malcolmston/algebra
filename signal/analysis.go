package signal

import "math"

// MovingAverage returns the causal (trailing) simple moving average of x with
// window length w: out[i] is the mean of x[max(0,i-w+1) … i]. At the start of
// the signal, where fewer than w samples are available, the average is taken
// over the samples seen so far. The result has the same length as x. It panics
// if w < 1.
func MovingAverage(x []float64, w int) []float64 {
	if w < 1 {
		panic("signal: MovingAverage window must be >= 1")
	}
	out := make([]float64, len(x))
	var sum float64
	for i := range x {
		sum += x[i]
		if i >= w {
			sum -= x[i-w]
		}
		denom := w
		if i+1 < w {
			denom = i + 1
		}
		out[i] = sum / float64(denom)
	}
	return out
}

// MovingAverageCentered returns the centered simple moving average of x with
// window length w. Each output is the mean of the w samples centred on the
// current index, with the window truncated at the signal boundaries so the
// output length matches the input. For an even w the window is biased one
// sample towards the past. It panics if w < 1.
func MovingAverageCentered(x []float64, w int) []float64 {
	if w < 1 {
		panic("signal: MovingAverageCentered window must be >= 1")
	}
	out := make([]float64, len(x))
	half := w / 2
	for i := range x {
		lo := i - half
		hi := i - half + w - 1
		if lo < 0 {
			lo = 0
		}
		if hi > len(x)-1 {
			hi = len(x) - 1
		}
		var sum float64
		for j := lo; j <= hi; j++ {
			sum += x[j]
		}
		out[i] = sum / float64(hi-lo+1)
	}
	return out
}

// ExponentialMovingAverage returns the exponentially-weighted moving average of
// x with smoothing factor alpha in [0, 1]: out[0] = x[0] and
// out[i] = alpha·x[i] + (1-alpha)·out[i-1]. Larger alpha weights recent samples
// more heavily and reacts faster. It panics if alpha is outside [0, 1]. An
// empty input yields an empty result.
func ExponentialMovingAverage(x []float64, alpha float64) []float64 {
	if alpha < 0 || alpha > 1 {
		panic("signal: ExponentialMovingAverage alpha must be in [0,1]")
	}
	out := make([]float64, len(x))
	if len(x) == 0 {
		return out
	}
	out[0] = x[0]
	for i := 1; i < len(x); i++ {
		out[i] = alpha*x[i] + (1-alpha)*out[i-1]
	}
	return out
}

// CumulativeSum returns the running (prefix) sum of x, where out[i] is the sum
// of x[0…i]. The result has the same length as x. The input is not modified.
func CumulativeSum(x []float64) []float64 {
	out := make([]float64, len(x))
	var sum float64
	for i, v := range x {
		sum += v
		out[i] = sum
	}
	return out
}

// Diff returns the first-order forward difference of x,
// out[i] = x[i+1] - x[i]. The result has length len(x)-1, or is empty when x
// has fewer than two elements. The input is not modified.
func Diff(x []float64) []float64 {
	if len(x) < 2 {
		return []float64{}
	}
	out := make([]float64, len(x)-1)
	for i := 0; i < len(x)-1; i++ {
		out[i] = x[i+1] - x[i]
	}
	return out
}

// RMS returns the root-mean-square value of x, √(Σ x[i]²/N). It is a measure of
// the effective amplitude or power of the signal. For an empty input it returns
// 0.
func RMS(x []float64) float64 {
	if len(x) == 0 {
		return 0
	}
	var sum float64
	for _, v := range x {
		sum += v * v
	}
	return math.Sqrt(sum / float64(len(x)))
}

// Energy returns the total energy of x, Σ x[i]². For an empty input it returns
// 0.
func Energy(x []float64) float64 {
	var sum float64
	for _, v := range x {
		sum += v * v
	}
	return sum
}

// ZeroPad returns a copy of x extended with trailing zeros to total length n.
// If n <= len(x) the leading n samples of x are returned unchanged (a copy),
// so the result always has length max(0, n). The input is not modified.
func ZeroPad(x []float64, n int) []float64 {
	if n < 0 {
		n = 0
	}
	out := make([]float64, n)
	m := len(x)
	if m > n {
		m = n
	}
	copy(out, x[:m])
	return out
}
