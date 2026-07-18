package plot

import (
	"math"
	"strconv"
)

// niceNum returns a "nice" number approximately equal to x. When round is true
// the result is rounded to the nearest nice value; otherwise it is the
// smallest nice value not less than x. Nice numbers are 1, 2, 5 and 10 scaled
// by a power of ten, which is the classic Wilkinson tick heuristic.
func niceNum(x float64, round bool) float64 {
	if x <= 0 || math.IsNaN(x) || math.IsInf(x, 0) {
		return 0
	}
	exp := math.Floor(math.Log10(x))
	frac := x / math.Pow(10, exp)
	var nf float64
	if round {
		switch {
		case frac < 1.5:
			nf = 1
		case frac < 3:
			nf = 2
		case frac < 7:
			nf = 5
		default:
			nf = 10
		}
	} else {
		switch {
		case frac <= 1:
			nf = 1
		case frac <= 2:
			nf = 2
		case frac <= 5:
			nf = 5
		default:
			nf = 10
		}
	}
	return nf * math.Pow(10, exp)
}

// ticksFor computes a set of evenly spaced "nice" tick locations covering the
// interval [lo, hi] using at most about maxTicks intervals. It returns the tick
// values together with the tick spacing. The returned ticks all lie within
// [lo, hi] (inclusive, within a small epsilon).
func ticksFor(lo, hi float64, maxTicks int) ([]float64, float64) {
	if maxTicks < 1 {
		maxTicks = 1
	}
	if hi < lo {
		lo, hi = hi, lo
	}
	span := hi - lo
	if span == 0 || math.IsNaN(span) || math.IsInf(span, 0) {
		// Degenerate range: emit a single tick at lo.
		return []float64{lo}, 0
	}
	rng := niceNum(span, false)
	step := niceNum(rng/float64(maxTicks), true)
	if step <= 0 {
		return []float64{lo}, 0
	}
	start := math.Ceil(lo/step) * step
	end := math.Floor(hi/step) * step
	var out []float64
	eps := step * 1e-9
	for v := start; v <= end+eps; v += step {
		// Snap values extremely close to zero to exactly zero for clean labels.
		if math.Abs(v) < eps {
			v = 0
		}
		out = append(out, v)
	}
	return out, step
}

// formatTick formats a tick value for display, choosing a fixed or scientific
// representation based on the tick spacing so labels stay compact and free of
// floating point noise.
func formatTick(v, step float64) string {
	if v == 0 {
		return "0"
	}
	abs := math.Abs(v)
	if abs >= 1e5 || abs < 1e-4 {
		return strconv.FormatFloat(v, 'g', 4, 64)
	}
	// Choose decimal places from the step magnitude.
	dec := 0
	if step > 0 && step < 1 {
		dec = int(math.Ceil(-math.Log10(step)))
		if dec > 6 {
			dec = 6
		}
	}
	s := strconv.FormatFloat(v, 'f', dec, 64)
	return s
}

// padRange expands [lo, hi] outward by frac on each side so that data does not
// touch the plot border. A degenerate range (lo == hi) is expanded to a unit
// window centred on the value.
func padRange(lo, hi, frac float64) (float64, float64) {
	if math.IsNaN(lo) || math.IsNaN(hi) {
		return 0, 1
	}
	if hi < lo {
		lo, hi = hi, lo
	}
	span := hi - lo
	if span == 0 {
		if lo == 0 {
			return -1, 1
		}
		d := math.Abs(lo) * 0.5
		return lo - d, hi + d
	}
	return lo - span*frac, hi + span*frac
}
