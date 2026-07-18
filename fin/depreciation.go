package fin

import "math"

// StraightLineDepreciation returns the constant per-period depreciation expense
// of an asset under the straight-line method: (cost − salvage)/life. It returns
// NaN for life ≤ 0.
func StraightLineDepreciation(cost, salvage, life float64) float64 {
	if life <= 0 {
		return math.NaN()
	}
	return (cost - salvage) / life
}

// StraightLineSchedule returns the per-period depreciation expense for each of
// the life periods under the straight-line method. Every element equals
// [StraightLineDepreciation]; the slice has one entry per period.
func StraightLineSchedule(cost, salvage float64, life int) []float64 {
	if life <= 0 {
		return nil
	}
	d := (cost - salvage) / float64(life)
	out := make([]float64, life)
	for i := range out {
		out[i] = d
	}
	return out
}

// DecliningBalance returns the depreciation expense in a single 1-based period
// under the declining-balance method at the given per-period rate (as a
// decimal). Depreciation is charged on the reducing book value but is capped so
// the book value never falls below salvage.
func DecliningBalance(cost, salvage, rate float64, period, life int) float64 {
	if period < 1 || period > life {
		return 0
	}
	book := cost
	var d float64
	for p := 1; p <= period; p++ {
		d = book * rate
		if book-d < salvage {
			d = book - salvage
		}
		if d < 0 {
			d = 0
		}
		book -= d
	}
	return d
}

// DoubleDecliningBalance returns the depreciation expense in a single 1-based
// period under the double-declining-balance method, which uses a rate of
// 2/life on the reducing book value, capped so book value never falls below
// salvage. It returns NaN for life ≤ 0.
func DoubleDecliningBalance(cost, salvage float64, period, life int) float64 {
	if life <= 0 {
		return math.NaN()
	}
	return DecliningBalance(cost, salvage, 2/float64(life), period, life)
}

// DoubleDecliningSchedule returns the depreciation expense for each of the life
// periods under the double-declining-balance method. The reported expenses sum
// to at most cost − salvage.
func DoubleDecliningSchedule(cost, salvage float64, life int) []float64 {
	if life <= 0 {
		return nil
	}
	rate := 2 / float64(life)
	out := make([]float64, life)
	book := cost
	for p := 0; p < life; p++ {
		d := book * rate
		if book-d < salvage {
			d = book - salvage
		}
		if d < 0 {
			d = 0
		}
		out[p] = d
		book -= d
	}
	return out
}

// SumOfYearsDigits returns the depreciation expense in a single 1-based period
// under the sum-of-years'-digits method. The remaining life (life−period+1) is
// weighted by the total of the digits 1..life, front-loading depreciation. It
// returns NaN for life ≤ 0.
func SumOfYearsDigits(cost, salvage float64, period, life int) float64 {
	if life <= 0 {
		return math.NaN()
	}
	if period < 1 || period > life {
		return 0
	}
	syd := float64(life) * float64(life+1) / 2
	return (cost - salvage) * float64(life-period+1) / syd
}

// SumOfYearsDigitsSchedule returns the depreciation expense for each of the
// life periods under the sum-of-years'-digits method. The expenses sum exactly
// to cost − salvage.
func SumOfYearsDigitsSchedule(cost, salvage float64, life int) []float64 {
	if life <= 0 {
		return nil
	}
	syd := float64(life) * float64(life+1) / 2
	out := make([]float64, life)
	for p := 1; p <= life; p++ {
		out[p-1] = (cost - salvage) * float64(life-p+1) / syd
	}
	return out
}
