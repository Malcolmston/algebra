package fin

import "math"

// NPV returns the net present value of a series of equally-spaced cash flows
// discounted at a constant per-period rate. The first element cashflows[0] is
// treated as occurring at time 0 (it is not discounted), cashflows[1] at time
// 1, and so on. This matches the textbook definition; note it differs from the
// spreadsheet NPV function, which discounts the first element by one period.
func NPV(rate float64, cashflows []float64) float64 {
	npv := 0.0
	for i, cf := range cashflows {
		npv += cf / math.Pow(1+rate, float64(i))
	}
	return npv
}

// IRR returns the internal rate of return of a series of equally-spaced cash
// flows: the per-period rate at which [NPV] equals zero. The first cash flow is
// at time 0. It uses Newton's method from guess and falls back to bisection
// over (−1, ∞); it returns NaN if no rate is found. A conventional series has
// exactly one sign change and a unique IRR. Pass guess = 0.1 for a typical
// starting point.
func IRR(cashflows []float64, guess float64) float64 {
	f := func(r float64) float64 { return NPV(r, cashflows) }
	df := func(r float64) float64 {
		d := 0.0
		for i, cf := range cashflows {
			d += -float64(i) * cf / math.Pow(1+r, float64(i)+1)
		}
		return d
	}
	r := guess
	for i := 0; i < 100; i++ {
		fr := f(r)
		if math.Abs(fr) < 1e-10 {
			return r
		}
		d := df(r)
		if d == 0 || math.IsNaN(d) {
			break
		}
		next := r - fr/d
		if math.IsNaN(next) || math.IsInf(next, 0) || next <= -1 {
			break
		}
		if math.Abs(next-r) < 1e-12 {
			return next
		}
		r = next
	}
	return IRRBisection(cashflows)
}

// IRRBisection returns the internal rate of return found purely by bisection
// over the bracket (−0.9999999, 1e6). It is slower but more robust than the
// Newton path in [IRR] and needs no initial guess. It returns NaN when the net
// present value does not change sign across the bracket.
func IRRBisection(cashflows []float64) float64 {
	f := func(r float64) float64 { return NPV(r, cashflows) }
	return finBisectRoot(f, -0.9999999, 1e6)
}

// MIRR returns the modified internal rate of return of a series of
// equally-spaced cash flows, where negative flows are financed at financeRate
// and positive flows are reinvested at reinvestRate. Unlike [IRR] it always has
// a unique real value for a conventional project. It returns NaN when there are
// fewer than two periods or when there are no negative flows to finance.
func MIRR(cashflows []float64, financeRate, reinvestRate float64) float64 {
	n := len(cashflows) - 1
	if n < 1 {
		return math.NaN()
	}
	var pvNeg, fvPos float64
	for i, cf := range cashflows {
		if cf < 0 {
			pvNeg += cf / math.Pow(1+financeRate, float64(i))
		} else if cf > 0 {
			fvPos += cf * math.Pow(1+reinvestRate, float64(n-i))
		}
	}
	if pvNeg == 0 || fvPos == 0 {
		return math.NaN()
	}
	return math.Pow(fvPos/-pvNeg, 1/float64(n)) - 1
}

// XNPV returns the net present value of cash flows that occur at arbitrary
// times, discounted at an annual rate. The times slice gives the time in years
// of each corresponding cash flow (times[0] is usually 0). The two slices must
// have equal length; otherwise it returns NaN.
func XNPV(rate float64, cashflows, times []float64) float64 {
	if len(cashflows) != len(times) {
		return math.NaN()
	}
	npv := 0.0
	for i := range cashflows {
		npv += cashflows[i] / math.Pow(1+rate, times[i])
	}
	return npv
}

// XIRR returns the annual internal rate of return for cash flows occurring at
// arbitrary times (in years), the rate at which [XNPV] is zero. It uses
// Newton's method from guess with a bisection fallback and returns NaN if no
// rate is found or the slice lengths differ.
func XIRR(cashflows, times []float64, guess float64) float64 {
	if len(cashflows) != len(times) || len(cashflows) < 2 {
		return math.NaN()
	}
	f := func(r float64) float64 { return XNPV(r, cashflows, times) }
	df := func(r float64) float64 {
		d := 0.0
		for i := range cashflows {
			d += -times[i] * cashflows[i] / math.Pow(1+r, times[i]+1)
		}
		return d
	}
	r := guess
	for i := 0; i < 100; i++ {
		fr := f(r)
		if math.Abs(fr) < 1e-9 {
			return r
		}
		d := df(r)
		if d == 0 || math.IsNaN(d) {
			break
		}
		next := r - fr/d
		if math.IsNaN(next) || math.IsInf(next, 0) || next <= -1 {
			break
		}
		if math.Abs(next-r) < 1e-12 {
			return next
		}
		r = next
	}
	return finBisectRoot(f, -0.9999999, 1e6)
}

// CAGR returns the compound annual growth rate that takes an initial value
// begin to a final value end over the given number of years:
// (end/begin)^(1/years) − 1. It returns NaN when begin ≤ 0 or years ≤ 0.
func CAGR(begin, end, years float64) float64 {
	if begin <= 0 || years <= 0 {
		return math.NaN()
	}
	return math.Pow(end/begin, 1/years) - 1
}
