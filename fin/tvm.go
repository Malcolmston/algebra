package fin

import "math"

// PaymentTiming selects whether the periodic payments of an annuity occur at
// the end of each period (an ordinary annuity) or at the beginning (an annuity
// due). It is the analogue of the "type" argument in spreadsheet TVM functions.
type PaymentTiming int

const (
	// EndOfPeriod indicates payments occur at the end of each period (an
	// ordinary annuity). This is the default in most financial software.
	EndOfPeriod PaymentTiming = 0
	// BeginningOfPeriod indicates payments occur at the start of each period
	// (an annuity due), such as a lease paid in advance.
	BeginningOfPeriod PaymentTiming = 1
)

// finTimingFactor returns 1 for BeginningOfPeriod and 0 for EndOfPeriod, the
// value that scales the per-period rate in the annuity cash-flow equation.
func finTimingFactor(when PaymentTiming) float64 {
	if when == BeginningOfPeriod {
		return 1
	}
	return 0
}

// FV returns the future value of an investment described by the standard
// cash-flow equation, given a constant per-period interest rate, a number of
// periods nper, a fixed per-period payment pmt, a present value pv, and the
// payment timing when. Following spreadsheet convention, money paid out is
// negative and money received is positive; for a loan or savings plan with
// positive deposits pmt and pv, the returned future value is negative.
func FV(rate float64, nper, pmt, pv float64, when PaymentTiming) float64 {
	t := finTimingFactor(when)
	if rate == 0 {
		return -(pv + pmt*nper)
	}
	g := math.Pow(1+rate, nper)
	return -(pv*g + pmt*(1+rate*t)*(g-1)/rate)
}

// PV returns the present value of a stream given a constant per-period rate,
// nper periods, a per-period payment pmt, a future value fv, and payment
// timing when. It inverts the same cash-flow equation as [FV]. With a positive
// future value the present value is typically negative (an outflow today).
func PV(rate float64, nper, pmt, fv float64, when PaymentTiming) float64 {
	t := finTimingFactor(when)
	if rate == 0 {
		return -(fv + pmt*nper)
	}
	g := math.Pow(1+rate, nper)
	return -(fv + pmt*(1+rate*t)*(g-1)/rate) / g
}

// PMT returns the constant per-period payment that amortises a present value
// pv to a future value fv over nper periods at the given per-period rate, with
// payment timing when. For a standard loan (positive pv, zero fv) the payment
// is negative, representing an outflow each period.
func PMT(rate float64, nper, pv, fv float64, when PaymentTiming) float64 {
	t := finTimingFactor(when)
	if nper == 0 {
		return math.NaN()
	}
	if rate == 0 {
		return -(pv + fv) / nper
	}
	g := math.Pow(1+rate, nper)
	return -(pv*g + fv) / ((1 + rate*t) * (g - 1) / rate)
}

// NPER returns the number of periods required to move from present value pv to
// future value fv with a constant per-period payment pmt at per-period rate,
// and payment timing when. It returns NaN when no positive solution exists
// (for example when the signs make the logarithm undefined).
func NPER(rate, pmt, pv, fv float64, when PaymentTiming) float64 {
	t := finTimingFactor(when)
	if rate == 0 {
		if pmt == 0 {
			return math.NaN()
		}
		return -(pv + fv) / pmt
	}
	a := pmt * (1 + rate*t)
	num := a - fv*rate
	den := a + pv*rate
	if num <= 0 || den <= 0 || den == 0 {
		return math.NaN()
	}
	return math.Log(num/den) / math.Log(1+rate)
}

// RATE returns the constant per-period interest rate that satisfies the
// cash-flow equation for nper periods, per-period payment pmt, present value
// pv, future value fv and payment timing when, starting from guess. It uses
// Newton's method with a numerical derivative and falls back to bisection over
// a wide bracket; it returns NaN if no rate can be found. Pass guess = 0.1 for
// a typical starting point.
func RATE(nper, pmt, pv, fv float64, when PaymentTiming, guess float64) float64 {
	t := finTimingFactor(when)
	f := func(r float64) float64 {
		if r == 0 {
			return pv + pmt*nper + fv
		}
		g := math.Pow(1+r, nper)
		return pv*g + pmt*(1+r*t)*(g-1)/r + fv
	}
	// Newton's method with a central-difference derivative.
	r := guess
	const h = 1e-7
	for i := 0; i < 100; i++ {
		fr := f(r)
		if math.Abs(fr) < 1e-10 {
			return r
		}
		d := (f(r+h) - f(r-h)) / (2 * h)
		if d == 0 || math.IsNaN(d) {
			break
		}
		next := r - fr/d
		if math.IsNaN(next) || math.IsInf(next, 0) {
			break
		}
		if math.Abs(next-r) < 1e-12 {
			return next
		}
		r = next
	}
	return finBisectRoot(f, -0.9999999, 1e6)
}

// FVLumpSum returns the future value of a single present amount pv compounded
// for nper periods at per-period rate: pv·(1+rate)^nper.
func FVLumpSum(pv, rate, nper float64) float64 {
	return pv * math.Pow(1+rate, nper)
}

// PVLumpSum returns the present value of a single future amount fv discounted
// over nper periods at per-period rate: fv/(1+rate)^nper.
func PVLumpSum(fv, rate, nper float64) float64 {
	return fv / math.Pow(1+rate, nper)
}

// FVAnnuity returns the future value of an ordinary annuity of nper level
// payments pmt at per-period rate: pmt·((1+rate)^nper−1)/rate. Payments are
// treated as positive inflows, so the result is positive.
func FVAnnuity(pmt, rate, nper float64) float64 {
	if rate == 0 {
		return pmt * nper
	}
	return pmt * (math.Pow(1+rate, nper) - 1) / rate
}

// PVAnnuity returns the present value of an ordinary annuity of nper level
// payments pmt at per-period rate: pmt·(1−(1+rate)^−nper)/rate.
func PVAnnuity(pmt, rate, nper float64) float64 {
	if rate == 0 {
		return pmt * nper
	}
	return pmt * (1 - math.Pow(1+rate, -nper)) / rate
}

// FVAnnuityDue returns the future value of an annuity due (payments at the
// beginning of each period) of nper level payments pmt at per-period rate. It
// equals [FVAnnuity] scaled by (1+rate).
func FVAnnuityDue(pmt, rate, nper float64) float64 {
	return FVAnnuity(pmt, rate, nper) * (1 + rate)
}

// PVAnnuityDue returns the present value of an annuity due (payments at the
// beginning of each period) of nper level payments pmt at per-period rate. It
// equals [PVAnnuity] scaled by (1+rate).
func PVAnnuityDue(pmt, rate, nper float64) float64 {
	return PVAnnuity(pmt, rate, nper) * (1 + rate)
}

// PVPerpetuity returns the present value of a perpetuity paying pmt every
// period forever at per-period rate: pmt/rate. It returns NaN for rate ≤ 0.
func PVPerpetuity(pmt, rate float64) float64 {
	if rate <= 0 {
		return math.NaN()
	}
	return pmt / rate
}

// PVGrowingPerpetuity returns the present value of a perpetuity whose first
// payment pmt grows at rate growth each period, discounted at rate: the
// Gordon growth value pmt/(rate−growth). It returns NaN when rate ≤ growth.
func PVGrowingPerpetuity(pmt, rate, growth float64) float64 {
	if rate <= growth {
		return math.NaN()
	}
	return pmt / (rate - growth)
}

// PVGrowingAnnuity returns the present value of a growing annuity: nper
// payments whose first payment is pmt and which grow at rate growth each
// period, discounted at rate. It handles the rate == growth case exactly.
func PVGrowingAnnuity(pmt, rate, growth, nper float64) float64 {
	if rate == growth {
		return pmt * nper / (1 + rate)
	}
	return pmt / (rate - growth) * (1 - math.Pow((1+growth)/(1+rate), nper))
}

// FVContinuous returns the future value of pv under continuous compounding at
// annual rate for time t years: pv·e^(rate·t).
func FVContinuous(pv, rate, t float64) float64 {
	return pv * math.Exp(rate*t)
}

// PVContinuous returns the present value of fv under continuous compounding at
// annual rate for time t years: fv·e^(−rate·t).
func PVContinuous(fv, rate, t float64) float64 {
	return fv * math.Exp(-rate*t)
}

// EffectiveRate converts a nominal annual rate compounded m times per year into
// the equivalent effective annual rate (EAR): (1+nominal/m)^m − 1. It returns
// NaN for m ≤ 0.
func EffectiveRate(nominal float64, m float64) float64 {
	if m <= 0 {
		return math.NaN()
	}
	return math.Pow(1+nominal/m, m) - 1
}

// NominalRate converts an effective annual rate into the equivalent nominal
// annual rate compounded m times per year: m·((1+effective)^(1/m) − 1). It
// returns NaN for m ≤ 0.
func NominalRate(effective float64, m float64) float64 {
	if m <= 0 {
		return math.NaN()
	}
	return m * (math.Pow(1+effective, 1/m) - 1)
}

// EffectiveFromContinuous converts a continuously compounded annual rate into
// the equivalent effective annual rate: e^rate − 1.
func EffectiveFromContinuous(rate float64) float64 {
	return math.Expm1(rate)
}

// ContinuousFromEffective converts an effective annual rate into the equivalent
// continuously compounded annual rate: ln(1+effective). It returns NaN when
// effective ≤ −1.
func ContinuousFromEffective(effective float64) float64 {
	if effective <= -1 {
		return math.NaN()
	}
	return math.Log1p(effective)
}

// DiscountFactor returns the single-period-to-nper discount factor
// 1/(1+rate)^nper, the present value of one unit of currency received nper
// periods from now.
func DiscountFactor(rate, nper float64) float64 {
	return math.Pow(1+rate, -nper)
}

// AnnuityFactor returns the present-value annuity factor
// (1−(1+rate)^−nper)/rate, the present value of a unit ordinary annuity. It is
// nper when rate == 0.
func AnnuityFactor(rate, nper float64) float64 {
	if rate == 0 {
		return nper
	}
	return (1 - math.Pow(1+rate, -nper)) / rate
}

// finBisectRoot returns a root of f in [lo, hi] by bisection, assuming f(lo)
// and f(hi) have opposite signs. It returns NaN if the interval is not
// bracketed. The iteration count is fixed so the result is deterministic.
func finBisectRoot(f func(float64) float64, lo, hi float64) float64 {
	flo, fhi := f(lo), f(hi)
	if math.IsNaN(flo) || math.IsNaN(fhi) || flo*fhi > 0 {
		return math.NaN()
	}
	for i := 0; i < 200; i++ {
		mid := (lo + hi) / 2
		fm := f(mid)
		if fm == 0 || (hi-lo)/2 < 1e-14 {
			return mid
		}
		if flo*fm < 0 {
			hi, fhi = mid, fm
		} else {
			lo, flo = mid, fm
		}
	}
	return (lo + hi) / 2
}
