package fin

import "math"

// BondPrice returns the fair price of a fixed-coupon bond given its face value,
// annual couponRate (as a decimal, e.g. 0.05 for 5%), annual yield to maturity
// yield, time to maturity in years, and the number of coupon payments per year
// freq (e.g. 2 for semi-annual). The price is the present value of all future
// coupons plus the redemption of face at maturity, discounted at yield/freq per
// period.
func BondPrice(face, couponRate, yield, years float64, freq int) float64 {
	if freq <= 0 {
		return math.NaN()
	}
	m := float64(freq)
	n := int(math.Round(years * m))
	i := yield / m
	coupon := face * couponRate / m
	var price float64
	for k := 1; k <= n; k++ {
		price += coupon / math.Pow(1+i, float64(k))
	}
	price += face / math.Pow(1+i, float64(n))
	return price
}

// ZeroCouponPrice returns the price of a zero-coupon bond: the present value of
// its face value discounted at annual yield compounded freq times per year over
// the given years to maturity. With freq = 1 this is face/(1+yield)^years.
func ZeroCouponPrice(face, yield, years float64, freq int) float64 {
	if freq <= 0 {
		return math.NaN()
	}
	m := float64(freq)
	return face / math.Pow(1+yield/m, years*m)
}

// BondYTM returns the annual yield to maturity that makes the [BondPrice] of a
// bond with the given face value, annual couponRate, years to maturity and
// coupon frequency freq equal to its observed market price. It solves by
// bisection over (0, 1] (0% to 100%) and returns NaN if no yield in that range
// reproduces the price.
func BondYTM(face, couponRate, price, years float64, freq int) float64 {
	if freq <= 0 {
		return math.NaN()
	}
	f := func(y float64) float64 { return BondPrice(face, couponRate, y, years, freq) - price }
	return finBisectRoot(f, -0.9999, 1.0)
}

// CurrentYield returns a bond's current yield: its annual coupon income
// divided by its market price, face·couponRate/price. It returns NaN for a
// non-positive price.
func CurrentYield(face, couponRate, price float64) float64 {
	if price <= 0 {
		return math.NaN()
	}
	return face * couponRate / price
}

// AccruedInterest returns the coupon interest accrued between coupon dates on a
// bond with the given face value and annual couponRate, where daysSinceCoupon
// days have elapsed of a coupon period of daysInPeriod days. It uses simple
// straight-line (actual/actual within period) accrual.
func AccruedInterest(face, couponRate float64, freq int, daysSinceCoupon, daysInPeriod float64) float64 {
	if freq <= 0 || daysInPeriod <= 0 {
		return math.NaN()
	}
	periodCoupon := face * couponRate / float64(freq)
	return periodCoupon * daysSinceCoupon / daysInPeriod
}

// MacaulayDuration returns the Macaulay duration (in years) of a fixed-coupon
// bond: the present-value-weighted average time to each cash flow, discounted
// at yield. Duration measures a bond's effective maturity and underlies its
// interest-rate sensitivity.
func MacaulayDuration(face, couponRate, yield, years float64, freq int) float64 {
	if freq <= 0 {
		return math.NaN()
	}
	m := float64(freq)
	n := int(math.Round(years * m))
	i := yield / m
	coupon := face * couponRate / m
	var price, weighted float64
	for k := 1; k <= n; k++ {
		cf := coupon
		if k == n {
			cf += face
		}
		pv := cf / math.Pow(1+i, float64(k))
		price += pv
		weighted += (float64(k) / m) * pv
	}
	if price == 0 {
		return math.NaN()
	}
	return weighted / price
}

// ModifiedDuration returns the modified duration (in years) of a fixed-coupon
// bond, the Macaulay duration divided by (1+yield/freq). It approximates the
// percentage change in price for a one-unit (100%) change in yield; multiply by
// a small yield change to estimate the price move.
func ModifiedDuration(face, couponRate, yield, years float64, freq int) float64 {
	m := float64(freq)
	mac := MacaulayDuration(face, couponRate, yield, years, freq)
	return mac / (1 + yield/m)
}

// Convexity returns the convexity (in years squared) of a fixed-coupon bond,
// the second-order sensitivity of price to yield. Together with
// [ModifiedDuration] it gives the price approximation
// ΔP/P ≈ −ModDur·Δy + ½·Convexity·Δy², improving the linear duration estimate
// for larger yield moves.
func Convexity(face, couponRate, yield, years float64, freq int) float64 {
	if freq <= 0 {
		return math.NaN()
	}
	m := float64(freq)
	n := int(math.Round(years * m))
	i := yield / m
	coupon := face * couponRate / m
	var price, c float64
	for k := 1; k <= n; k++ {
		cf := coupon
		if k == n {
			cf += face
		}
		pv := cf / math.Pow(1+i, float64(k))
		price += pv
		c += float64(k) * float64(k+1) * pv
	}
	if price == 0 {
		return math.NaN()
	}
	// Second derivative of price w.r.t. per-period yield, annualised by m².
	return c / (price * math.Pow(1+i, 2) * m * m)
}
