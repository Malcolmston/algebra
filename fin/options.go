package fin

import "math"

// OptionType distinguishes a call option from a put option in the
// Black–Scholes–Merton pricing functions.
type OptionType int

const (
	// Call is a call option, giving the right to buy the underlying at the
	// strike price.
	Call OptionType = iota
	// Put is a put option, giving the right to sell the underlying at the
	// strike price.
	Put
)

// Greeks bundles the first-order option sensitivities (and second-order Gamma)
// returned by [ComputeGreeks]. All values are per unit change in the underlying
// parameter: Vega per 1.00 of volatility, Theta per year, Rho per 1.00 of the
// interest rate.
type Greeks struct {
	// Delta is ∂price/∂spot, the sensitivity to the underlying price.
	Delta float64
	// Gamma is ∂²price/∂spot², the rate of change of Delta.
	Gamma float64
	// Vega is ∂price/∂volatility, per 1.00 (100%) of volatility.
	Vega float64
	// Theta is ∂price/∂time, per year (typically negative).
	Theta float64
	// Rho is ∂price/∂rate, per 1.00 (100%) of the risk-free rate.
	Rho float64
}

// finNormCDF returns the standard normal cumulative distribution function Φ(x)
// using the complementary error function for accuracy in the tails.
func finNormCDF(x float64) float64 {
	return 0.5 * math.Erfc(-x/math.Sqrt2)
}

// finNormPDF returns the standard normal probability density φ(x).
func finNormPDF(x float64) float64 {
	return math.Exp(-x*x/2) / math.Sqrt(2*math.Pi)
}

// finD1D2 returns the Black–Scholes d1 and d2 terms for spot s, strike k, time
// t, risk-free rate r, dividend yield q and volatility sigma.
func finD1D2(s, k, t, r, q, sigma float64) (d1, d2 float64) {
	vsqt := sigma * math.Sqrt(t)
	d1 = (math.Log(s/k) + (r-q+sigma*sigma/2)*t) / vsqt
	d2 = d1 - vsqt
	return d1, d2
}

// BlackScholesPrice returns the Black–Scholes–Merton price of a European option
// of the given OptionType, for spot price s, strike k, time to expiry t (in
// years), risk-free rate r, dividend yield q and volatility sigma (all rates
// annualised, continuously compounded). See [BlackScholesCall] and
// [BlackScholesPut] for type-specific wrappers.
func BlackScholesPrice(otype OptionType, s, k, t, r, q, sigma float64) float64 {
	if otype == Put {
		return BlackScholesPut(s, k, t, r, q, sigma)
	}
	return BlackScholesCall(s, k, t, r, q, sigma)
}

// BlackScholesCall returns the price of a European call option under the
// Black–Scholes–Merton model with continuous dividend yield q. The parameters
// are spot s, strike k, time to expiry t in years, risk-free rate r and
// volatility sigma.
func BlackScholesCall(s, k, t, r, q, sigma float64) float64 {
	d1, d2 := finD1D2(s, k, t, r, q, sigma)
	return s*math.Exp(-q*t)*finNormCDF(d1) - k*math.Exp(-r*t)*finNormCDF(d2)
}

// BlackScholesPut returns the price of a European put option under the
// Black–Scholes–Merton model with continuous dividend yield q, for spot s,
// strike k, time to expiry t in years, risk-free rate r and volatility sigma.
func BlackScholesPut(s, k, t, r, q, sigma float64) float64 {
	d1, d2 := finD1D2(s, k, t, r, q, sigma)
	return k*math.Exp(-r*t)*finNormCDF(-d2) - s*math.Exp(-q*t)*finNormCDF(-d1)
}

// Delta returns the option delta ∂price/∂spot for an option of the given
// OptionType under Black–Scholes–Merton, with the usual parameters. A call
// delta lies in [0, e^(−qt)] and a put delta in [−e^(−qt), 0].
func Delta(otype OptionType, s, k, t, r, q, sigma float64) float64 {
	d1, _ := finD1D2(s, k, t, r, q, sigma)
	if otype == Put {
		return -math.Exp(-q*t) * finNormCDF(-d1)
	}
	return math.Exp(-q*t) * finNormCDF(d1)
}

// Gamma returns the option gamma ∂²price/∂spot², the rate of change of [Delta].
// Gamma is identical for calls and puts with the same parameters.
func Gamma(s, k, t, r, q, sigma float64) float64 {
	d1, _ := finD1D2(s, k, t, r, q, sigma)
	return math.Exp(-q*t) * finNormPDF(d1) / (s * sigma * math.Sqrt(t))
}

// Vega returns the option vega ∂price/∂volatility, the sensitivity to a 1.00
// (100%) change in volatility. Divide by 100 for the change per 1% (percentage
// point). Vega is identical for calls and puts with the same parameters.
func Vega(s, k, t, r, q, sigma float64) float64 {
	d1, _ := finD1D2(s, k, t, r, q, sigma)
	return s * math.Exp(-q*t) * finNormPDF(d1) * math.Sqrt(t)
}

// Theta returns the option theta ∂price/∂time for an option of the given
// OptionType, expressed per year (usually negative, reflecting time decay).
// Divide by 365 for the approximate one-day theta.
func Theta(otype OptionType, s, k, t, r, q, sigma float64) float64 {
	d1, d2 := finD1D2(s, k, t, r, q, sigma)
	term1 := -s * math.Exp(-q*t) * finNormPDF(d1) * sigma / (2 * math.Sqrt(t))
	if otype == Put {
		return term1 + r*k*math.Exp(-r*t)*finNormCDF(-d2) - q*s*math.Exp(-q*t)*finNormCDF(-d1)
	}
	return term1 - r*k*math.Exp(-r*t)*finNormCDF(d2) + q*s*math.Exp(-q*t)*finNormCDF(d1)
}

// Rho returns the option rho ∂price/∂rate for an option of the given
// OptionType, the sensitivity to a 1.00 (100%) change in the risk-free rate.
// Divide by 100 for the change per 1% (percentage point).
func Rho(otype OptionType, s, k, t, r, q, sigma float64) float64 {
	_, d2 := finD1D2(s, k, t, r, q, sigma)
	if otype == Put {
		return -k * t * math.Exp(-r*t) * finNormCDF(-d2)
	}
	return k * t * math.Exp(-r*t) * finNormCDF(d2)
}

// ComputeGreeks returns all option sensitivities in a single [Greeks] value for
// an option of the given OptionType, avoiding repeated evaluation of the shared
// d1/d2 terms. Field conventions match the individual Greek functions.
func ComputeGreeks(otype OptionType, s, k, t, r, q, sigma float64) Greeks {
	return Greeks{
		Delta: Delta(otype, s, k, t, r, q, sigma),
		Gamma: Gamma(s, k, t, r, q, sigma),
		Vega:  Vega(s, k, t, r, q, sigma),
		Theta: Theta(otype, s, k, t, r, q, sigma),
		Rho:   Rho(otype, s, k, t, r, q, sigma),
	}
}

// ImpliedVolatility returns the volatility that makes the Black–Scholes–Merton
// price of an option of the given OptionType equal the observed market price,
// given spot s, strike k, time to expiry t, rate r and dividend yield q. It
// uses Newton's method seeded at 0.2 with a bisection fallback over (0, 5], and
// returns NaN when the price is not attainable (for example below intrinsic
// value).
func ImpliedVolatility(otype OptionType, price, s, k, t, r, q float64) float64 {
	f := func(sig float64) float64 {
		return BlackScholesPrice(otype, s, k, t, r, q, sig) - price
	}
	sig := 0.2
	for i := 0; i < 100; i++ {
		fv := f(sig)
		if math.Abs(fv) < 1e-8 {
			return sig
		}
		v := Vega(s, k, t, r, q, sig)
		if v < 1e-12 {
			break
		}
		next := sig - fv/v
		if next <= 0 || math.IsNaN(next) || math.IsInf(next, 0) {
			break
		}
		if math.Abs(next-sig) < 1e-10 {
			return next
		}
		sig = next
	}
	return finBisectRoot(f, 1e-6, 5.0)
}
