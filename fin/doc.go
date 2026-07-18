// Package fin is a small, dependency-free financial mathematics toolkit
// written entirely with the Go standard library (only the math package). It is
// completely standalone: it does not import the parent algebra package, and it
// performs no I/O, no concurrency and no random sampling, so every result is
// deterministic.
//
// The package is organised into five groups.
//
// # Time value of money
//
// The Excel-compatible core solvers [FV], [PV], [PMT], [NPER] and [RATE] all
// operate on the standard cash-flow equation
//
//	pv·(1+r)ⁿ + pmt·(1+r·t)·((1+r)ⁿ−1)/r + fv = 0
//
// where t selects [EndOfPeriod] (an ordinary annuity) or [BeginningOfPeriod]
// (an annuity due). Convenience helpers cover single sums ([FVLumpSum],
// [PVLumpSum]), level annuities ([FVAnnuity], [PVAnnuity], [FVAnnuityDue],
// [PVAnnuityDue]), perpetuities ([PVPerpetuity], [PVGrowingPerpetuity]),
// growing annuities ([PVGrowingAnnuity]), and continuous compounding
// ([FVContinuous], [PVContinuous]). Rate conversions are provided by
// [EffectiveRate], [NominalRate], [EffectiveFromContinuous] and
// [ContinuousFromEffective].
//
// # Cash-flow analysis
//
// [NPV], [IRR], [IRRBisection], [MIRR], [XNPV] and [XIRR] evaluate and solve
// discounted-cash-flow problems. [CAGR] reports a compound annual growth rate.
//
// # Fixed income
//
// [BondPrice], [BondYTM], [ZeroCouponPrice], [CurrentYield], [AccruedInterest],
// [MacaulayDuration], [ModifiedDuration] and [Convexity] price coupon bonds and
// measure their interest-rate sensitivity.
//
// # Options
//
// [BlackScholesPrice] (with [BlackScholesCall], [BlackScholesPut]) prices
// European options under the Black–Scholes–Merton model with a continuous
// dividend yield. The Greeks [Delta], [Gamma], [Vega], [Theta] and [Rho] are
// available individually or bundled by [ComputeGreeks]; [ImpliedVolatility]
// inverts the price.
//
// # Loans, depreciation and performance
//
// [AmortizationSchedule], [RemainingBalance] and [TotalInterest] analyse
// fully-amortising loans. [StraightLineDepreciation], [DecliningBalance],
// [DoubleDecliningBalance] and [SumOfYearsDigits] (with their *Schedule
// variants) model asset depreciation. [SharpeRatio], [SortinoRatio],
// [Volatility], [AnnualizedReturn] and [HoldingPeriodReturn] summarise an
// investment return series.
//
// # Conventions
//
// Following spreadsheet convention, cash paid out is negative and cash received
// is positive, so [PV] and [PMT] of a positive future value are typically
// negative. Functions that have no meaningful result for their inputs (an empty
// series, a non-converging root search, a non-positive argument to a logarithm)
// return NaN rather than panicking, so callers can test with [math.IsNaN].
package fin
