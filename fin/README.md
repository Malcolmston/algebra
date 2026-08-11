# fin — small, dependency-free financial mathematics toolkit written entirely with the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/fin.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/fin)

Package fin is a small, dependency-free financial mathematics toolkit written
entirely with the Go standard library (only the math package). It is
completely standalone: it does not import the parent algebra package, and it
performs no I/O, no concurrency and no random sampling, so every result is
deterministic.

The package is organised into five groups.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/fin"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AccruedInterest(face, couponRate float64, freq int, daysSinceCoupon, daysInPeriod float64) float64` | AccruedInterest returns the coupon interest accrued between coupon dates on a bond with the given face value and annual couponRate, where… |
| `func Annualize(periodVol float64, periodsPerYear float64) float64` | Annualize scales a per-period volatility to an annual figure given the number of periods per year, multiplying by √periodsPerYear (the… |
| `func AnnualizedReturn(returns []float64) float64` | AnnualizedReturn returns the geometric mean per-period return of a series compounded to a single period figure: (∏(1+rᵢ))^(1/n) − 1. |
| `func AnnuityFactor(rate, nper float64) float64` | AnnuityFactor returns the present-value annuity factor (1−(1+rate)^−nper)/rate, the present value of a unit ordinary annuity. |
| `func BlackScholesCall(s, k, t, r, q, sigma float64) float64` | BlackScholesCall returns the price of a European call option under the Black–Scholes–Merton model with continuous dividend yield q. |
| `func BlackScholesPrice(otype OptionType, s, k, t, r, q, sigma float64) float64` | BlackScholesPrice returns the Black–Scholes–Merton price of a European option of the given OptionType, for spot price s, strike k, time to expiry… |
| `func BlackScholesPut(s, k, t, r, q, sigma float64) float64` | BlackScholesPut returns the price of a European put option under the Black–Scholes–Merton model with continuous dividend yield q, for spot s,… |
| `func BondPrice(face, couponRate, yield, years float64, freq int) float64` | BondPrice returns the fair price of a fixed-coupon bond given its face value, annual couponRate (as a decimal, e.g. |
| `func BondYTM(face, couponRate, price, years float64, freq int) float64` | BondYTM returns the annual yield to maturity that makes the `BondPrice` of a bond with the given face value, annual couponRate, years to maturity and… |
| `func CAGR(begin, end, years float64) float64` | CAGR returns the compound annual growth rate that takes an initial value begin to a final value end over the given number of years:… |
| `func ContinuousFromEffective(effective float64) float64` | ContinuousFromEffective converts an effective annual rate into the equivalent continuously compounded annual rate: ln(1+effective). |
| `func Convexity(face, couponRate, yield, years float64, freq int) float64` | Convexity returns the convexity (in years squared) of a fixed-coupon bond, the second-order sensitivity of price to yield. |
| `func CurrentYield(face, couponRate, price float64) float64` | CurrentYield returns a bond's current yield: its annual coupon income divided by its market price, face·couponRate/price. |
| `func DecliningBalance(cost, salvage, rate float64, period, life int) float64` | DecliningBalance returns the depreciation expense in a single 1-based period under the declining-balance method at the given per-period rate (as a… |
| `func Delta(otype OptionType, s, k, t, r, q, sigma float64) float64` | Delta returns the option delta ∂price/∂spot for an option of the given OptionType under Black–Scholes–Merton, with the usual parameters. |
| `func DiscountFactor(rate, nper float64) float64` | DiscountFactor returns the single-period-to-nper discount factor 1/(1+rate)^nper, the present value of one unit of currency received nper periods… |
| `func DoubleDecliningBalance(cost, salvage float64, period, life int) float64` | DoubleDecliningBalance returns the depreciation expense in a single 1-based period under the double-declining-balance method, which uses a rate of… |
| `func DoubleDecliningSchedule(cost, salvage float64, life int) []float64` | DoubleDecliningSchedule returns the depreciation expense for each of the life periods under the double-declining-balance method. |
| `func EffectiveFromContinuous(rate float64) float64` | EffectiveFromContinuous converts a continuously compounded annual rate into the equivalent effective annual rate: e^rate − 1. |
| `func EffectiveRate(nominal float64, m float64) float64` | EffectiveRate converts a nominal annual rate compounded m times per year into the equivalent effective annual rate (EAR): (1+nominal/m)^m − 1. |
| `func FV(rate float64, nper, pmt, pv float64, when PaymentTiming) float64` | FV returns the future value of an investment described by the standard cash-flow equation, given a constant per-period interest rate, a number of… |
| `func FVAnnuity(pmt, rate, nper float64) float64` | FVAnnuity returns the future value of an ordinary annuity of nper level payments pmt at per-period rate: pmt·((1+rate)^nper−1)/rate. |
| `func FVAnnuityDue(pmt, rate, nper float64) float64` | FVAnnuityDue returns the future value of an annuity due (payments at the beginning of each period) of nper level payments pmt at per-period rate. |
| `func FVContinuous(pv, rate, t float64) float64` | FVContinuous returns the future value of pv under continuous compounding at annual rate for time t years: pv·e^(rate·t). |
| `func FVLumpSum(pv, rate, nper float64) float64` | FVLumpSum returns the future value of a single present amount pv compounded for nper periods at per-period rate: pv·(1+rate)^nper. |
| `func Gamma(s, k, t, r, q, sigma float64) float64` | Gamma returns the option gamma ∂²price/∂spot², the rate of change of `Delta`. |
| `func HoldingPeriodReturn(returns []float64) float64` | HoldingPeriodReturn returns the total return of chaining a series of periodic returns (each expressed as a decimal, e.g. |
| `func IRR(cashflows []float64, guess float64) float64` | IRR returns the internal rate of return of a series of equally-spaced cash flows: the per-period rate at which `NPV` equals zero. |
| `func IRRBisection(cashflows []float64) float64` | IRRBisection returns the internal rate of return found purely by bisection over the bracket (−0.9999999, 1e6). |
| `func ImpliedVolatility(otype OptionType, price, s, k, t, r, q float64) float64` | ImpliedVolatility returns the volatility that makes the Black–Scholes–Merton price of an option of the given OptionType equal the observed market… |
| `func LogReturn(begin, end float64) float64` | LogReturn returns the continuously compounded (log) return from an opening value begin to a closing value end: ln(end/begin). |
| `func MIRR(cashflows []float64, financeRate, reinvestRate float64) float64` | MIRR returns the modified internal rate of return of a series of equally-spaced cash flows, where negative flows are financed at financeRate and… |
| `func MacaulayDuration(face, couponRate, yield, years float64, freq int) float64` | MacaulayDuration returns the Macaulay duration (in years) of a fixed-coupon bond: the present-value-weighted average time to each cash flow,… |
| `func ModifiedDuration(face, couponRate, yield, years float64, freq int) float64` | ModifiedDuration returns the modified duration (in years) of a fixed-coupon bond, the Macaulay duration divided by (1+yield/freq). |
| `func NPER(rate, pmt, pv, fv float64, when PaymentTiming) float64` | NPER returns the number of periods required to move from present value pv to future value fv with a constant per-period payment pmt at per-period… |
| `func NPV(rate float64, cashflows []float64) float64` | NPV returns the net present value of a series of equally-spaced cash flows discounted at a constant per-period rate. |
| `func NominalRate(effective float64, m float64) float64` | NominalRate converts an effective annual rate into the equivalent nominal annual rate compounded m times per year: m·((1+effective)^(1/m) − 1). |
| `func PMT(rate float64, nper, pv, fv float64, when PaymentTiming) float64` | PMT returns the constant per-period payment that amortises a present value pv to a future value fv over nper periods at the given per-period rate,… |
| `func PV(rate float64, nper, pmt, fv float64, when PaymentTiming) float64` | PV returns the present value of a stream given a constant per-period rate, nper periods, a per-period payment pmt, a future value fv, and payment… |
| `func PVAnnuity(pmt, rate, nper float64) float64` | PVAnnuity returns the present value of an ordinary annuity of nper level payments pmt at per-period rate: pmt·(1−(1+rate)^−nper)/rate. |
| `func PVAnnuityDue(pmt, rate, nper float64) float64` | PVAnnuityDue returns the present value of an annuity due (payments at the beginning of each period) of nper level payments pmt at per-period rate. |
| `func PVContinuous(fv, rate, t float64) float64` | PVContinuous returns the present value of fv under continuous compounding at annual rate for time t years: fv·e^(−rate·t). |
| `func PVGrowingAnnuity(pmt, rate, growth, nper float64) float64` | PVGrowingAnnuity returns the present value of a growing annuity: nper payments whose first payment is pmt and which grow at rate growth each period,… |
| `func PVGrowingPerpetuity(pmt, rate, growth float64) float64` | PVGrowingPerpetuity returns the present value of a perpetuity whose first payment pmt grows at rate growth each period, discounted at rate: the… |
| `func PVLumpSum(fv, rate, nper float64) float64` | PVLumpSum returns the present value of a single future amount fv discounted over nper periods at per-period rate: fv/(1+rate)^nper. |
| `func PVPerpetuity(pmt, rate float64) float64` | PVPerpetuity returns the present value of a perpetuity paying pmt every period forever at per-period rate: pmt/rate. |
| `func RATE(nper, pmt, pv, fv float64, when PaymentTiming, guess float64) float64` | RATE returns the constant per-period interest rate that satisfies the cash-flow equation for nper periods, per-period payment pmt, present value pv,… |
| `func RemainingBalance(principal, rate float64, nper, periodsPaid int) float64` | RemainingBalance returns the outstanding principal on a fully-amortising loan of the given original principal at per-period rate over nper periods,… |
| `func Rho(otype OptionType, s, k, t, r, q, sigma float64) float64` | Rho returns the option rho ∂price/∂rate for an option of the given OptionType, the sensitivity to a 1.00 (100%) change in the risk-free rate. |
| `func SharpeRatio(returns []float64, riskFree float64) float64` | SharpeRatio returns the Sharpe ratio of a series of periodic returns given a per-period riskFree rate: the mean excess return divided by the standard… |
| `func SimpleReturn(begin, end float64) float64` | SimpleReturn returns the simple (arithmetic) return from an opening value begin to a closing value end: (end − begin)/begin. |
| `func SortinoRatio(returns []float64, riskFree float64) float64` | SortinoRatio returns the Sortino ratio of a series of periodic returns given a per-period riskFree rate: the mean excess return divided by the… |
| `func StraightLineDepreciation(cost, salvage, life float64) float64` | StraightLineDepreciation returns the constant per-period depreciation expense of an asset under the straight-line method: (cost − salvage)/life. |
| `func StraightLineSchedule(cost, salvage float64, life int) []float64` | StraightLineSchedule returns the per-period depreciation expense for each of the life periods under the straight-line method. |
| `func SumOfYearsDigits(cost, salvage float64, period, life int) float64` | SumOfYearsDigits returns the depreciation expense in a single 1-based period under the sum-of-years'-digits method. |
| `func SumOfYearsDigitsSchedule(cost, salvage float64, life int) []float64` | SumOfYearsDigitsSchedule returns the depreciation expense for each of the life periods under the sum-of-years'-digits method. |
| `func Theta(otype OptionType, s, k, t, r, q, sigma float64) float64` | Theta returns the option theta ∂price/∂time for an option of the given OptionType, expressed per year (usually negative, reflecting time decay). |
| `func TotalInterest(principal, rate float64, nper int) float64` | TotalInterest returns the total interest paid over the life of a fully-amortising loan of the given principal at per-period rate over nper periods:… |
| `func Vega(s, k, t, r, q, sigma float64) float64` | Vega returns the option vega ∂price/∂volatility, the sensitivity to a 1.00 (100%) change in volatility. |
| `func Volatility(returns []float64) float64` | Volatility returns the sample standard deviation of a series of periodic returns, the standard measure of investment risk. |
| `func XIRR(cashflows, times []float64, guess float64) float64` | XIRR returns the annual internal rate of return for cash flows occurring at arbitrary times (in years), the rate at which `XNPV` is zero. |
| `func XNPV(rate float64, cashflows, times []float64) float64` | XNPV returns the net present value of cash flows that occur at arbitrary times, discounted at an annual rate. |
| `func ZeroCouponPrice(face, yield, years float64, freq int) float64` | ZeroCouponPrice returns the price of a zero-coupon bond: the present value of its face value discounted at annual yield compounded freq times per… |

### Types

| Type | What it is |
| --- | --- |
| `AmortizationRow` | AmortizationRow describes one period of a fully-amortising loan repayment schedule. |
| `Greeks` | Greeks bundles the first-order option sensitivities (and second-order Gamma) returned by `ComputeGreeks`. |
| `OptionType` | OptionType distinguishes a call option from a put option in the Black–Scholes–Merton pricing functions. |
| `PaymentTiming` | PaymentTiming selects whether the periodic payments of an annuity occur at the end of each period (an ordinary annuity) or at the beginning (an… |

<details>
<summary><code>AmortizationRow</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AmortizationSchedule(principal, rate float64, nper int) []AmortizationRow` | AmortizationSchedule returns the period-by-period repayment schedule of a fully-amortising loan of the given principal at per-period rate over nper… |

</details>

<details>
<summary><code>Greeks</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ComputeGreeks(otype OptionType, s, k, t, r, q, sigma float64) Greeks` | ComputeGreeks returns all option sensitivities in a single `Greeks` value for an option of the given OptionType, avoiding repeated evaluation of the… |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/fin).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
