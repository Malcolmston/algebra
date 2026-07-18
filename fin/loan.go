package fin

import "math"

// AmortizationRow describes one period of a fully-amortising loan repayment
// schedule.
type AmortizationRow struct {
	// Period is the 1-based payment number.
	Period int
	// Payment is the total payment made in this period (interest + principal).
	Payment float64
	// Interest is the portion of the payment covering interest.
	Interest float64
	// Principal is the portion of the payment reducing the balance.
	Principal float64
	// Balance is the remaining principal after this payment.
	Balance float64
}

// AmortizationSchedule returns the period-by-period repayment schedule of a
// fully-amortising loan of the given principal at per-period rate over nper
// periods. Each row splits the level payment into its interest and principal
// components and reports the remaining balance. Payment and its components are
// reported as positive amounts. The final balance is driven to zero exactly.
func AmortizationSchedule(principal, rate float64, nper int) []AmortizationRow {
	if nper <= 0 {
		return nil
	}
	pmt := -PMT(rate, float64(nper), principal, 0, EndOfPeriod)
	rows := make([]AmortizationRow, 0, nper)
	balance := principal
	for p := 1; p <= nper; p++ {
		interest := balance * rate
		prin := pmt - interest
		balance -= prin
		if p == nper {
			// Absorb any residual rounding into the last payment.
			prin += balance
			pmt += balance
			balance = 0
		}
		rows = append(rows, AmortizationRow{
			Period:    p,
			Payment:   pmt,
			Interest:  interest,
			Principal: prin,
			Balance:   balance,
		})
	}
	return rows
}

// RemainingBalance returns the outstanding principal on a fully-amortising loan
// of the given original principal at per-period rate over nper periods, after
// periodsPaid level payments have been made. It equals the present value of the
// remaining payments.
func RemainingBalance(principal, rate float64, nper, periodsPaid int) float64 {
	if periodsPaid <= 0 {
		return principal
	}
	if periodsPaid >= nper {
		return 0
	}
	pmt := -PMT(rate, float64(nper), principal, 0, EndOfPeriod)
	if rate == 0 {
		return principal - pmt*float64(periodsPaid)
	}
	g := math.Pow(1+rate, float64(periodsPaid))
	return principal*g - pmt*(g-1)/rate
}

// TotalInterest returns the total interest paid over the life of a
// fully-amortising loan of the given principal at per-period rate over nper
// periods: the sum of all level payments minus the principal.
func TotalInterest(principal, rate float64, nper int) float64 {
	if nper <= 0 {
		return math.NaN()
	}
	pmt := -PMT(rate, float64(nper), principal, 0, EndOfPeriod)
	return pmt*float64(nper) - principal
}
