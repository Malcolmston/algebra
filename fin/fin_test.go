package fin

import (
	"math"
	"testing"
)

const tol = 1e-6

func approx(t *testing.T, got, want, eps float64, name string) {
	t.Helper()
	if math.IsNaN(want) {
		if !math.IsNaN(got) {
			t.Errorf("%s = %v, want NaN", name, got)
		}
		return
	}
	if math.Abs(got-want) > eps {
		t.Errorf("%s = %v, want %v (diff %g)", name, got, want, math.Abs(got-want))
	}
}

func TestTVMCore(t *testing.T) {
	// Future value of an ordinary annuity of -100 per period, 5% for 10 periods.
	approx(t, FV(0.05, 10, -100, 0, EndOfPeriod), 1257.789253554883, 1e-6, "FV")
	// Present value equals the annuity factor times payment.
	approx(t, PV(0.05, 10, -100, 0, EndOfPeriod), 772.1734929184818, 1e-6, "PV")
	// Level payment on a 1000 loan, 5% over 10 periods (Excel PMT).
	approx(t, PMT(0.05, 10, 1000, 0, EndOfPeriod), -129.50457496545658, 1e-6, "PMT")
	// Periods to repay a 1000 loan with 100 payments at 5%.
	approx(t, NPER(0.05, -100, 1000, 0, EndOfPeriod), 14.206699082890463, 1e-6, "NPER")
	// RATE inverts PMT.
	approx(t, RATE(10, -129.50457496545658, 1000, 0, EndOfPeriod, 0.1), 0.05, 1e-7, "RATE")

	// Zero-rate degenerate cases.
	approx(t, FV(0, 10, -100, 0, EndOfPeriod), 1000, 1e-9, "FV rate=0")
	approx(t, PMT(0, 10, 1000, 0, EndOfPeriod), -100, 1e-9, "PMT rate=0")
	approx(t, NPER(0, -100, 1000, 0, EndOfPeriod), 10, 1e-9, "NPER rate=0")

	// Annuity due is ordinary annuity times (1+rate).
	approx(t, FVAnnuityDue(100, 0.05, 10), FVAnnuity(100, 0.05, 10)*1.05, 1e-9, "FVAnnuityDue")
	approx(t, PVAnnuityDue(100, 0.05, 10), PVAnnuity(100, 0.05, 10)*1.05, 1e-9, "PVAnnuityDue")
}

func TestLumpAndAnnuity(t *testing.T) {
	approx(t, FVLumpSum(1000, 0.05, 10), 1628.894626777442, 1e-6, "FVLumpSum")
	approx(t, PVLumpSum(1628.894626777442, 0.05, 10), 1000, 1e-6, "PVLumpSum")
	approx(t, FVAnnuity(100, 0.05, 10), 1257.789253554883, 1e-6, "FVAnnuity")
	approx(t, PVAnnuity(100, 0.05, 10), 772.1734929184818, 1e-6, "PVAnnuity")
	approx(t, PVPerpetuity(100, 0.05), 2000, 1e-9, "PVPerpetuity")
	approx(t, PVGrowingPerpetuity(100, 0.05, 0.02), 3333.3333333333, 1e-6, "PVGrowingPerpetuity")
	approx(t, PVGrowingAnnuity(100, 0.05, 0.05, 10), 952.3809523809, 1e-6, "PVGrowingAnnuity (rate==growth)")
	approx(t, DiscountFactor(0.05, 10), 0.6139132535407591, 1e-9, "DiscountFactor")
	approx(t, AnnuityFactor(0.05, 10), 7.721734929184818, 1e-9, "AnnuityFactor")
	approx(t, AnnuityFactor(0, 10), 10, 1e-9, "AnnuityFactor rate=0")
}

func TestContinuousAndRateConv(t *testing.T) {
	approx(t, FVContinuous(1000, 0.05, 3), 1000*math.Exp(0.15), 1e-9, "FVContinuous")
	approx(t, PVContinuous(1000, 0.05, 3), 1000*math.Exp(-0.15), 1e-9, "PVContinuous")
	approx(t, EffectiveRate(0.12, 12), 0.12682503013196977, 1e-9, "EffectiveRate")
	approx(t, NominalRate(0.12682503013196977, 12), 0.12, 1e-9, "NominalRate")
	approx(t, EffectiveFromContinuous(0.10), math.Expm1(0.10), 1e-12, "EffectiveFromContinuous")
	approx(t, ContinuousFromEffective(math.Expm1(0.10)), 0.10, 1e-12, "ContinuousFromEffective")
}

func TestCashflow(t *testing.T) {
	cf := []float64{-1000, 500, 500, 500}
	approx(t, NPV(0.1, cf), 243.42599549211099, 1e-6, "NPV")

	irr := IRR(cf, 0.1)
	approx(t, NPV(irr, cf), 0, 1e-8, "NPV at IRR")
	approx(t, irr, 0.23375192852825868, 1e-6, "IRR")

	irrb := IRRBisection(cf)
	approx(t, irrb, irr, 1e-6, "IRRBisection matches IRR")

	mirr := MIRR([]float64{-120000, 39000, 30000, 21000, 37000, 46000}, 0.10, 0.12)
	approx(t, mirr, 0.12609413036, 1e-6, "MIRR")

	approx(t, CAGR(100, 200, 10), math.Pow(2, 0.1)-1, 1e-9, "CAGR")

	times := []float64{0, 1, 2, 3}
	approx(t, XNPV(0.1, cf, times), NPV(0.1, cf), 1e-9, "XNPV matches NPV")
	approx(t, XIRR(cf, times, 0.1), irr, 1e-6, "XIRR matches IRR")
}

func TestBonds(t *testing.T) {
	// Semiannual 5% coupon, 6% yield, 10 years.
	price := BondPrice(1000, 0.05, 0.06, 10, 2)
	approx(t, price, 925.6126251217, 1e-6, "BondPrice")
	approx(t, BondYTM(1000, 0.05, price, 10, 2), 0.06, 1e-6, "BondYTM")

	// Par bond: coupon == yield => price == face.
	approx(t, BondPrice(1000, 0.05, 0.05, 10, 2), 1000, 1e-6, "BondPrice par")

	approx(t, ZeroCouponPrice(1000, 0.05, 10, 1), 1000/math.Pow(1.05, 10), 1e-9, "ZeroCouponPrice")
	approx(t, CurrentYield(1000, 0.05, 950), 50.0/950, 1e-9, "CurrentYield")
	approx(t, AccruedInterest(1000, 0.06, 2, 90, 180), 15, 1e-9, "AccruedInterest")

	// Zero-coupon: Macaulay duration equals maturity.
	approx(t, MacaulayDuration(1000, 0, 0.05, 5, 1), 5, 1e-9, "MacaulayDuration zero-coupon")
	approx(t, ModifiedDuration(1000, 0, 0.05, 5, 1), 5/1.05, 1e-9, "ModifiedDuration zero-coupon")
	// Zero-coupon annual convexity = n(n+1)/(1+i)^2.
	approx(t, Convexity(1000, 0, 0.05, 5, 1), 30/(1.05*1.05), 1e-9, "Convexity zero-coupon")

	// Coupon bond: duration is below maturity.
	md := MacaulayDuration(1000, 0.05, 0.06, 10, 2)
	if md <= 0 || md >= 10 {
		t.Errorf("MacaulayDuration coupon = %v, want in (0,10)", md)
	}
}

func TestBlackScholes(t *testing.T) {
	s, k, tt, r, q, sig := 100.0, 100.0, 1.0, 0.05, 0.0, 0.2
	call := BlackScholesCall(s, k, tt, r, q, sig)
	put := BlackScholesPut(s, k, tt, r, q, sig)
	approx(t, call, 10.450583572185565, 1e-6, "BlackScholesCall")
	approx(t, put, 5.573526022256971, 1e-6, "BlackScholesPut")
	// Put-call parity: C - P = S e^-qT - K e^-rT.
	approx(t, call-put, s*math.Exp(-q*tt)-k*math.Exp(-r*tt), 1e-9, "put-call parity")
	approx(t, BlackScholesPrice(Call, s, k, tt, r, q, sig), call, 1e-12, "BlackScholesPrice call")
	approx(t, BlackScholesPrice(Put, s, k, tt, r, q, sig), put, 1e-12, "BlackScholesPrice put")

	approx(t, Delta(Call, s, k, tt, r, q, sig), 0.6368306511756191, 1e-6, "Delta call")
	approx(t, Delta(Put, s, k, tt, r, q, sig), 0.6368306511756191-1, 1e-6, "Delta put")
	approx(t, Gamma(s, k, tt, r, q, sig), 0.018762017345846895, 1e-6, "Gamma")
	approx(t, Vega(s, k, tt, r, q, sig), 37.52403469169379, 1e-6, "Vega")
	approx(t, Theta(Call, s, k, tt, r, q, sig), -6.414027546438197, 1e-6, "Theta call")
	approx(t, Rho(Call, s, k, tt, r, q, sig), 53.232481545376345, 1e-6, "Rho call")

	g := ComputeGreeks(Call, s, k, tt, r, q, sig)
	approx(t, g.Delta, Delta(Call, s, k, tt, r, q, sig), 1e-12, "Greeks.Delta")
	approx(t, g.Vega, Vega(s, k, tt, r, q, sig), 1e-12, "Greeks.Vega")

	iv := ImpliedVolatility(Call, call, s, k, tt, r, q)
	approx(t, iv, 0.2, 1e-6, "ImpliedVolatility")
}

func TestLoan(t *testing.T) {
	sched := AmortizationSchedule(1000, 0.1, 3)
	if len(sched) != 3 {
		t.Fatalf("schedule length = %d, want 3", len(sched))
	}
	approx(t, sched[0].Payment, 402.11480362537765, 1e-6, "payment")
	approx(t, sched[0].Interest, 100, 1e-6, "period 1 interest")
	approx(t, sched[0].Principal, 302.11480362537765, 1e-6, "period 1 principal")
	approx(t, sched[2].Balance, 0, 1e-9, "final balance")

	// Principal repaid across all rows equals the original loan.
	var totPrin float64
	for _, row := range sched {
		totPrin += row.Principal
	}
	approx(t, totPrin, 1000, 1e-6, "sum of principal")

	approx(t, RemainingBalance(1000, 0.1, 3, 1), 697.8851963746223, 1e-6, "RemainingBalance")
	approx(t, RemainingBalance(1000, 0.1, 3, 3), 0, 1e-9, "RemainingBalance final")
	approx(t, TotalInterest(1000, 0.1, 3), 402.11480362537765*3-1000, 1e-6, "TotalInterest")
}

func TestDepreciation(t *testing.T) {
	approx(t, StraightLineDepreciation(10000, 1000, 5), 1800, 1e-9, "StraightLine")
	sl := StraightLineSchedule(10000, 1000, 5)
	approx(t, sl[0], 1800, 1e-9, "StraightLineSchedule[0]")

	ddb := DoubleDecliningSchedule(10000, 1000, 5)
	want := []float64{4000, 2400, 1440, 864, 296}
	for i, w := range want {
		approx(t, ddb[i], w, 1e-6, "DDB schedule")
	}
	approx(t, DoubleDecliningBalance(10000, 1000, 1, 5), 4000, 1e-6, "DoubleDecliningBalance p1")
	approx(t, DoubleDecliningBalance(10000, 1000, 5, 5), 296, 1e-6, "DoubleDecliningBalance p5")

	syd := SumOfYearsDigitsSchedule(10000, 1000, 5)
	sydWant := []float64{3000, 2400, 1800, 1200, 600}
	var sydSum float64
	for i, w := range sydWant {
		approx(t, syd[i], w, 1e-6, "SYD schedule")
		sydSum += syd[i]
	}
	approx(t, sydSum, 9000, 1e-6, "SYD schedule sum")
	approx(t, SumOfYearsDigits(10000, 1000, 1, 5), 3000, 1e-6, "SumOfYearsDigits p1")

	dbSum := 0.0
	for p := 1; p <= 5; p++ {
		dbSum += DecliningBalance(10000, 1000, 0.3, p, 5)
	}
	if dbSum <= 0 || dbSum > 9000+1e-6 {
		t.Errorf("DecliningBalance total = %v, want in (0, 9000]", dbSum)
	}
}

func TestPerformance(t *testing.T) {
	rets := []float64{0.1, 0.2, 0.3}
	approx(t, Volatility(rets), 0.1, 1e-9, "Volatility")
	approx(t, SharpeRatio(rets, 0.05), 1.5, 1e-9, "SharpeRatio")
	approx(t, Annualize(0.05, 12), 0.05*math.Sqrt(12), 1e-9, "Annualize")

	approx(t, HoldingPeriodReturn([]float64{0.1, 0.1}), 0.21, 1e-9, "HoldingPeriodReturn")
	approx(t, AnnualizedReturn([]float64{0.1, 0.1}), 0.1, 1e-9, "AnnualizedReturn")
	approx(t, SimpleReturn(100, 110), 0.1, 1e-9, "SimpleReturn")
	approx(t, LogReturn(100, 110), math.Log(1.1), 1e-9, "LogReturn")

	approx(t, SortinoRatio([]float64{0.1, -0.05, 0.2}, 0), 2.886751345948, 1e-6, "SortinoRatio")
}

func TestNaNGuards(t *testing.T) {
	if !math.IsNaN(PVPerpetuity(100, 0)) {
		t.Error("PVPerpetuity(100,0) should be NaN")
	}
	if !math.IsNaN(PVGrowingPerpetuity(100, 0.02, 0.05)) {
		t.Error("PVGrowingPerpetuity with rate<=growth should be NaN")
	}
	if !math.IsNaN(CAGR(-1, 100, 5)) {
		t.Error("CAGR with begin<=0 should be NaN")
	}
	if !math.IsNaN(Volatility([]float64{1})) {
		t.Error("Volatility of single value should be NaN")
	}
	if !math.IsNaN(XNPV(0.1, []float64{1, 2}, []float64{0})) {
		t.Error("XNPV with mismatched lengths should be NaN")
	}
}

func BenchmarkBondYTM(b *testing.B) {
	price := BondPrice(1000, 0.05, 0.06, 30, 2)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = BondYTM(1000, 0.05, price, 30, 2)
	}
}
