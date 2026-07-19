package reliability_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/malcolmston/algebra/reliability"
)

const tol = 1e-9

func approx(t *testing.T, name string, got, want, eps float64) {
	t.Helper()
	if math.IsNaN(want) {
		if !math.IsNaN(got) {
			t.Errorf("%s = %v, want NaN", name, got)
		}
		return
	}
	if math.Abs(got-want) > eps {
		t.Errorf("%s = %.12g, want %.12g (diff %.3g)", name, got, want, math.Abs(got-want))
	}
}

func TestExponential(t *testing.T) {
	lambda := 0.25
	approx(t, "R(4)", reliability.ExponentialReliability(4, lambda), math.Exp(-1), tol)
	approx(t, "F(4)", reliability.ExponentialCDF(4, lambda), 1-math.Exp(-1), tol)
	approx(t, "pdf(0)", reliability.ExponentialPDF(0, lambda), lambda, tol)
	approx(t, "hazard", reliability.ExponentialHazard(10, lambda), lambda, tol)
	approx(t, "cumhaz", reliability.ExponentialCumulativeHazard(4, lambda), 1, tol)
	approx(t, "mean", reliability.ExponentialMean(lambda), 4, tol)
	approx(t, "mttf", reliability.ExponentialMTTF(lambda), 4, tol)
	approx(t, "var", reliability.ExponentialVariance(lambda), 16, tol)
	approx(t, "std", reliability.ExponentialStdDev(lambda), 4, tol)
	approx(t, "median", reliability.ExponentialMedian(lambda), math.Ln2/lambda, tol)
	approx(t, "quantile", reliability.ExponentialQuantile(1-math.Exp(-1), lambda), 4, 1e-9)
	// memoryless conditional reliability
	approx(t, "cond", reliability.ExponentialConditionalReliability(100, 4, lambda), math.Exp(-1), tol)
	approx(t, "mrl", reliability.ExponentialMeanResidualLife(100, lambda), 4, tol)
	approx(t, "rateFromR", reliability.ExponentialRateFromReliability(4, math.Exp(-1)), lambda, tol)
	approx(t, "rateFromMTTF", reliability.ExponentialRateFromMTTF(4), lambda, tol)
	// domain errors
	if !math.IsNaN(reliability.ExponentialReliability(1, -1)) {
		t.Error("expected NaN for negative lambda")
	}
}

func TestWeibull(t *testing.T) {
	k, lam := 2.0, 100.0
	approx(t, "R(50)", reliability.WeibullReliability(50, k, lam), 0.7788007830714049, 1e-12)
	approx(t, "F(50)", reliability.WeibullCDF(50, k, lam), 1-0.7788007830714049, 1e-12)
	approx(t, "cumhaz(50)", reliability.WeibullCumulativeHazard(50, k, lam), 0.25, tol)
	approx(t, "hazard(50)", reliability.WeibullHazard(50, k, lam), (k/lam)*math.Pow(0.5, k-1), tol)
	approx(t, "mean", reliability.WeibullMean(k, lam), 88.62269254527581, 1e-9)
	approx(t, "var", reliability.WeibullVariance(k, lam), 2146.018366025516, 1e-6)
	approx(t, "std", reliability.WeibullStdDev(k, lam), math.Sqrt(2146.018366025516), 1e-7)
	approx(t, "median", reliability.WeibullMedian(k, lam), lam*math.Pow(math.Ln2, 1/k), tol)
	approx(t, "mode", reliability.WeibullMode(k, lam), lam*math.Pow(0.5, 0.5), tol)
	approx(t, "charlife", reliability.WeibullCharacteristicLife(lam), 100, tol)
	// quantile round-trips against CDF
	q := reliability.WeibullQuantile(0.3, k, lam)
	approx(t, "quantile-roundtrip", reliability.WeibullCDF(q, k, lam), 0.3, 1e-9)
	// exponential special case k=1
	approx(t, "weib-exp", reliability.WeibullReliability(50, 1, lam), math.Exp(-0.5), tol)
	// conditional reliability
	approx(t, "cond", reliability.WeibullConditionalReliability(50, 0, k, lam), 1, tol)
	// MRL by numeric integration, cross-checked against mean at t=0
	approx(t, "mrl0", reliability.WeibullMeanResidualLife(0, k, lam), reliability.WeibullMean(k, lam), 1e-4)
}

func TestLognormal(t *testing.T) {
	mu, sig := 3.0, 0.5
	approx(t, "cdf(20)", reliability.LognormalCDF(20, mu, sig), 0.49659488830502574, 1e-9)
	approx(t, "R(20)", reliability.LognormalReliability(20, mu, sig), 1-0.49659488830502574, 1e-9)
	approx(t, "mean", reliability.LognormalMean(mu, sig), 22.75989509352673, 1e-8)
	approx(t, "median", reliability.LognormalMedian(mu, sig), math.Exp(mu), tol)
	approx(t, "mode", reliability.LognormalMode(mu, sig), math.Exp(mu-sig*sig), tol)
	approx(t, "var", reliability.LognormalVariance(mu, sig), (math.Exp(sig*sig)-1)*math.Exp(2*mu+sig*sig), 1e-6)
	// quantile round trip
	q := reliability.LognormalQuantile(0.4, mu, sig)
	approx(t, "quantile-roundtrip", reliability.LognormalCDF(q, mu, sig), 0.4, 1e-9)
	if h := reliability.LognormalHazard(20, mu, sig); !(h > 0) {
		t.Errorf("lognormal hazard should be positive, got %v", h)
	}
}

func TestGamma(t *testing.T) {
	k, theta := 3.0, 2.0
	approx(t, "cdf(6)", reliability.GammaCDF(6, k, theta), 0.5768099188731566, 1e-9)
	approx(t, "R(6)", reliability.GammaReliability(6, k, theta), 1-0.5768099188731566, 1e-9)
	approx(t, "mean", reliability.GammaMean(k, theta), 6, tol)
	approx(t, "var", reliability.GammaVariance(k, theta), 12, tol)
	approx(t, "std", reliability.GammaStdDev(k, theta), math.Sqrt(12), tol)
	approx(t, "mode", reliability.GammaMode(k, theta), 4, tol)
	// quantile inverts CDF
	q := reliability.GammaQuantile(0.5768099188731566, k, theta)
	approx(t, "quantile", q, 6, 1e-6)
	// exponential special case k=1: R(t)=exp(-t/theta)
	approx(t, "gamma-exp", reliability.GammaReliability(2, 1, theta), math.Exp(-1), 1e-9)
	// pdf integrates via hazard = pdf/R consistency
	approx(t, "hazard", reliability.GammaHazard(6, k, theta),
		reliability.GammaPDF(6, k, theta)/reliability.GammaReliability(6, k, theta), tol)
}

func TestRayleigh(t *testing.T) {
	sig := 10.0
	approx(t, "R", reliability.RayleighReliability(10, sig), math.Exp(-0.5), tol)
	approx(t, "hazard", reliability.RayleighHazard(10, sig), 10.0/(sig*sig), tol)
	approx(t, "cumhaz", reliability.RayleighCumulativeHazard(10, sig), 0.5, tol)
	approx(t, "mean", reliability.RayleighMean(sig), 12.533141373155, 1e-9)
	approx(t, "var", reliability.RayleighVariance(sig), (2-math.Pi/2)*100, 1e-9)
	approx(t, "median", reliability.RayleighMedian(sig), sig*math.Sqrt(2*math.Ln2), tol)
	approx(t, "mode", reliability.RayleighMode(sig), sig, tol)
	// Rayleigh is Weibull(k=2, lambda=sig*sqrt2)
	lam := sig * math.Sqrt2
	approx(t, "rayl=weib", reliability.RayleighReliability(15, sig),
		reliability.WeibullReliability(15, 2, lam), 1e-12)
	q := reliability.RayleighQuantile(0.6, sig)
	approx(t, "quantile-roundtrip", reliability.RayleighCDF(q, sig), 0.6, 1e-9)
}

func TestNormal(t *testing.T) {
	mu, sig := 1000.0, 100.0
	approx(t, "cdf", reliability.NormalCDF(1000, mu, sig), 0.5, tol)
	approx(t, "R", reliability.NormalReliability(1000, mu, sig), 0.5, tol)
	approx(t, "probit975", reliability.StandardNormalQuantile(0.975), 1.9599639845400536, 1e-9)
	approx(t, "phi0", reliability.StandardNormalCDF(0), 0.5, tol)
	q := reliability.NormalQuantile(0.9, mu, sig)
	approx(t, "quantile-roundtrip", reliability.NormalCDF(q, mu, sig), 0.9, 1e-9)
	approx(t, "mttf", reliability.NormalMTTF(mu, sig), 1000, tol)
}

func TestGompertz(t *testing.T) {
	a, b := 0.01, 0.1
	approx(t, "R(10)", reliability.GompertzReliability(10, a, b), 0.8421238520626434, 1e-12)
	approx(t, "hazard", reliability.GompertzHazard(10, a, b), a*math.Exp(b*10), tol)
	approx(t, "pdf", reliability.GompertzPDF(10, a, b),
		reliability.GompertzHazard(10, a, b)*reliability.GompertzReliability(10, a, b), tol)
	q := reliability.GompertzQuantile(0.3, a, b)
	approx(t, "quantile-roundtrip", reliability.GompertzCDF(q, a, b), 0.3, 1e-9)
	// Makeham reduces to Gompertz when lambda=0
	approx(t, "makeham0", reliability.GompertzMakehamReliability(10, 0, a, b),
		reliability.GompertzReliability(10, a, b), 1e-12)
	// added constant term lowers reliability
	if reliability.GompertzMakehamReliability(10, 0.02, a, b) >= reliability.GompertzReliability(10, a, b) {
		t.Error("Makeham term should reduce reliability")
	}
}

func TestBathtubAndHjorth(t *testing.T) {
	// bathtub hazard is U-shaped: decreasing then increasing
	a, b, c, e := 1.0, 1.0, 0.05, 0.02
	h0 := reliability.BathtubHazard(0, a, b, c, e)
	hmid := reliability.BathtubHazard(3, a, b, c, e)
	hlate := reliability.BathtubHazard(30, a, b, c, e)
	if !(h0 > hmid && hlate > hmid) {
		t.Errorf("bathtub not U-shaped: %.4f %.4f %.4f", h0, hmid, hlate)
	}
	// reliability consistent with cumulative hazard
	approx(t, "bathtub-R", reliability.BathtubReliability(3, a, b, c, e),
		math.Exp(-reliability.BathtubCumulativeHazard(3, a, b, c, e)), 1e-12)
	approx(t, "bathtub-pdf", reliability.BathtubPDF(3, a, b, c, e),
		reliability.BathtubHazard(3, a, b, c, e)*reliability.BathtubReliability(3, a, b, c, e), tol)

	approx(t, "hjorth-R", reliability.HjorthReliability(5, 0.01, 0.5, 0.2), 0.15600488604842286, 1e-12)
	approx(t, "hjorth-cdf", reliability.HjorthCDF(5, 0.01, 0.5, 0.2), 1-0.15600488604842286, 1e-12)
	approx(t, "hjorth-pdf", reliability.HjorthPDF(5, 0.01, 0.5, 0.2),
		reliability.HjorthHazard(5, 0.01, 0.5, 0.2)*reliability.HjorthReliability(5, 0.01, 0.5, 0.2), tol)
}

func TestSystemReliability(t *testing.T) {
	approx(t, "series", reliability.SeriesReliability(0.9, 0.8, 0.95), 0.9*0.8*0.95, tol)
	approx(t, "parallel", reliability.ParallelReliability(0.9, 0.8), 1-0.1*0.2, tol)
	approx(t, "series-fail", reliability.SeriesFailureProbability(0.9, 0.9), 1-0.81, tol)
	approx(t, "parallel-fail", reliability.ParallelFailureProbability(0.9, 0.8), 0.1*0.2, tol)
	approx(t, "series-id", reliability.SeriesReliabilityIdentical(0.99, 5), math.Pow(0.99, 5), tol)
	approx(t, "parallel-id", reliability.ParallelReliabilityIdentical(0.9, 3), 1-math.Pow(0.1, 3), tol)
	approx(t, "kofn", reliability.KofNReliability(2, 3, 0.9), 0.972, 1e-12)
	approx(t, "kofn-all", reliability.KofNReliability(3, 3, 0.9), math.Pow(0.9, 3), tol)
	approx(t, "kofn-any", reliability.KofNReliability(1, 3, 0.9), reliability.ParallelReliabilityIdentical(0.9, 3), tol)

	// general k-of-n with identical components matches the closed form
	rg, err := reliability.KofNReliabilityGeneral(2, []float64{0.9, 0.9, 0.9})
	if err != nil {
		t.Fatal(err)
	}
	approx(t, "kofn-general", rg, 0.972, 1e-12)
	// heterogeneous general k-of-n, 2-of-3 with distinct reliabilities
	rh, _ := reliability.KofNReliabilityGeneral(2, []float64{0.9, 0.8, 0.7})
	// exact: sum of exactly-2 and exactly-3 working
	p := []float64{0.9, 0.8, 0.7}
	e2 := p[0]*p[1]*(1-p[2]) + p[0]*(1-p[1])*p[2] + (1-p[0])*p[1]*p[2]
	e3 := p[0] * p[1] * p[2]
	approx(t, "kofn-het", rh, e2+e3, 1e-12)

	approx(t, "bridge", reliability.BridgeReliability(0.9, 0.9, 0.9, 0.9, 0.9), 0.97848, 1e-9)
	// bridge with perfectly reliable bridge reduces to two parallel pairs in series
	approx(t, "bridge-perfect", reliability.BridgeReliability(0.8, 0.8, 0.8, 0.8, 1.0),
		reliability.ParallelReliability(0.8, 0.8)*reliability.ParallelReliability(0.8, 0.8), 1e-12)

	approx(t, "series-parallel", reliability.SeriesParallelReliability([][]float64{{0.9, 0.9}, {0.8}}),
		(1-0.01)*0.8, tol)
	approx(t, "parallel-series", reliability.ParallelSeriesReliability([][]float64{{0.9, 0.9}, {0.8, 0.8}}),
		1-(1-0.81)*(1-0.64), tol)
}

func TestStandbyAndRates(t *testing.T) {
	approx(t, "standby", reliability.StandbyRedundancyReliability(5, 0.1, 2), 0.9097959895689501, 1e-12)
	approx(t, "standby1", reliability.StandbyRedundancyReliability(5, 0.1, 1),
		reliability.ExponentialReliability(5, 0.1), 1e-12)
	approx(t, "standby-mttf", reliability.StandbyRedundancyMTTF(0.1, 3), 30, tol)
	approx(t, "series-rate", reliability.SeriesSystemFailureRate(0.1, 0.2, 0.3), 0.6, tol)
	approx(t, "series-mttf", reliability.SeriesSystemMTTF(0.1, 0.2, 0.3), 1.0/0.6, tol)
	approx(t, "parallel-mttf", reliability.ParallelSystemMTTF(0.1, 2), (1+0.5)/0.1, tol)
}

func TestRedundancy(t *testing.T) {
	if got := reliability.MinParallelForReliability(0.999, 0.9); got != 3 {
		t.Errorf("MinParallelForReliability = %d, want 3", got)
	}
	if got := reliability.MaxSeriesForReliability(0.9, 0.99); got != 10 {
		t.Errorf("MaxSeriesForReliability = %d, want 10", got)
	}
	approx(t, "gain", reliability.RedundancyGain(0.9), 0.09, tol)
	approx(t, "birnbaum", reliability.ReliabilityImportanceBirnbaum(0.99, 0.9), 0.09, tol)
	approx(t, "criticality", reliability.CriticalityImportance(0.5, 0.9, 0.99),
		0.5*(1-0.9)/(1-0.99), tol)
	approx(t, "alloc-series", reliability.AllocateReliabilitySeriesEqual(0.9, 3), math.Pow(0.9, 1.0/3), tol)
	approx(t, "alloc-rate", reliability.AllocateFailureRateSeriesEqual(0.6, 3), 0.2, tol)
}

func TestMetrics(t *testing.T) {
	approx(t, "failrate", reliability.FailureRate(0.2, 0.8), 0.25, tol)
	approx(t, "H-from-R", reliability.HazardFromReliability(math.Exp(-1)), 1, tol)
	approx(t, "R-from-H", reliability.ReliabilityFromHazard(1), math.Exp(-1), tol)
	approx(t, "fit", reliability.FIT(1e-6), 1000, tol)
	approx(t, "fit-round", reliability.FITToRate(reliability.FIT(2e-7)), 2e-7, 1e-18)
	approx(t, "fpmh", reliability.FailuresPerMillionHours(1e-6), 1, tol)
	approx(t, "mtbf-rate", reliability.MTBFFromFailureRate(0.001), 1000, tol)
	approx(t, "rate-mtbf", reliability.FailureRateFromMTBF(1000), 0.001, tol)
	approx(t, "mtbf", reliability.MTBF(900, 100), 1000, tol)
	approx(t, "mtbf-data", reliability.MTBFFromData(10000, 5), 2000, tol)
	approx(t, "mttr", reliability.MTTR(0.2), 5, tol)
	approx(t, "reprate", reliability.RepairRateFromMTTR(5), 0.2, tol)

	m, err := reliability.MTTFFromData([]float64{10, 20, 30})
	if err != nil {
		t.Fatal(err)
	}
	approx(t, "mttf-data", m, 20, tol)

	// MTTF from reliability integral for an exponential should be 1/lambda
	approx(t, "mttf-integral",
		reliability.MTTFFromReliability(func(x float64) float64 { return math.Exp(-0.5 * x) }, 2),
		2, 1e-4)

	approx(t, "medrank", reliability.MedianRank(3, 10), (3-0.3)/(10+0.4), tol)
	approx(t, "meanrank", reliability.MeanRank(3, 10), 3.0/11, tol)

	data := []float64{5, 10, 15, 20, 25}
	approx(t, "emp-R", reliability.EmpiricalReliability(data, 12), 3.0/5, tol)
	approx(t, "emp-F", reliability.EmpiricalCDF(data, 12), 2.0/5, tol)
	approx(t, "design-life", reliability.DesignLifeForReliability(math.Exp(-1), 0.25), 4, tol)
}

func TestAvailability(t *testing.T) {
	approx(t, "avail", reliability.Availability(900, 100), 0.9, tol)
	approx(t, "inherent", reliability.InherentAvailability(950, 50), 0.95, tol)
	approx(t, "operational", reliability.OperationalAvailability(90, 10), 0.9, tol)
	approx(t, "achieved", reliability.AchievedAvailability(80, 20), 0.8, tol)
	approx(t, "unavail", reliability.Unavailability(900, 100), 0.1, tol)
	approx(t, "steady", reliability.SteadyStateAvailability(0.01, 0.1), 0.1/0.11, tol)
	approx(t, "steady-alias", reliability.AvailabilityFromRates(0.01, 0.1), 0.1/0.11, tol)
	approx(t, "instant", reliability.InstantaneousAvailability(10, 0.01, 0.1), 0.9393519166998255, 1e-12)
	// instantaneous availability at t=0 is 1
	approx(t, "instant0", reliability.InstantaneousAvailability(0, 0.01, 0.1), 1, tol)
	// instantaneous tends to steady state for large t
	approx(t, "instant-inf", reliability.InstantaneousAvailability(1e6, 0.01, 0.1), 0.1/0.11, 1e-9)
	approx(t, "avail-series", reliability.AvailabilitySeriesSystem(0.9, 0.95), 0.855, tol)
	approx(t, "avail-parallel", reliability.AvailabilityParallelSystem(0.9, 0.9), 0.99, tol)
	approx(t, "avail-kofn", reliability.AvailabilityKofN(2, 3, 0.9), 0.972, 1e-12)
	approx(t, "downtime", reliability.ExpectedDowntime(0.99, 8760), 87.6, 1e-9)
	approx(t, "uptime", reliability.ExpectedUptime(0.99, 8760), 8672.4, 1e-9)
}

func TestKaplanMeier(t *testing.T) {
	times := []float64{1, 2, 3, 4, 5}
	events := []bool{true, false, true, false, true}
	km, err := reliability.KaplanMeier(times, events)
	if err != nil {
		t.Fatal(err)
	}
	if len(km.Times) != 3 {
		t.Fatalf("expected 3 event times, got %d", len(km.Times))
	}
	approx(t, "km-t0", km.Times[0], 1, tol)
	approx(t, "km-atrisk0", float64(km.AtRisk[0]), 5, tol)
	approx(t, "km-atrisk1", float64(km.AtRisk[1]), 3, tol)
	approx(t, "km-S0", km.Survival[0], 0.8, tol)
	approx(t, "km-S1", km.Survival[1], 0.5333333333333333, 1e-12)
	approx(t, "km-S2", km.Survival[2], 0, tol)
	approx(t, "km-var1", km.Variance[1], 0.06162962962962963, 1e-12)
	approx(t, "km-se1", km.StdErr[1], 0.24825315633367007, 1e-12)
	// step function
	approx(t, "km-at-2.5", km.SurvivalAt(2.5), 0.8, tol)
	approx(t, "km-at-0.5", km.SurvivalAt(0.5), 1, tol)
	approx(t, "km-at-3", km.SurvivalAt(3), 0.5333333333333333, 1e-12)
	// median reached at t=5 (first time S<=0.5; S(3)=0.533 is still above 0.5)
	approx(t, "km-median", km.MedianSurvival(), 5, tol)

	// complete data: survival matches empirical fractions
	km2, _ := reliability.KaplanMeier([]float64{1, 2, 3, 4, 5}, []bool{true, true, true, true, true})
	approx(t, "km-complete-S1", km2.Survival[0], 0.8, tol)
	approx(t, "km-complete-S3", km2.Survival[2], 0.4, tol)
	// RMST over [0,5] for complete uniform failures: 1*1+0.8*1+0.6*1+0.4*1+0.2*1 = 3
	approx(t, "km-rmst", km2.RestrictedMeanSurvival(5), 3, tol)
}

func TestNelsonAalen(t *testing.T) {
	times := []float64{1, 2, 3, 4, 5}
	events := []bool{true, false, true, false, true}
	na, err := reliability.NelsonAalen(times, events)
	if err != nil {
		t.Fatal(err)
	}
	approx(t, "na-H0", na.CumHazard[0], 0.2, tol)
	approx(t, "na-H1", na.CumHazard[1], 0.2+1.0/3, 1e-12)
	approx(t, "na-H2", na.CumHazard[2], 0.2+1.0/3+1, 1e-12)
	approx(t, "na-var0", na.Variance[0], 1.0/25, tol)
	approx(t, "na-at-2.5", na.CumulativeHazardAt(2.5), 0.2, tol)
	approx(t, "na-surv", na.SurvivalAt(1), math.Exp(-0.2), 1e-12)
}

func TestLogRank(t *testing.T) {
	g1t := []float64{1, 2}
	g1e := []bool{true, true}
	g2t := []float64{3, 4}
	g2e := []bool{true, true}
	res, err := reliability.LogRankTest(g1t, g1e, g2t, g2e)
	if err != nil {
		t.Fatal(err)
	}
	approx(t, "lr-O1", res.Observed1, 2, tol)
	approx(t, "lr-E1", res.Expected1, 0.8333333333333333, 1e-12)
	approx(t, "lr-V", res.Variance, 0.4722222222222222, 1e-12)
	approx(t, "lr-chi2", res.ChiSquare, 2.8823529411764715, 1e-12)
	approx(t, "lr-p", res.PValue, 0.08955507441364256, 1e-12)
	approx(t, "lr-hr", res.HazardRatio, 11.82940327795448, 1e-9)

	// identical groups -> no difference
	res2, _ := reliability.LogRankTest(g1t, g1e, g1t, g1e)
	approx(t, "lr-identical-chi2", res2.ChiSquare, 0, 1e-12)
	approx(t, "lr-identical-p", res2.PValue, 1, tol)
}

func TestHazardHelpers(t *testing.T) {
	approx(t, "hr", reliability.HazardRatio(0.02, 0.01), 2, tol)
	approx(t, "cox", reliability.ProportionalHazardsReliability(0.8, 2), 0.64, tol)
	approx(t, "chi2p", reliability.ChiSquarePValue1DF(2.8823529411764715), 0.08955507441364256, 1e-12)
	approx(t, "chi2p0", reliability.ChiSquarePValue1DF(0), 1, tol)
}

func TestErrorPaths(t *testing.T) {
	if _, err := reliability.KaplanMeier([]float64{1, 2}, []bool{true}); err == nil {
		t.Error("expected length-mismatch error")
	}
	if _, err := reliability.KaplanMeier([]float64{1, 2}, []bool{false, false}); err == nil {
		t.Error("expected no-failures error")
	}
	if _, err := reliability.KofNReliabilityGeneral(5, []float64{0.9, 0.9}); err == nil {
		t.Error("expected k>n error")
	}
	if _, err := reliability.MTTFFromData(nil); err == nil {
		t.Error("expected empty-data error")
	}
	if n := reliability.MinParallelForReliability(0.9, 1.5); n != -1 {
		t.Error("expected -1 for invalid reliability")
	}
}

// Example functions -------------------------------------------------------

func ExampleWeibullReliability() {
	// A component with Weibull shape 2 (wear-out) and scale 100 hours.
	r := reliability.WeibullReliability(50, 2, 100)
	fmt.Printf("R(50) = %.4f\n", r)
	// Output: R(50) = 0.7788
}

func ExampleSeriesReliability() {
	// Three independent components in series.
	fmt.Printf("%.4f\n", reliability.SeriesReliability(0.99, 0.98, 0.97))
	// Output: 0.9411
}

func ExampleKofNReliability() {
	// A 2-out-of-3 voting system of identical 0.9-reliable units.
	fmt.Printf("%.3f\n", reliability.KofNReliability(2, 3, 0.9))
	// Output: 0.972
}

func ExampleAvailability() {
	// Availability from a 900-hour MTBF and 100-hour MTTR.
	fmt.Printf("%.2f\n", reliability.Availability(900, 100))
	// Output: 0.90
}

func ExampleKaplanMeier() {
	times := []float64{1, 2, 3, 4, 5}
	events := []bool{true, false, true, false, true}
	km, _ := reliability.KaplanMeier(times, events)
	fmt.Printf("S(3) = %.4f\n", km.SurvivalAt(3))
	// Output: S(3) = 0.5333
}
