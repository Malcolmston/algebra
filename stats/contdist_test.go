package stats

import (
	"math"
	"testing"
)

func TestBetaContDist(t *testing.T) {
	b := Beta{Alpha: 2, Beta: 3}
	// B(2,3) = 1/12, so PDF(0.5) = 0.5·0.25·12 = 1.5.
	approx(t, b.PDF(0.5), 1.5, tol, "beta pdf")
	approx(t, b.CDF(0.5), 0.6875, 1e-9, "beta cdf")
	approx(t, b.Mean(), 0.4, tol, "beta mean")
	approx(t, b.Variance(), 0.04, tol, "beta var")
	// Quantile inverts CDF.
	approx(t, b.Quantile(0.6875), 0.5, 1e-6, "beta quantile")
	approx(t, b.CDF(b.Quantile(0.3)), 0.3, 1e-6, "beta roundtrip")
	// Out-of-support and bad-parameter handling.
	approx(t, b.PDF(-0.1), 0, 0, "beta pdf below")
	approx(t, b.PDF(1.1), 0, 0, "beta pdf above")
	approx(t, b.CDF(0), 0, 0, "beta cdf 0")
	approx(t, b.CDF(1), 1, 0, "beta cdf 1")
	if !math.IsNaN(b.Quantile(1.5)) {
		t.Fatalf("beta quantile out-of-range: want NaN")
	}
	if !math.IsNaN((Beta{Alpha: 0, Beta: 3}).PDF(0.5)) {
		t.Fatalf("beta bad param: want NaN")
	}
}

func TestFDistContDist(t *testing.T) {
	approx(t, FDist{D1: 1, D2: 1}.CDF(1), 0.5, 1e-9, "F(1,1) cdf1")
	f := FDist{D1: 3, D2: 10}
	approx(t, f.Mean(), 1.25, tol, "F mean")
	// Variance = 2·100·11 / (3·64·6) = 2200/1152.
	approx(t, f.Variance(), 2200.0/1152.0, 1e-9, "F var")
	// Quantile inverts CDF; known critical value F(3,10) at 0.95 is ~3.7083.
	q := f.Quantile(0.95)
	approx(t, f.CDF(q), 0.95, 1e-6, "F roundtrip")
	approx(t, q, 3.7083, 1e-3, "F crit")
	approx(t, f.CDF(-1), 0, 0, "F cdf below")
	// Mean/variance existence boundaries.
	if !math.IsNaN((FDist{D1: 3, D2: 2}).Mean()) {
		t.Fatalf("F mean D2=2: want NaN")
	}
	if !math.IsInf((FDist{D1: 3, D2: 3}).Variance(), 1) {
		t.Fatalf("F var 2<D2<=4: want +Inf")
	}
	if !math.IsNaN((FDist{D1: 3, D2: 2}).Variance()) {
		t.Fatalf("F var D2<=2: want NaN")
	}
}

func TestLogNormalContDist(t *testing.T) {
	l := LogNormal{Mu: 0, Sigma: 1}
	approx(t, l.CDF(1), 0.5, tol, "lognorm cdf1")
	approx(t, l.PDF(1), 1/math.Sqrt(2*math.Pi), tol, "lognorm pdf1")
	approx(t, l.Mean(), math.Exp(0.5), tol, "lognorm mean")
	approx(t, l.Variance(), math.Exp(1)*(math.Exp(1)-1), tol, "lognorm var")
	approx(t, l.Quantile(0.5), 1, 1e-9, "lognorm median")
	approx(t, l.CDF(l.Quantile(0.9)), 0.9, 1e-9, "lognorm roundtrip")
	approx(t, l.PDF(-1), 0, 0, "lognorm pdf neg")
	approx(t, l.CDF(0), 0, 0, "lognorm cdf 0")
	if !math.IsNaN((LogNormal{Mu: 0, Sigma: 0}).PDF(1)) {
		t.Fatalf("lognorm bad sigma: want NaN")
	}
}

func TestWeibullContDist(t *testing.T) {
	// Shape=1 reduces to Exponential(rate 1/Scale).
	w := Weibull{Shape: 1, Scale: 1}
	approx(t, w.CDF(1), 1-math.Exp(-1), tol, "weibull cdf")
	approx(t, w.PDF(0.5), math.Exp(-0.5), tol, "weibull pdf")
	approx(t, w.Mean(), 1, tol, "weibull mean")
	// Weibull(2,1): mean = Gamma(1.5) = sqrt(pi)/2, var = 1 - pi/4.
	w2 := Weibull{Shape: 2, Scale: 1}
	approx(t, w2.Mean(), math.Sqrt(math.Pi)/2, 1e-9, "weibull2 mean")
	approx(t, w2.Variance(), 1-math.Pi/4, 1e-9, "weibull2 var")
	approx(t, w2.CDF(w2.Quantile(0.7)), 0.7, 1e-9, "weibull roundtrip")
	approx(t, w.CDF(-1), 0, 0, "weibull cdf neg")
	if !math.IsNaN((Weibull{Shape: 0, Scale: 1}).PDF(1)) {
		t.Fatalf("weibull bad shape: want NaN")
	}
}

func BenchmarkBetaContCDF(b *testing.B) {
	d := Beta{Alpha: 2.5, Beta: 3.5}
	for i := 0; i < b.N; i++ {
		_ = d.CDF(0.4)
	}
}

func BenchmarkBetaContQuantile(b *testing.B) {
	d := Beta{Alpha: 2.5, Beta: 3.5}
	for i := 0; i < b.N; i++ {
		_ = d.Quantile(0.4)
	}
}

func BenchmarkFDistContCDF(b *testing.B) {
	d := FDist{D1: 5, D2: 12}
	for i := 0; i < b.N; i++ {
		_ = d.CDF(2.3)
	}
}

func BenchmarkFDistContQuantile(b *testing.B) {
	d := FDist{D1: 5, D2: 12}
	for i := 0; i < b.N; i++ {
		_ = d.Quantile(0.95)
	}
}
