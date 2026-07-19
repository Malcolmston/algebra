package queueing

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool {
	if math.IsInf(a, 0) || math.IsInf(b, 0) {
		return a == b
	}
	return math.Abs(a-b) <= eps
}

func TestMM1(t *testing.T) {
	q, err := NewMM1(1, 2)
	if err != nil {
		t.Fatalf("NewMM1: %v", err)
	}
	tests := []struct {
		name string
		got  float64
		want float64
	}{
		{"Rho", q.Rho(), 0.5},
		{"P0", q.P0(), 0.5},
		{"P1", q.Pn(1), 0.25},
		{"PAtLeast2", q.PAtLeast(2), 0.25},
		{"L", q.L(), 1.0},
		{"Lq", q.Lq(), 0.5},
		{"W", q.W(), 1.0},
		{"Wq", q.Wq(), 0.5},
		{"VarianceN", q.VarianceN(), 2.0},
		{"WaitTail0", q.WaitTailProb(0), 1.0},
		{"WaitCDFinf", q.WaitCDF(1e9), 1.0},
		{"MeanBusyPeriod", q.MeanBusyPeriod(), 1.0},
	}
	for _, tc := range tests {
		if !approx(tc.got, tc.want, tol) {
			t.Errorf("MM1 %s = %v, want %v", tc.name, tc.got, tc.want)
		}
	}
	// WaitTailProb(t) = exp(-(mu-lambda)t) = exp(-t)
	if got, want := q.WaitTailProb(1), math.Exp(-1); !approx(got, want, tol) {
		t.Errorf("WaitTailProb(1) = %v, want %v", got, want)
	}
	// median sojourn time: -ln(0.5)/(mu-lambda)
	if got, want := q.WaitPercentile(0.5), math.Ln2; !approx(got, want, tol) {
		t.Errorf("WaitPercentile(0.5) = %v, want %v", got, want)
	}
	if _, err := NewMM1(3, 2); err != ErrUnstable {
		t.Errorf("expected ErrUnstable, got %v", err)
	}
	if _, err := NewMM1(-1, 2); err != ErrNonPositiveRate {
		t.Errorf("expected ErrNonPositiveRate, got %v", err)
	}
}

func TestMMc(t *testing.T) {
	q, err := NewMMc(2, 1, 3)
	if err != nil {
		t.Fatalf("NewMMc: %v", err)
	}
	tests := []struct {
		name string
		got  float64
		want float64
	}{
		{"OfferedLoad", q.OfferedLoad(), 2.0},
		{"Rho", q.Rho(), 2.0 / 3.0},
		{"P0", q.P0(), 1.0 / 9.0},
		{"ErlangC", q.ErlangC(), 4.0 / 9.0},
		{"Lq", q.Lq(), 8.0 / 9.0},
		{"L", q.L(), 2.0 + 8.0/9.0},
		{"Wq", q.Wq(), 4.0 / 9.0},
		{"W", q.W(), 1.0 + 4.0/9.0},
		{"MeanBusyServers", q.MeanBusyServers(), 2.0},
	}
	for _, tc := range tests {
		if !approx(tc.got, tc.want, 1e-9) {
			t.Errorf("MMc %s = %v, want %v", tc.name, tc.got, tc.want)
		}
	}
	// Single server M/M/c must reduce to M/M/1.
	q1, _ := NewMMc(1, 2, 1)
	m1, _ := NewMM1(1, 2)
	if !approx(q1.L(), m1.L(), 1e-9) {
		t.Errorf("MMc(c=1).L=%v, MM1.L=%v", q1.L(), m1.L())
	}
	// Distribution sums to ~1.
	sum := 0.0
	for n := 0; n < 200; n++ {
		sum += q.Pn(n)
	}
	if !approx(sum, 1, 1e-6) {
		t.Errorf("MMc distribution sums to %v", sum)
	}
	// Waiting-time tail at 0 equals delay probability.
	if !approx(q.WaitqTailProb(0), q.ErlangC(), 1e-9) {
		t.Errorf("WaitqTailProb(0)=%v, ErlangC=%v", q.WaitqTailProb(0), q.ErlangC())
	}
}

func TestMM1K(t *testing.T) {
	q, _ := NewMM1K(1, 2, 2) // rho=0.5, K=2
	if !approx(q.P0(), 0.5/0.875, tol) {
		t.Errorf("P0=%v", q.P0())
	}
	if !approx(q.BlockingProb(), 0.142857142857, 1e-9) {
		t.Errorf("BlockingProb=%v", q.BlockingProb())
	}
	if !approx(q.L(), 0.571428571428, 1e-9) {
		t.Errorf("L=%v", q.L())
	}
	// rho == 1 special case.
	qe, _ := NewMM1K(1, 1, 3)
	if !approx(qe.L(), 1.5, tol) {
		t.Errorf("L(rho=1)=%v, want 1.5", qe.L())
	}
	if !approx(qe.BlockingProb(), 0.25, tol) {
		t.Errorf("Blocking(rho=1)=%v, want 0.25", qe.BlockingProb())
	}
	// Distribution sums to one.
	s := 0.0
	for n := 0; n <= q.K; n++ {
		s += q.Pn(n)
	}
	if !approx(s, 1, tol) {
		t.Errorf("MM1K dist sum=%v", s)
	}
}

func TestMMcK(t *testing.T) {
	// M/M/1/K via M/M/c/K with c=1 must match closed form.
	a, _ := NewMMcK(1, 2, 1, 2)
	b, _ := NewMM1K(1, 2, 2)
	if !approx(a.BlockingProb(), b.BlockingProb(), 1e-9) {
		t.Errorf("MMcK block=%v, MM1K block=%v", a.BlockingProb(), b.BlockingProb())
	}
	if !approx(a.L(), b.L(), 1e-9) {
		t.Errorf("MMcK L=%v, MM1K L=%v", a.L(), b.L())
	}
	q, _ := NewMMcK(3, 2, 2, 5)
	s := 0.0
	for n := 0; n <= q.K; n++ {
		s += q.Pn(n)
	}
	if !approx(s, 1, 1e-9) {
		t.Errorf("MMcK dist sum=%v", s)
	}
	if q.BlockingProb() <= 0 || q.BlockingProb() >= 1 {
		t.Errorf("MMcK blocking out of range: %v", q.BlockingProb())
	}
	// Little's law consistency.
	if !approx(q.L(), q.EffectiveArrivalRate()*q.W(), 1e-9) {
		t.Errorf("Little's law fails: L=%v lambdaEff*W=%v", q.L(), q.EffectiveArrivalRate()*q.W())
	}
}

func TestMMInf(t *testing.T) {
	q, _ := NewMMInf(2, 1)
	if !approx(q.L(), 2, tol) {
		t.Errorf("L=%v", q.L())
	}
	if !approx(q.P0(), math.Exp(-2), tol) {
		t.Errorf("P0=%v", q.P0())
	}
	if !approx(q.Pn(2), math.Exp(-2)*2, tol) {
		t.Errorf("Pn(2)=%v", q.Pn(2))
	}
	if q.Wq() != 0 || q.Lq() != 0 {
		t.Errorf("infinite-server should have no wait")
	}
}

func TestMG1AndMD1(t *testing.T) {
	// M/G/1 with exponential service reduces to M/M/1.
	mg, _ := NewMG1(1, 0.5, 0.25) // var of exp(mu=2) is 0.25
	if !approx(mg.Lq(), 0.5, 1e-9) {
		t.Errorf("MG1 Lq=%v, want 0.5", mg.Lq())
	}
	if !approx(mg.ServiceSCV(), 1, 1e-9) {
		t.Errorf("MG1 SCV=%v, want 1", mg.ServiceSCV())
	}
	// From SCV constructor.
	mg2, _ := NewMG1FromSCV(1, 0.5, 1)
	if !approx(mg2.Lq(), mg.Lq(), 1e-12) {
		t.Errorf("NewMG1FromSCV mismatch")
	}
	// M/D/1.
	md, _ := NewMD1(1, 2)
	if !approx(md.Lq(), 0.25, 1e-9) {
		t.Errorf("MD1 Lq=%v, want 0.25", md.Lq())
	}
	if !approx(md.Wq(), 0.25, 1e-9) {
		t.Errorf("MD1 Wq=%v, want 0.25", md.Wq())
	}
	// PK free functions.
	if !approx(PollaczekKhinchineLq(0.5, 1), 0.5, 1e-9) {
		t.Errorf("PK Lq mismatch")
	}
	if !approx(PollaczekKhinchineLq(0.5, 0), 0.25, 1e-9) {
		t.Errorf("PK Lq (D) mismatch")
	}
}

func TestGG1(t *testing.T) {
	// Kingman with ca2=cs2=1 reproduces M/M/1 waiting time.
	g, _ := NewGG1(1, 0.5, 1, 1)
	if !approx(g.Wq(), 0.5, 1e-9) {
		t.Errorf("GG1 Wq=%v, want 0.5", g.Wq())
	}
	if !approx(KingmanWq(0.5, 1, 1, 0.5), 0.5, 1e-9) {
		t.Errorf("KingmanWq mismatch")
	}
	// Allen-Cunneen with ca2=cs2=1 reproduces M/M/c waiting time.
	mmc, _ := NewMMc(2, 1, 3)
	if !approx(AllenCunneenWq(2, 1, 3, 1, 1), mmc.Wq(), 1e-9) {
		t.Errorf("AllenCunneen != MMc Wq")
	}
}

func TestErlang(t *testing.T) {
	tests := []struct {
		c    int
		a    float64
		want float64
	}{
		{1, 1, 0.5},
		{2, 1, 0.2},
		{3, 2, 0.210526315789},
	}
	for _, tc := range tests {
		if got := ErlangB(tc.c, tc.a); !approx(got, tc.want, 1e-9) {
			t.Errorf("ErlangB(%d,%v)=%v, want %v", tc.c, tc.a, got, tc.want)
		}
	}
	if got := ErlangC(3, 2); !approx(got, 4.0/9.0, 1e-9) {
		t.Errorf("ErlangC(3,2)=%v, want %v", got, 4.0/9.0)
	}
	// Series consistency.
	ser := ErlangBSeries(3, 2)
	if !approx(ser[3], ErlangB(3, 2), 1e-12) {
		t.Errorf("ErlangBSeries mismatch")
	}
	// Server sizing.
	if n := ErlangBServersFor(2, 0.05, 20); n != 5 {
		t.Errorf("ErlangBServersFor=%d, want 5", n)
	}
	// Carried load.
	if !approx(CarriedLoad(3, 2), 2*(1-ErlangB(3, 2)), 1e-12) {
		t.Errorf("CarriedLoad mismatch")
	}
}

func TestBirthDeath(t *testing.T) {
	c, err := NewBirthDeath([]float64{1, 1, 1}, []float64{2, 2, 2})
	if err != nil {
		t.Fatalf("NewBirthDeath: %v", err)
	}
	p := c.StationaryDistribution()
	// M/M/1/3 with rho=0.5: unnormalized 1,0.5,0.25,0.125 -> total 1.875.
	if !approx(p[0], 1/1.875, 1e-9) {
		t.Errorf("p0=%v", p[0])
	}
	if !approx(c.MeanNumber(), 1.375/1.875, 1e-9) {
		t.Errorf("MeanNumber=%v", c.MeanNumber())
	}
	if !c.DetailedBalance(p, 1e-9) {
		t.Errorf("detailed balance should hold")
	}
	// Matches MM1K.
	q, _ := NewMM1K(1, 2, 3)
	if !approx(c.MeanNumber(), q.L(), 1e-9) {
		t.Errorf("BirthDeath L=%v, MM1K L=%v", c.MeanNumber(), q.L())
	}
}

func TestJackson(t *testing.T) {
	net, err := NewJacksonNetwork(
		[]float64{1, 0},
		[][]float64{{0, 1}, {0, 0}},
		[]float64{2, 2},
		nil,
	)
	if err != nil {
		t.Fatalf("NewJacksonNetwork: %v", err)
	}
	lam := net.ArrivalRates()
	if !approx(lam[0], 1, 1e-9) || !approx(lam[1], 1, 1e-9) {
		t.Errorf("arrival rates=%v", lam)
	}
	if !net.Stable() {
		t.Errorf("network should be stable")
	}
	if !approx(net.TotalMean(), 2, 1e-9) {
		t.Errorf("TotalMean=%v, want 2", net.TotalMean())
	}
	if !approx(net.Throughput(), 1, 1e-9) {
		t.Errorf("Throughput=%v, want 1", net.Throughput())
	}
	if !approx(net.SojournTime(), 2, 1e-9) {
		t.Errorf("SojournTime=%v, want 2", net.SojournTime())
	}
	// Feedback network: node routes back to itself with prob 0.5.
	net2, _ := NewJacksonNetwork(
		[]float64{1},
		[][]float64{{0.5}},
		[]float64{4},
		nil,
	)
	lam2 := net2.ArrivalRates()
	if !approx(lam2[0], 2, 1e-9) { // lambda = 1 + 0.5 lambda => lambda=2
		t.Errorf("feedback arrival=%v, want 2", lam2[0])
	}
}

func TestPriorityMM1(t *testing.T) {
	q, err := NewPriorityMM1([]float64{0.5, 0.5}, []float64{2, 2})
	if err != nil {
		t.Fatalf("NewPriorityMM1: %v", err)
	}
	if !approx(q.ClassWq(0), 1.0/3.0, 1e-9) {
		t.Errorf("ClassWq(0)=%v, want 1/3", q.ClassWq(0))
	}
	if !approx(q.ClassWq(1), 2.0/3.0, 1e-9) {
		t.Errorf("ClassWq(1)=%v, want 2/3", q.ClassWq(1))
	}
	// Conservation: aggregate weighted wait equals the M/M/1 wait of the
	// pooled stream (work-conserving, non-preemptive).
	pooled, _ := NewMM1(1, 2)
	if !approx(q.MeanWait(), pooled.Wq(), 1e-9) {
		t.Errorf("priority MeanWait=%v, pooled Wq=%v", q.MeanWait(), pooled.Wq())
	}
}

func TestLittle(t *testing.T) {
	if !approx(LittleL(2, 3), 6, tol) {
		t.Errorf("LittleL")
	}
	if !approx(LittleW(6, 2), 3, tol) {
		t.Errorf("LittleW")
	}
	if !approx(LittleLambda(6, 3), 2, tol) {
		t.Errorf("LittleLambda")
	}
}

func TestHelpers(t *testing.T) {
	if Factorial(5) != 120 {
		t.Errorf("Factorial(5)=%v", Factorial(5))
	}
	if !approx(PoissonPMF(0, 2), math.Exp(-2), 1e-12) {
		t.Errorf("PoissonPMF")
	}
	if !approx(PoissonCDF(1, 1), 2*math.Exp(-1), 1e-12) {
		t.Errorf("PoissonCDF")
	}
	if !approx(ErlangCDF(1, 2, 1), ExponentialCDF(2, 1), 1e-12) {
		t.Errorf("ErlangCDF(k=1) should equal exponential CDF")
	}
	if !approx(SquaredCoefficientOfVariation(2, 8), 2, 1e-12) {
		t.Errorf("SCV")
	}
}

func TestFlatAPI(t *testing.T) {
	if !approx(MM1L(1, 2), 1, tol) {
		t.Errorf("MM1L")
	}
	if !approx(MMcLq(2, 1, 3), 8.0/9.0, 1e-9) {
		t.Errorf("MMcLq")
	}
	if !approx(MD1Lq(1, 2), 0.25, 1e-9) {
		t.Errorf("MD1Lq")
	}
}

func ExampleMM1() {
	q, _ := NewMM1(8, 10) // arrivals 8/hr, service 10/hr
	fmt.Printf("rho=%.2f L=%.2f W=%.2f\n", q.Rho(), q.L(), q.W())
	// Output: rho=0.80 L=4.00 W=0.50
}

func ExampleErlangC() {
	// 10 erlangs of offered load, 14 trunks: probability of delay.
	fmt.Printf("%.4f\n", ErlangC(14, 10))
	// Output: 0.1741
}
