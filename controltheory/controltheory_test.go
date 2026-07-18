package controltheory

import (
	"math"
	"math/cmplx"
	"sort"
	"testing"
)

func approx(a, b, tol float64) bool {
	return math.Abs(a-b) <= tol
}

func sortComplex(z []complex128) {
	sort.Slice(z, func(i, j int) bool {
		if real(z[i]) != real(z[j]) {
			return real(z[i]) < real(z[j])
		}
		return imag(z[i]) < imag(z[j])
	})
}

func TestPolyEvalAndMul(t *testing.T) {
	p := NewPoly(6, 5, 1) // 6 + 5s + s^2
	if got := p.Eval(2); !approx(got, 6+10+4, 1e-12) {
		t.Errorf("Eval(2)=%v want 20", got)
	}
	if p.Degree() != 2 {
		t.Errorf("Degree=%d want 2", p.Degree())
	}
	// (s+1)(s+2) = s^2+3s+2
	prod := NewPoly(1, 1).Mul(NewPoly(2, 1))
	want := []float64{2, 3, 1}
	for i, w := range want {
		if !approx(prod[i], w, 1e-12) {
			t.Errorf("Mul coeff %d = %v want %v", i, prod[i], w)
		}
	}
}

func TestPolyFromRootsAndRoots(t *testing.T) {
	p := PolyFromRoots(-2, -3) // s^2+5s+6
	want := []float64{6, 5, 1}
	for i, w := range want {
		if !approx(p[i], w, 1e-9) {
			t.Errorf("coeff %d = %v want %v", i, p[i], w)
		}
	}
	roots := NewPoly(6, 5, 1).Roots()
	sortComplex(roots)
	if len(roots) != 2 || !approx(real(roots[0]), -3, 1e-6) || !approx(real(roots[1]), -2, 1e-6) {
		t.Errorf("roots=%v want -3,-2", roots)
	}
	// Complex conjugate roots: s^2+2s+5 -> -1 +- 2i
	cr := NewPoly(5, 2, 1).Roots()
	sortComplex(cr)
	if !approx(real(cr[0]), -1, 1e-6) || !approx(imag(cr[0]), -2, 1e-6) {
		t.Errorf("complex roots=%v want -1-2i,-1+2i", cr)
	}
}

func TestPolyDivMod(t *testing.T) {
	// (s^2+5s+6)/(s+2) = s+3 remainder 0
	q, r := NewPoly(6, 5, 1).DivMod(NewPoly(2, 1))
	if !approx(q[0], 3, 1e-12) || !approx(q[1], 1, 1e-12) {
		t.Errorf("quotient=%v want {3,1}", q)
	}
	if r.Degree() != 0 || !approx(r.Eval(0), 0, 1e-12) {
		t.Errorf("remainder=%v want 0", r)
	}
}

func TestTransferFunctionBasics(t *testing.T) {
	g := NewTransferFunction([]float64{1}, []float64{2, 3, 1}) // 1/(s^2+3s+2)
	if !approx(g.DCGain(), 0.5, 1e-12) {
		t.Errorf("DCGain=%v want 0.5", g.DCGain())
	}
	if !g.IsStable() {
		t.Error("expected stable")
	}
	poles := g.Poles()
	sortComplex(poles)
	if !approx(real(poles[0]), -2, 1e-6) || !approx(real(poles[1]), -1, 1e-6) {
		t.Errorf("poles=%v want -2,-1", poles)
	}
}

func TestFrequencyResponse(t *testing.T) {
	g := NewTransferFunction([]float64{1}, []float64{1, 1}) // 1/(s+1)
	resp := g.FrequencyResponse(1)
	if !approx(real(resp), 0.5, 1e-12) || !approx(imag(resp), -0.5, 1e-12) {
		t.Errorf("G(j1)=%v want 0.5-0.5i", resp)
	}
	if !approx(MagnitudeDB(resp), -3.0103, 1e-3) {
		t.Errorf("magDB=%v want -3.0103", MagnitudeDB(resp))
	}
	if !approx(PhaseDeg(resp), -45, 1e-6) {
		t.Errorf("phase=%v want -45", PhaseDeg(resp))
	}
}

func TestSeriesParallelFeedback(t *testing.T) {
	g1 := NewTransferFunction([]float64{1}, []float64{1, 1}) // 1/(s+1)
	g2 := NewTransferFunction([]float64{1}, []float64{2, 1}) // 1/(s+2)
	ser := Series(g1, g2)                                    // 1/(s^2+3s+2)
	if !approx(ser.Den[0], 2, 1e-12) || !approx(ser.Den[1], 3, 1e-12) || !approx(ser.Den[2], 1, 1e-12) {
		t.Errorf("series den=%v want {2,3,1}", ser.Den)
	}
	par := Parallel(g1, g2) // (2s+3)/(s^2+3s+2)
	if !approx(par.Num[0], 3, 1e-12) || !approx(par.Num[1], 2, 1e-12) {
		t.Errorf("parallel num=%v want {3,2}", par.Num)
	}
	// Unity feedback of 1/s -> 1/(s+1)
	fb := UnityFeedback(NewTransferFunction([]float64{1}, []float64{0, 1}))
	if !approx(fb.Den[0], 1, 1e-12) || !approx(fb.Den[1], 1, 1e-12) {
		t.Errorf("feedback den=%v want {1,1}", fb.Den)
	}
}

func TestStateSpaceRoundTrip(t *testing.T) {
	g := NewTransferFunction([]float64{1}, []float64{2, 3, 1}) // 1/(s^2+3s+2)
	ss := TransferFunctionToStateSpace(g)
	if ss.Order() != 2 {
		t.Fatalf("order=%d want 2", ss.Order())
	}
	if !ss.IsControllable() || !ss.IsObservable() {
		t.Error("expected controllable and observable")
	}
	back := ss.TransferFunction()
	if !approx(back.Num[0], 1, 1e-9) {
		t.Errorf("roundtrip num=%v want {1}", back.Num)
	}
	wantDen := []float64{2, 3, 1}
	for i, w := range wantDen {
		if !approx(back.Den[i], w, 1e-9) {
			t.Errorf("roundtrip den[%d]=%v want %v", i, back.Den[i], w)
		}
	}
	cp := ss.CharacteristicPolynomial()
	for i, w := range wantDen {
		if !approx(cp[i], w, 1e-9) {
			t.Errorf("charpoly[%d]=%v want %v", i, cp[i], w)
		}
	}
}

func TestControllabilityObservabilityMatrices(t *testing.T) {
	// A = [[0,1],[-2,-3]], B=[0,1], C=[1,0]
	ss := NewStateSpace([][]float64{{0, 1}, {-2, -3}}, []float64{0, 1}, []float64{1, 0}, 0)
	if r := ss.ControllabilityRank(); r != 2 {
		t.Errorf("controllability rank=%d want 2", r)
	}
	if r := ss.ObservabilityRank(); r != 2 {
		t.Errorf("observability rank=%d want 2", r)
	}
	// Known uncontrollable system: diagonal A with B reaching only one mode.
	un := NewStateSpace([][]float64{{-1, 0}, {0, -2}}, []float64{1, 0}, []float64{1, 1}, 0)
	if un.IsControllable() {
		t.Error("expected uncontrollable")
	}
}

func TestStepResponse(t *testing.T) {
	g := NewTransferFunction([]float64{1}, []float64{1, 1}) // 1/(s+1): y=1-e^-t
	times := make([]float64, 201)
	for i := range times {
		times[i] = float64(i) * 0.01
	}
	y := g.StepResponse(times)
	if !approx(y[100], 1-math.Exp(-1), 1e-4) {
		t.Errorf("y(1)=%v want %v", y[100], 1-math.Exp(-1))
	}
	if !approx(y[200], 1-math.Exp(-2), 1e-4) {
		t.Errorf("y(2)=%v want %v", y[200], 1-math.Exp(-2))
	}
}

func TestImpulseResponse(t *testing.T) {
	g := NewTransferFunction([]float64{1}, []float64{1, 1}) // impulse resp = e^-t
	times := make([]float64, 101)
	for i := range times {
		times[i] = float64(i) * 0.01
	}
	y := g.ImpulseResponse(times)
	if !approx(y[0], 1, 1e-9) {
		t.Errorf("y(0)=%v want 1", y[0])
	}
	if !approx(y[100], math.Exp(-1), 1e-4) {
		t.Errorf("y(1)=%v want %v", y[100], math.Exp(-1))
	}
}

func TestRouthHurwitz(t *testing.T) {
	// (s+1)(s+2)(s+3) = s^3+6s^2+11s+6 stable
	stable := RouthHurwitz(NewPoly(6, 11, 6, 1))
	if !stable.Stable || stable.SignChanges != 0 {
		t.Errorf("expected stable, got %+v", stable)
	}
	// s^3+s^2+2s+8 has 2 RHP roots
	unst := RouthHurwitz(NewPoly(8, 2, 1, 1))
	if unst.Stable || unst.SignChanges != 2 {
		t.Errorf("expected 2 RHP roots, got %+v", unst)
	}
	if !IsHurwitzStable(NewPoly(6, 11, 6, 1)) {
		t.Error("IsHurwitzStable should be true")
	}
	if NumRightHalfPlaneRoots(NewPoly(8, 2, 1, 1)) != 2 {
		t.Error("NumRightHalfPlaneRoots should be 2")
	}
}

func TestBode(t *testing.T) {
	g := NewTransferFunction([]float64{1}, []float64{1, 1}) // 1/(s+1)
	pts := g.Bode([]float64{1})
	if !approx(pts[0].MagnitudeDB, -3.0103, 1e-3) {
		t.Errorf("mag=%v want -3.0103", pts[0].MagnitudeDB)
	}
	if !approx(pts[0].PhaseDeg, -45, 1e-6) {
		t.Errorf("phase=%v want -45", pts[0].PhaseDeg)
	}
}

func TestNyquist(t *testing.T) {
	g := NewTransferFunction([]float64{1}, []float64{1, 1})
	pts := g.Nyquist([]float64{1})
	if !approx(pts[0].Real, 0.5, 1e-9) || !approx(pts[0].Imag, -0.5, 1e-9) {
		t.Errorf("nyquist=%v want 0.5,-0.5", pts[0])
	}
}

func TestPhaseMargin(t *testing.T) {
	// G = 1/(s(s+1)); expected phase margin about 51.83 deg
	g := NewTransferFunction([]float64{1}, []float64{0, 1, 1})
	omegas := LogSpace(-2, 2, 4000)
	pm, ok := g.PhaseMargin(omegas)
	if !ok {
		t.Fatal("expected gain crossover")
	}
	if !approx(pm, 51.83, 0.5) {
		t.Errorf("phase margin=%v want ~51.83", pm)
	}
}

func TestGainMargin(t *testing.T) {
	// G = 1/(s+1)^3 -> gain margin 20*log10(8) = 18.06 dB
	den := NewPoly(1, 1).Mul(NewPoly(1, 1)).Mul(NewPoly(1, 1)) // (s+1)^3
	g := TransferFunction{Num: Poly{1}, Den: den}
	omegas := LogSpace(-2, 2, 4000)
	gm, ok := g.GainMargin(omegas)
	if !ok {
		t.Fatal("expected phase crossover")
	}
	if !approx(gm, 18.06, 0.2) {
		t.Errorf("gain margin=%v want ~18.06", gm)
	}
}

func TestSecondOrderSystem(t *testing.T) {
	g := SecondOrderSystem(2, 0.5)
	if !approx(g.Num[0], 4, 1e-12) || !approx(g.Den[0], 4, 1e-12) || !approx(g.Den[1], 2, 1e-12) {
		t.Errorf("2nd order tf wrong: %+v", g)
	}
	poles := g.Poles()
	if !approx(NaturalFrequency(poles[0]), 2, 1e-6) {
		t.Errorf("wn=%v want 2", NaturalFrequency(poles[0]))
	}
	if !approx(DampingRatio(poles[0]), 0.5, 1e-6) {
		t.Errorf("zeta=%v want 0.5", DampingRatio(poles[0]))
	}
}

func TestPIDUpdate(t *testing.T) {
	c := NewPIDController(2, 0, 0)
	if got := c.Update(3, 1); !approx(got, 6, 1e-12) {
		t.Errorf("P output=%v want 6", got)
	}
	ci := NewPIDController(0, 1, 0)
	_ = ci.Update(1, 1)
	if got := ci.Update(1, 1); !approx(got, 2, 1e-12) {
		t.Errorf("I output=%v want 2", got)
	}
	if !approx(ci.Integral(), 2, 1e-12) {
		t.Errorf("integral=%v want 2", ci.Integral())
	}
	tf := NewPIDController(1, 2, 3).TransferFunction()
	if !approx(tf.Num[0], 2, 1e-12) || !approx(tf.Num[1], 1, 1e-12) || !approx(tf.Num[2], 3, 1e-12) {
		t.Errorf("PID tf num=%v want {2,1,3}", tf.Num)
	}
}

func TestZieglerNichols(t *testing.T) {
	c := ZieglerNicholsPID(2, 4)
	if !approx(c.Kp, 1.2, 1e-9) || !approx(c.Ki, 0.6, 1e-9) || !approx(c.Kd, 0.6, 1e-9) {
		t.Errorf("ZN gains=%+v want Kp1.2 Ki0.6 Kd0.6", c)
	}
}

func TestEvaluatePoleInfinite(t *testing.T) {
	g := NewTransferFunction([]float64{1}, []float64{1, 1})
	if !cmplx.IsInf(g.Evaluate(complex(-1, 0))) {
		t.Error("expected infinite response at pole")
	}
}

func BenchmarkStepResponse(b *testing.B) {
	g := SecondOrderSystem(2, 0.5)
	times := make([]float64, 1000)
	for i := range times {
		times[i] = float64(i) * 0.01
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = g.StepResponse(times)
	}
}
