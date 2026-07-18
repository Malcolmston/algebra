package physics

import (
	"math"
	"testing"
)

func TestMeasurementConstructAndRelative(t *testing.T) {
	// Negative sigma is normalized to its magnitude.
	if m := NewMeasurement(1, -0.2); m.Sigma != 0.2 {
		t.Errorf("NewMeasurement sigma = %v, want 0.2", m.Sigma)
	}
	assertApprox(t, "Relative", NewMeasurement(10, 0.3).Relative(), 0.03, 1e-12)
}

func TestMeasurementAddSub(t *testing.T) {
	a := NewMeasurement(10, 0.3)
	b := NewMeasurement(20, 0.4)

	sum := a.Add(b)
	assertApprox(t, "Add value", sum.Value, 30, 1e-12)
	assertApprox(t, "Add sigma", sum.Sigma, 0.5, 1e-12) // √(0.09+0.16)

	diff := b.Sub(a)
	assertApprox(t, "Sub value", diff.Value, 10, 1e-12)
	assertApprox(t, "Sub sigma", diff.Sigma, 0.5, 1e-12)
}

func TestMeasurementMulDiv(t *testing.T) {
	// (2 ± 0.1)·(3 ± 0.3): rel = √(0.05² + 0.10²) = 0.111803…, σ = 6·rel.
	prod := NewMeasurement(2, 0.1).Mul(NewMeasurement(3, 0.3))
	assertApprox(t, "Mul value", prod.Value, 6, 1e-12)
	assertApprox(t, "Mul sigma", prod.Sigma, 6*0.11180339887498949, 1e-12)

	// (6 ± 0.6)/(3 ± 0.3): rel = √(0.1² + 0.1²) = 0.141421…, σ = 2·rel.
	q := NewMeasurement(6, 0.6).Div(NewMeasurement(3, 0.3))
	assertApprox(t, "Div value", q.Value, 2, 1e-12)
	assertApprox(t, "Div sigma", q.Sigma, 2*0.14142135623730951, 1e-12)
}

func TestMeasurementScalePow(t *testing.T) {
	s := NewMeasurement(2, 0.1).Scale(-3)
	assertApprox(t, "Scale value", s.Value, -6, 1e-12)
	assertApprox(t, "Scale sigma", s.Sigma, 0.3, 1e-12) // |k| scales sigma

	// (2 ± 0.1)² = 4 ± 0.4 (relative uncertainty doubles).
	p := NewMeasurement(2, 0.1).Pow(2)
	assertApprox(t, "Pow value", p.Value, 4, 1e-12)
	assertApprox(t, "Pow sigma", p.Sigma, 0.4, 1e-12)
}

func TestMeasurementApply(t *testing.T) {
	cases := []struct {
		name      string
		f         func(float64) float64
		in        Measurement
		wantValue float64
		wantSigma float64
	}{
		// f(x)=x², f'(x)=2x. At x=3: value 9, σ=|6|·0.1=0.6.
		{"square", func(x float64) float64 { return x * x }, NewMeasurement(3, 0.1), 9, 0.6},
		// f(x)=eˣ, f'(x)=eˣ. At x=1: value e, σ=e·0.1.
		{"exp", math.Exp, NewMeasurement(1, 0.1), math.E, math.E * 0.1},
		// f(x)=ln x, f'(x)=1/x. At x=2: value ln2, σ=(1/2)·0.2=0.1.
		{"log", math.Log, NewMeasurement(2, 0.2), math.Ln2, 0.1},
	}
	for _, c := range cases {
		got := c.in.Apply(c.f)
		assertApprox(t, c.name+" value", got.Value, c.wantValue, 1e-12)
		assertApprox(t, c.name+" sigma", got.Sigma, c.wantSigma, 1e-6)
	}
}

func TestPropagate(t *testing.T) {
	// f(x,y)=x·y at (2,3): ∂f/∂x=y=3, ∂f/∂y=x=2.
	// σ = √((3·0.1)² + (2·0.3)²) = √0.45 = 0.670820…
	m, err := Propagate(
		func(v []float64) float64 { return v[0] * v[1] },
		[]float64{2, 3}, []float64{0.1, 0.3},
	)
	if err != nil {
		t.Fatalf("Propagate returned error: %v", err)
	}
	assertApprox(t, "Propagate value", m.Value, 6, 1e-9)
	assertApprox(t, "Propagate sigma", m.Sigma, math.Sqrt(0.45), 1e-6)

	// f(x,y)=x+y: σ = √(σx²+σy²), independent of the point.
	sumM, err := Propagate(
		func(v []float64) float64 { return v[0] + v[1] },
		[]float64{10, 20}, []float64{0.3, 0.4},
	)
	if err != nil {
		t.Fatalf("Propagate(sum) returned error: %v", err)
	}
	assertApprox(t, "Propagate sum value", sumM.Value, 30, 1e-9)
	assertApprox(t, "Propagate sum sigma", sumM.Sigma, 0.5, 1e-6)
}

func TestPropagateLengthMismatch(t *testing.T) {
	_, err := Propagate(
		func(v []float64) float64 { return v[0] },
		[]float64{1, 2}, []float64{0.1},
	)
	if err == nil {
		t.Fatal("Propagate: expected error on length mismatch, got nil")
	}
}

func TestMeasurementString(t *testing.T) {
	if got := NewMeasurement(3, 0.5).String(); got != "3 ± 0.5" {
		t.Errorf("String = %q, want %q", got, "3 ± 0.5")
	}
}

func BenchmarkMeasurementMul(b *testing.B) {
	x := NewMeasurement(2, 0.1)
	y := NewMeasurement(3, 0.3)
	var acc float64
	for i := 0; i < b.N; i++ {
		acc += x.Mul(y).Sigma
	}
	_ = acc
}

func BenchmarkPropagate(b *testing.B) {
	f := func(v []float64) float64 { return v[0] * v[1] * v[2] }
	means := []float64{2, 3, 5}
	sigmas := []float64{0.1, 0.2, 0.3}
	var acc float64
	for i := 0; i < b.N; i++ {
		m, _ := Propagate(f, means, sigmas)
		acc += m.Sigma
	}
	_ = acc
}
