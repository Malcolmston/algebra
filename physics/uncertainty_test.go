package physics

import "testing"

func TestUncertaintyAddSub(t *testing.T) {
	a := NewMeasurement(10, 0.3)
	b := NewMeasurement(20, 0.4)

	sum := a.Add(b)
	assertApprox(t, "Add value", sum.Value, 30, 1e-12)
	assertApprox(t, "Add unc", sum.Uncertainty, 0.5, 1e-12) // √(0.09+0.16)

	diff := b.Sub(a)
	assertApprox(t, "Sub value", diff.Value, 10, 1e-12)
	assertApprox(t, "Sub unc", diff.Uncertainty, 0.5, 1e-12)

	// Negative uncertainty is normalized to its magnitude.
	if NewMeasurement(1, -0.2).Uncertainty != 0.2 {
		t.Error("uncertainty should be stored non-negative")
	}
}

func TestUncertaintyMulDiv(t *testing.T) {
	a := NewMeasurement(2, 0.1) // 5% relative
	b := NewMeasurement(3, 0.3) // 10% relative

	prod := a.Mul(b)
	assertApprox(t, "Mul value", prod.Value, 6, 1e-12)
	// rel = √(0.05² + 0.10²) = 0.111803; σ = 6·rel.
	assertApprox(t, "Mul unc", prod.Uncertainty, 6*0.11180339887498949, 1e-9)

	q := NewMeasurement(6, 0.6).Div(NewMeasurement(3, 0.3))
	assertApprox(t, "Div value", q.Value, 2, 1e-12)
	assertApprox(t, "Div unc", q.Uncertainty, 2*0.14142135623730951, 1e-9) // rel=√(0.1²+0.1²)
}

func TestUncertaintyScalePowRel(t *testing.T) {
	assertApprox(t, "RelativeUncertainty", NewMeasurement(10, 0.3).RelativeUncertainty(), 0.03, 1e-12)
	if NewMeasurement(0, 0.5).RelativeUncertainty() != 0 {
		t.Error("RelativeUncertainty of zero value should be 0")
	}

	s := NewMeasurement(2, 0.1).Scale(3)
	assertApprox(t, "Scale value", s.Value, 6, 1e-12)
	assertApprox(t, "Scale unc", s.Uncertainty, 0.3, 1e-12)

	// (2 ± 0.1)² = 4 ± 0.4  (relative uncertainty doubles).
	p := NewMeasurement(2, 0.1).Pow(2)
	assertApprox(t, "Pow value", p.Value, 4, 1e-12)
	assertApprox(t, "Pow unc", p.Uncertainty, 0.4, 1e-12)
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
		acc += x.Mul(y).Uncertainty
	}
	_ = acc
}
