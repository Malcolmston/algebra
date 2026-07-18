package physics

import (
	"math"
	"testing"
)

// physicsClose reports whether a and b agree to a small relative/absolute
// tolerance, so floating-point rounding does not make the tests flaky.
func physicsClose(a, b, tol float64) bool {
	d := math.Abs(a - b)
	if d <= tol {
		return true
	}
	return d <= tol*math.Max(math.Abs(a), math.Abs(b))
}

func TestDimString(t *testing.T) {
	tests := []struct {
		name string
		dim  Dim
		want string
	}{
		{"length", DimLength, "m"},
		{"mass", DimMass, "kg"},
		{"time", DimTime, "s"},
		{"velocity", DimVelocity, "m·s⁻¹"},
		{"accel", DimAccel, "m·s⁻²"},
		{"force", DimForce, "kg·m·s⁻²"},
		{"energy", DimEnergy, "kg·m²·s⁻²"},
		{"power", DimPower, "kg·m²·s⁻³"},
		{"pressure", DimPressure, "kg·m⁻¹·s⁻²"},
		{"charge", DimCharge, "s·A"},
		{"voltage", DimVoltage, "kg·m²·s⁻³·A⁻¹"},
		{"dimensionless", Dim{}, "1"},
	}
	for _, tc := range tests {
		if got := tc.dim.String(); got != tc.want {
			t.Errorf("%s: String() = %q, want %q", tc.name, got, tc.want)
		}
	}
}

func TestQuantityMulDiv(t *testing.T) {
	m := NewQuantity(10, DimLength) // 10 m
	s := NewQuantity(2, DimTime)    // 2 s
	kg := NewQuantity(3, DimMass)   // 3 kg
	a := NewQuantity(5, DimAccel)   // 5 m/s²

	// velocity = length / time
	v := m.Div(s)
	if v.Dim != DimVelocity {
		t.Errorf("Div dim = %v, want %v", v.Dim, DimVelocity)
	}
	if !physicsClose(v.Value, 5, 1e-12) {
		t.Errorf("Div value = %v, want 5", v.Value)
	}

	// force = mass * accel
	f := kg.Mul(a)
	if f.Dim != DimForce {
		t.Errorf("Mul dim = %v, want %v", f.Dim, DimForce)
	}
	if !physicsClose(f.Value, 15, 1e-12) {
		t.Errorf("Mul value = %v, want 15", f.Value)
	}
}

func TestQuantityAddSub(t *testing.T) {
	a := NewQuantity(3, DimLength)
	b := NewQuantity(4, DimLength)
	sum, err := a.Add(b)
	if err != nil {
		t.Fatalf("Add returned error: %v", err)
	}
	if !physicsClose(sum.Value, 7, 1e-12) || sum.Dim != DimLength {
		t.Errorf("Add = %+v, want value 7 dim length", sum)
	}
	diff, err := b.Sub(a)
	if err != nil {
		t.Fatalf("Sub returned error: %v", err)
	}
	if !physicsClose(diff.Value, 1, 1e-12) {
		t.Errorf("Sub value = %v, want 1", diff.Value)
	}

	// Mismatched dimensions must error.
	if _, err := a.Add(NewQuantity(1, DimMass)); err == nil {
		t.Error("Add of length+mass: expected error, got nil")
	}
	if _, err := a.Sub(NewQuantity(1, DimTime)); err == nil {
		t.Error("Sub of length-time: expected error, got nil")
	}
}

func TestQuantityPow(t *testing.T) {
	m := NewQuantity(3, DimLength)
	area := m.Pow(2)
	wantArea := Dim{physicsIdxLength: 2}
	if area.Dim != wantArea {
		t.Errorf("Pow(2) dim = %v, want %v", area.Dim, wantArea)
	}
	if !physicsClose(area.Value, 9, 1e-12) {
		t.Errorf("Pow(2) value = %v, want 9", area.Value)
	}

	inv := m.Pow(-1)
	wantInv := Dim{physicsIdxLength: -1}
	if inv.Dim != wantInv {
		t.Errorf("Pow(-1) dim = %v, want %v", inv.Dim, wantInv)
	}
	if !physicsClose(inv.Value, 1.0/3.0, 1e-12) {
		t.Errorf("Pow(-1) value = %v, want %v", inv.Value, 1.0/3.0)
	}

	zero := m.Pow(0)
	if !zero.IsDimensionless() {
		t.Errorf("Pow(0) should be dimensionless, dim = %v", zero.Dim)
	}
	if !physicsClose(zero.Value, 1, 1e-12) {
		t.Errorf("Pow(0) value = %v, want 1", zero.Value)
	}
}

func TestIsDimensionless(t *testing.T) {
	if !NewQuantity(2, Dim{}).IsDimensionless() {
		t.Error("empty dim should be dimensionless")
	}
	if NewQuantity(2, DimLength).IsDimensionless() {
		t.Error("length should not be dimensionless")
	}
}

func TestParseQuantity(t *testing.T) {
	tests := []struct {
		value     float64
		symbol    string
		wantValue float64
		wantDim   Dim
	}{
		{1, "m", 1, DimLength},
		{1, "km", 1000, DimLength},
		{1, "mi", 1609.344, DimLength},
		{1, "yd", 0.9144, DimLength},
		{1, "lb", 0.45359237, DimMass},
		{1, "hr", 3600, DimTime},
		{1, "J", 1, DimEnergy},
		{1, "BTU", 1055.05585262, DimEnergy},
		{1, "therm", 1.05505585262e8, DimEnergy},
		{180, "deg", math.Pi, Dim{}},
		{1, "lbf", 4.4482216152605, DimForce},
		{1, "hp", 745.6998715822702, DimPower},
		{1, "psi", 6894.757293168361, DimPressure},
		{1, "atm", 101325, DimPressure},
		{1, "bar", 100000, DimPressure},
		{1, "mph", 0.44704, DimVelocity},
		{1, "knot", 0.5144444444444445, DimVelocity},
		{1, "gal", 0.003785411784, Dim{physicsIdxLength: 3}},
	}
	for _, tc := range tests {
		q, err := ParseQuantity(tc.value, tc.symbol)
		if err != nil {
			t.Errorf("ParseQuantity(%v, %q) error: %v", tc.value, tc.symbol, err)
			continue
		}
		if q.Dim != tc.wantDim {
			t.Errorf("ParseQuantity(%v, %q) dim = %v, want %v", tc.value, tc.symbol, q.Dim, tc.wantDim)
		}
		if !physicsClose(q.Value, tc.wantValue, 1e-9) {
			t.Errorf("ParseQuantity(%v, %q) value = %v, want %v", tc.value, tc.symbol, q.Value, tc.wantValue)
		}
	}

	if _, err := ParseQuantity(1, "nope"); err == nil {
		t.Error("ParseQuantity of unknown symbol: expected error")
	}
	// Affine temperature units are intentionally unsupported here.
	if _, err := ParseQuantity(1, "C"); err == nil {
		t.Error("ParseQuantity of Celsius: expected error")
	}
}

func TestQuantityTo(t *testing.T) {
	// 1609.344 m expressed in miles is exactly 1.
	q := NewQuantity(1609.344, DimLength)
	if got, err := q.To("mi"); err != nil || !physicsClose(got, 1, 1e-12) {
		t.Errorf("To(mi) = %v, err %v, want 1", got, err)
	}
	// 1000 m in km.
	if got, err := NewQuantity(1000, DimLength).To("km"); err != nil || !physicsClose(got, 1, 1e-12) {
		t.Errorf("To(km) = %v, err %v, want 1", got, err)
	}
	// A force expressed in lbf.
	f := NewQuantity(4.4482216152605, DimForce)
	if got, err := f.To("lbf"); err != nil || !physicsClose(got, 1, 1e-12) {
		t.Errorf("To(lbf) = %v, err %v, want 1", got, err)
	}
	// Round trip: parse then convert back.
	p, err := ParseQuantity(60, "mph")
	if err != nil {
		t.Fatalf("ParseQuantity(mph): %v", err)
	}
	if got, err := p.To("mph"); err != nil || !physicsClose(got, 60, 1e-12) {
		t.Errorf("round trip mph = %v, err %v, want 60", got, err)
	}

	// Dimension mismatch must error.
	if _, err := NewQuantity(1, DimLength).To("kg"); err == nil {
		t.Error("To(kg) on a length: expected dimension error")
	}
	// Unknown symbol must error.
	if _, err := NewQuantity(1, DimLength).To("nope"); err == nil {
		t.Error("To(nope): expected error")
	}
}

func TestComposedFromBase(t *testing.T) {
	// Build a force from base quantities and confirm it converts to lbf.
	kg := NewQuantity(0.45359237, DimMass)
	g := NewQuantity(StandardGravity, DimAccel)
	f := kg.Mul(g)
	if f.Dim != DimForce {
		t.Fatalf("mass*accel dim = %v, want force", f.Dim)
	}
	got, err := f.To("lbf")
	if err != nil {
		t.Fatalf("To(lbf): %v", err)
	}
	if !physicsClose(got, 1, 1e-12) {
		t.Errorf("1 lb * g in lbf = %v, want 1", got)
	}
}

func BenchmarkParseQuantity(b *testing.B) {
	b.ReportAllocs()
	var sink Quantity
	for i := 0; i < b.N; i++ {
		q, err := ParseQuantity(60, "mph")
		if err != nil {
			b.Fatal(err)
		}
		sink = q
	}
	_ = sink
}

func BenchmarkQuantityTo(b *testing.B) {
	q := NewQuantity(26.8224, DimVelocity)
	b.ReportAllocs()
	var sink float64
	for i := 0; i < b.N; i++ {
		v, err := q.To("mph")
		if err != nil {
			b.Fatal(err)
		}
		sink = v
	}
	_ = sink
}

func BenchmarkQuantityMul(b *testing.B) {
	kg := NewQuantity(3, DimMass)
	a := NewQuantity(5, DimAccel)
	b.ReportAllocs()
	var sink Quantity
	for i := 0; i < b.N; i++ {
		sink = kg.Mul(a)
	}
	_ = sink
}
