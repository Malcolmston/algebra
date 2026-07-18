package physics

import (
	"math"
	"testing"
)

// physicsVecTestClose reports whether every component of a and b differs by at
// most tol. It is named distinctively to avoid colliding with helpers defined
// in sibling test files of the physics package.
func physicsVecTestClose(a, b Vec3, tol float64) bool {
	return math.Abs(a.X-b.X) <= tol &&
		math.Abs(a.Y-b.Y) <= tol &&
		math.Abs(a.Z-b.Z) <= tol
}

// physicsScalarTestClose reports whether a and b differ by at most tol.
func physicsScalarTestClose(a, b, tol float64) bool {
	return math.Abs(a-b) <= tol
}

func TestUnit(t *testing.T) {
	tests := []struct {
		name string
		in   Vec3
		want Vec3
	}{
		{"axis-x", Vec3{3, 0, 0}, Vec3{1, 0, 0}},
		{"3-4-5", Vec3{0, 3, 4}, Vec3{0, 0.6, 0.8}},
		{"negative", Vec3{-5, 0, 0}, Vec3{-1, 0, 0}},
		{"zero", Vec3{0, 0, 0}, Vec3{0, 0, 0}},
	}
	for _, tc := range tests {
		got := tc.in.Unit()
		if !physicsVecTestClose(got, tc.want, 1e-12) {
			t.Errorf("%s: Unit() = %v, want %v", tc.name, got, tc.want)
		}
		if n := tc.in.Norm(); n != 0 {
			if l := got.Norm(); !physicsScalarTestClose(l, 1, 1e-12) {
				t.Errorf("%s: |Unit()| = %v, want 1", tc.name, l)
			}
		}
	}
}

func TestAngleBetween(t *testing.T) {
	tests := []struct {
		name string
		a, b Vec3
		want float64
	}{
		{"orthogonal", Vec3{1, 0, 0}, Vec3{0, 1, 0}, math.Pi / 2},
		{"parallel", Vec3{2, 0, 0}, Vec3{5, 0, 0}, 0},
		{"antiparallel", Vec3{1, 0, 0}, Vec3{-1, 0, 0}, math.Pi},
		{"45deg", Vec3{1, 0, 0}, Vec3{1, 1, 0}, math.Pi / 4},
		{"zero-operand", Vec3{0, 0, 0}, Vec3{1, 1, 1}, 0},
	}
	for _, tc := range tests {
		got := tc.a.AngleBetween(tc.b)
		if !physicsScalarTestClose(got, tc.want, 1e-12) {
			t.Errorf("%s: AngleBetween() = %v, want %v", tc.name, got, tc.want)
		}
	}
}

func TestAddInto(t *testing.T) {
	var dst Vec3
	AddInto(&dst, Vec3{1, 2, 3}, Vec3{4, 5, 6})
	if want := (Vec3{5, 7, 9}); !physicsVecTestClose(dst, want, 0) {
		t.Errorf("AddInto = %v, want %v", dst, want)
	}
	// Aliasing: dst == a is allowed.
	dst = Vec3{1, 1, 1}
	AddInto(&dst, dst, Vec3{2, 3, 4})
	if want := (Vec3{3, 4, 5}); !physicsVecTestClose(dst, want, 0) {
		t.Errorf("AddInto aliased = %v, want %v", dst, want)
	}
}

func TestScaleInto(t *testing.T) {
	var dst Vec3
	ScaleInto(&dst, Vec3{1, -2, 3}, 2)
	if want := (Vec3{2, -4, 6}); !physicsVecTestClose(dst, want, 0) {
		t.Errorf("ScaleInto = %v, want %v", dst, want)
	}
	dst = Vec3{2, 4, 6}
	ScaleInto(&dst, dst, 0.5)
	if want := (Vec3{1, 2, 3}); !physicsVecTestClose(dst, want, 0) {
		t.Errorf("ScaleInto aliased = %v, want %v", dst, want)
	}
}

func TestSumSlice(t *testing.T) {
	tests := []struct {
		name string
		in   []Vec3
		want Vec3
	}{
		{"nil", nil, Vec3{0, 0, 0}},
		{"empty", []Vec3{}, Vec3{0, 0, 0}},
		{"single", []Vec3{{1, 2, 3}}, Vec3{1, 2, 3}},
		{"several", []Vec3{{1, 2, 3}, {4, 5, 6}, {-5, -7, -9}}, Vec3{0, 0, 0}},
	}
	for _, tc := range tests {
		got := SumSlice(tc.in)
		if !physicsVecTestClose(got, tc.want, 0) {
			t.Errorf("%s: SumSlice = %v, want %v", tc.name, got, tc.want)
		}
	}
}

func TestGradient(t *testing.T) {
	// f(x,y,z) = x^2 + y^2 + z^2 has gradient (2x, 2y, 2z).
	f := func(v Vec3) float64 { return v.X*v.X + v.Y*v.Y + v.Z*v.Z }
	p := Vec3{1, 2, 3}
	got := Gradient(f, p, 1e-4)
	want := Vec3{2, 4, 6}
	if !physicsVecTestClose(got, want, 1e-6) {
		t.Errorf("Gradient = %v, want %v", got, want)
	}
	// Non-positive step must fall back to the default and still be accurate.
	got = Gradient(f, p, 0)
	if !physicsVecTestClose(got, want, 1e-3) {
		t.Errorf("Gradient (default step) = %v, want %v", got, want)
	}
}

func TestDivergence(t *testing.T) {
	// F(x,y,z) = (x, y, z) has divergence 3 everywhere.
	f := func(v Vec3) Vec3 { return v }
	got := Divergence(f, Vec3{1, 2, 3}, 1e-4)
	if !physicsScalarTestClose(got, 3, 1e-6) {
		t.Errorf("Divergence = %v, want 3", got)
	}
	// F(x,y,z) = (x^2, y^2, z^2) has divergence 2x+2y+2z.
	g := func(v Vec3) Vec3 { return Vec3{v.X * v.X, v.Y * v.Y, v.Z * v.Z} }
	p := Vec3{1, 2, 3}
	got = Divergence(g, p, 1e-4)
	want := 2*p.X + 2*p.Y + 2*p.Z
	if !physicsScalarTestClose(got, want, 1e-5) {
		t.Errorf("Divergence quadratic = %v, want %v", got, want)
	}
}

func TestCurl(t *testing.T) {
	// Solid-body rotation F = (-y, x, 0) has curl (0, 0, 2).
	f := func(v Vec3) Vec3 { return Vec3{-v.Y, v.X, 0} }
	got := Curl(f, Vec3{3, -4, 5}, 1e-4)
	want := Vec3{0, 0, 2}
	if !physicsVecTestClose(got, want, 1e-6) {
		t.Errorf("Curl rotation = %v, want %v", got, want)
	}
	// A gradient field F = (x, y, z) is irrotational: curl is zero.
	g := func(v Vec3) Vec3 { return v }
	got = Curl(g, Vec3{1, 2, 3}, 1e-4)
	if !physicsVecTestClose(got, Vec3{0, 0, 0}, 1e-6) {
		t.Errorf("Curl irrotational = %v, want zero", got)
	}
}

func TestLaplacian(t *testing.T) {
	// f = x^2 + y^2 + z^2 has Laplacian 6 everywhere.
	f := func(v Vec3) float64 { return v.X*v.X + v.Y*v.Y + v.Z*v.Z }
	got := Laplacian(f, Vec3{1, 2, 3}, 1e-3)
	if !physicsScalarTestClose(got, 6, 1e-4) {
		t.Errorf("Laplacian = %v, want 6", got)
	}
	// A linear field has zero Laplacian.
	g := func(v Vec3) float64 { return 2*v.X - 3*v.Y + v.Z }
	got = Laplacian(g, Vec3{5, 6, 7}, 1e-3)
	if !physicsScalarTestClose(got, 0, 1e-4) {
		t.Errorf("Laplacian linear = %v, want 0", got)
	}
}

func TestStepFallback(t *testing.T) {
	if got := physicsStep(-1); got != physicsDefaultStep {
		t.Errorf("physicsStep(-1) = %v, want %v", got, physicsDefaultStep)
	}
	if got := physicsStep(0); got != physicsDefaultStep {
		t.Errorf("physicsStep(0) = %v, want %v", got, physicsDefaultStep)
	}
	if got := physicsStep(0.25); got != 0.25 {
		t.Errorf("physicsStep(0.25) = %v, want 0.25", got)
	}
}

func BenchmarkAddInto(b *testing.B) {
	var dst Vec3
	x, y := Vec3{1, 2, 3}, Vec3{4, 5, 6}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		AddInto(&dst, x, y)
	}
	_ = dst
}

func BenchmarkScaleInto(b *testing.B) {
	var dst Vec3
	x := Vec3{1, 2, 3}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		ScaleInto(&dst, x, 2)
	}
	_ = dst
}

func BenchmarkSumSlice(b *testing.B) {
	v := make([]Vec3, 1024)
	for i := range v {
		v[i] = Vec3{float64(i), float64(2 * i), float64(3 * i)}
	}
	b.ReportAllocs()
	var sink Vec3
	for i := 0; i < b.N; i++ {
		sink = SumSlice(v)
	}
	_ = sink
}

func BenchmarkGradient(b *testing.B) {
	f := func(v Vec3) float64 { return v.X*v.X + v.Y*v.Y + v.Z*v.Z }
	p := Vec3{1, 2, 3}
	b.ReportAllocs()
	var sink Vec3
	for i := 0; i < b.N; i++ {
		sink = Gradient(f, p, 1e-4)
	}
	_ = sink
}
