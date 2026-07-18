package physics

import (
	"math"
	"testing"
)

func TestVec3Arithmetic(t *testing.T) {
	a := NewVec3(1, 2, 3)
	b := NewVec3(4, 5, 6)

	if got := a.Add(b); got != (Vec3{5, 7, 9}) {
		t.Errorf("Add = %+v", got)
	}
	if got := b.Sub(a); got != (Vec3{3, 3, 3}) {
		t.Errorf("Sub = %+v", got)
	}
	if got := a.Scale(2); got != (Vec3{2, 4, 6}) {
		t.Errorf("Scale = %+v", got)
	}
	if got := a.Neg(); got != (Vec3{-1, -2, -3}) {
		t.Errorf("Neg = %+v", got)
	}
	assertApprox(t, "Dot", a.Dot(b), 32, 1e-12) // 4+10+18
}

func TestVec3Cross(t *testing.T) {
	x := NewVec3(1, 0, 0)
	y := NewVec3(0, 1, 0)
	// x × y = z (right-handed).
	if got := x.Cross(y); got != (Vec3{0, 0, 1}) {
		t.Errorf("x×y = %+v, want (0,0,1)", got)
	}
	// Cross product is orthogonal to both operands.
	a := NewVec3(2, -3, 1)
	b := NewVec3(-1, 4, 5)
	c := a.Cross(b)
	assertApprox(t, "c·a", c.Dot(a), 0, 1e-12)
	assertApprox(t, "c·b", c.Dot(b), 0, 1e-12)
}

func TestVec3NormNormalize(t *testing.T) {
	v := NewVec3(3, 4, 0)
	assertApprox(t, "Norm", v.Norm(), 5, 1e-12)
	assertApprox(t, "Norm2", v.Norm2(), 25, 1e-12)

	u := v.Normalize()
	assertApprox(t, "unit len", u.Norm(), 1, 1e-12)
	assertApprox(t, "unit x", u.X, 0.6, 1e-12)

	// Zero vector normalizes to itself.
	if z := (Vec3{}).Normalize(); z != (Vec3{}) {
		t.Errorf("zero Normalize = %+v", z)
	}
}

func TestVec3DistanceLerpAngle(t *testing.T) {
	a := NewVec3(0, 0, 0)
	b := NewVec3(0, 3, 4)
	assertApprox(t, "Distance", a.Distance(b), 5, 1e-12)

	mid := a.Lerp(b, 0.5)
	if mid != (Vec3{0, 1.5, 2}) {
		t.Errorf("Lerp mid = %+v", mid)
	}

	// Orthogonal vectors → π/2.
	assertApprox(t, "Angle ortho", NewVec3(1, 0, 0).Angle(NewVec3(0, 1, 0)), math.Pi/2, 1e-12)
	// Parallel vectors → 0.
	assertApprox(t, "Angle parallel", NewVec3(1, 1, 1).Angle(NewVec3(2, 2, 2)), 0, 1e-9)
	// Anti-parallel vectors → π.
	assertApprox(t, "Angle anti", NewVec3(1, 0, 0).Angle(NewVec3(-1, 0, 0)), math.Pi, 1e-9)
}

func BenchmarkVec3Cross(b *testing.B) {
	u := NewVec3(1.5, -2.5, 3.5)
	v := NewVec3(-0.5, 4.0, 1.0)
	var acc Vec3
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		acc = acc.Add(u.Cross(v))
	}
	_ = acc
}

func BenchmarkVec3Normalize(b *testing.B) {
	v := NewVec3(1.5, -2.5, 3.5)
	var acc float64
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		acc += v.Normalize().X
	}
	_ = acc
}
