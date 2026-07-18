package diffgeo

import (
	"math"
	"testing"
)

const testTol = 1e-9

func TestVec3Arithmetic(t *testing.T) {
	a := NewVec3(1, 2, 3)
	b := NewVec3(4, -5, 6)
	if got := a.Add(b); !got.Equal(Vec3{5, -3, 9}, testTol) {
		t.Errorf("Add = %v", got)
	}
	if got := a.Sub(b); !got.Equal(Vec3{-3, 7, -3}, testTol) {
		t.Errorf("Sub = %v", got)
	}
	if got := a.Scale(2); !got.Equal(Vec3{2, 4, 6}, testTol) {
		t.Errorf("Scale = %v", got)
	}
	if got := a.Neg(); !got.Equal(Vec3{-1, -2, -3}, testTol) {
		t.Errorf("Neg = %v", got)
	}
	if got := a.Lerp(b, 0.5); !got.Equal(Vec3{2.5, -1.5, 4.5}, testTol) {
		t.Errorf("Lerp = %v", got)
	}
}

func TestVec3Products(t *testing.T) {
	a := NewVec3(1, 2, 3)
	b := NewVec3(4, -5, 6)
	if got := a.Dot(b); math.Abs(got-12) > testTol { // 4 -10 +18
		t.Errorf("Dot = %v want 12", got)
	}
	// Standard basis cross product.
	x := Vec3{1, 0, 0}
	y := Vec3{0, 1, 0}
	if got := x.Cross(y); !got.Equal(Vec3{0, 0, 1}, testTol) {
		t.Errorf("x cross y = %v want z", got)
	}
	// Cross is orthogonal to both operands.
	c := a.Cross(b)
	if math.Abs(c.Dot(a)) > testTol || math.Abs(c.Dot(b)) > testTol {
		t.Errorf("cross not orthogonal: %v", c)
	}
}

func TestVec3NormAndAngle(t *testing.T) {
	v := Vec3{3, 4, 0}
	if got := v.Norm(); math.Abs(got-5) > testTol {
		t.Errorf("Norm = %v want 5", got)
	}
	if got := v.Norm2(); math.Abs(got-25) > testTol {
		t.Errorf("Norm2 = %v want 25", got)
	}
	if got := v.Normalize().Norm(); math.Abs(got-1) > testTol {
		t.Errorf("normalized norm = %v want 1", got)
	}
	if got := (Vec3{}).Normalize(); !got.IsZero(testTol) {
		t.Errorf("zero normalize = %v", got)
	}
	// Right angle between axes.
	if got := (Vec3{1, 0, 0}).Angle(Vec3{0, 1, 0}); math.Abs(got-math.Pi/2) > testTol {
		t.Errorf("angle = %v want pi/2", got)
	}
	// 45 degrees.
	if got := (Vec3{1, 0, 0}).Angle(Vec3{1, 1, 0}); math.Abs(got-math.Pi/4) > 1e-9 {
		t.Errorf("angle = %v want pi/4", got)
	}
	if got := (Vec3{1, 2, 2}).Distance(Vec3{1, 2, 2}); got != 0 {
		t.Errorf("self distance = %v", got)
	}
}
