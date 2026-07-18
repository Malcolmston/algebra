package geom3d

import (
	"math"
	"testing"
)

const testEps = 1e-9

func vecClose(a, b Vec3) bool { return a.Equal(b, 1e-9) }

func TestVec3Arithmetic(t *testing.T) {
	a := NewVec3(1, 2, 3)
	b := NewVec3(4, -5, 6)
	if got := a.Add(b); !vecClose(got, Vec3{5, -3, 9}) {
		t.Errorf("Add = %v", got)
	}
	if got := a.Sub(b); !vecClose(got, Vec3{-3, 7, -3}) {
		t.Errorf("Sub = %v", got)
	}
	if got := a.Scale(2); !vecClose(got, Vec3{2, 4, 6}) {
		t.Errorf("Scale = %v", got)
	}
	if got := a.Neg(); !vecClose(got, Vec3{-1, -2, -3}) {
		t.Errorf("Neg = %v", got)
	}
	if got := b.Div(2); !vecClose(got, Vec3{2, -2.5, 3}) {
		t.Errorf("Div = %v", got)
	}
	if got := a.Hadamard(b); !vecClose(got, Vec3{4, -10, 18}) {
		t.Errorf("Hadamard = %v", got)
	}
}

func TestVec3DotCross(t *testing.T) {
	a := NewVec3(1, 2, 3)
	b := NewVec3(4, -5, 6)
	if got := a.Dot(b); math.Abs(got-12) > testEps {
		t.Errorf("Dot = %v want 12", got)
	}
	// x cross y = z
	if got := NewVec3(1, 0, 0).Cross(NewVec3(0, 1, 0)); !vecClose(got, Vec3{0, 0, 1}) {
		t.Errorf("Cross basis = %v", got)
	}
	if got := a.Cross(b); !vecClose(got, Vec3{27, 6, -13}) {
		t.Errorf("Cross = %v", got)
	}
	// a x a = 0
	if got := a.Cross(a); !got.IsZero() {
		t.Errorf("a x a = %v want 0", got)
	}
}

func TestVec3Norms(t *testing.T) {
	v := NewVec3(3, 4, 0)
	if got := v.Length(); math.Abs(got-5) > testEps {
		t.Errorf("Length = %v want 5", got)
	}
	if got := v.LengthSq(); math.Abs(got-25) > testEps {
		t.Errorf("LengthSq = %v want 25", got)
	}
	u, n := v.Normalize()
	if math.Abs(n-5) > testEps || !vecClose(u, Vec3{0.6, 0.8, 0}) {
		t.Errorf("Normalize = %v, %v", u, n)
	}
	if math.Abs(u.Length()-1) > testEps {
		t.Errorf("unit length = %v", u.Length())
	}
	if z, zn := (Vec3{}).Normalize(); !z.IsZero() || zn != 0 {
		t.Errorf("zero normalize = %v, %v", z, zn)
	}
}

func TestVec3Distance(t *testing.T) {
	a := NewVec3(1, 1, 1)
	b := NewVec3(1, 1, 4)
	if got := a.Distance(b); math.Abs(got-3) > testEps {
		t.Errorf("Distance = %v want 3", got)
	}
	if got := a.DistanceSq(b); math.Abs(got-9) > testEps {
		t.Errorf("DistanceSq = %v want 9", got)
	}
}

func TestVec3Angle(t *testing.T) {
	tests := []struct {
		a, b Vec3
		want float64
	}{
		{Vec3{1, 0, 0}, Vec3{0, 1, 0}, math.Pi / 2},
		{Vec3{1, 0, 0}, Vec3{1, 0, 0}, 0},
		{Vec3{1, 0, 0}, Vec3{-1, 0, 0}, math.Pi},
		{Vec3{1, 1, 0}, Vec3{1, 0, 0}, math.Pi / 4},
	}
	for i, tc := range tests {
		if got := tc.a.Angle(tc.b); math.Abs(got-tc.want) > 1e-9 {
			t.Errorf("case %d: Angle = %v want %v", i, got, tc.want)
		}
	}
}

func TestVec3Lerp(t *testing.T) {
	a := NewVec3(0, 0, 0)
	b := NewVec3(10, 20, 30)
	if got := a.Lerp(b, 0.5); !vecClose(got, Vec3{5, 10, 15}) {
		t.Errorf("Lerp .5 = %v", got)
	}
	if got := a.Lerp(b, 0); !vecClose(got, a) {
		t.Errorf("Lerp 0 = %v", got)
	}
	if got := a.Lerp(b, 1); !vecClose(got, b) {
		t.Errorf("Lerp 1 = %v", got)
	}
}

func TestVec3ProjectRejectReflect(t *testing.T) {
	v := NewVec3(2, 3, 0)
	onto := NewVec3(1, 0, 0)
	if got := v.Project(onto); !vecClose(got, Vec3{2, 0, 0}) {
		t.Errorf("Project = %v", got)
	}
	if got := v.Reject(onto); !vecClose(got, Vec3{0, 3, 0}) {
		t.Errorf("Reject = %v", got)
	}
	// project + reject reconstructs v
	if got := v.Project(onto).Add(v.Reject(onto)); !vecClose(got, v) {
		t.Errorf("proj+rej = %v", got)
	}
	// reflect (1,-1,0) across plane with unit normal (0,1,0) -> (1,1,0)
	if got := NewVec3(1, -1, 0).Reflect(NewVec3(0, 1, 0)); !vecClose(got, Vec3{1, 1, 0}) {
		t.Errorf("Reflect = %v", got)
	}
}

func TestVec3MinMaxAbs(t *testing.T) {
	a := NewVec3(1, -2, 3)
	b := NewVec3(-1, 2, 3)
	if got := a.Min(b); !vecClose(got, Vec3{-1, -2, 3}) {
		t.Errorf("Min = %v", got)
	}
	if got := a.Max(b); !vecClose(got, Vec3{1, 2, 3}) {
		t.Errorf("Max = %v", got)
	}
	if got := a.Abs(); !vecClose(got, Vec3{1, 2, 3}) {
		t.Errorf("Abs = %v", got)
	}
}

func TestScalarTriple(t *testing.T) {
	// Volume of unit cube basis = 1.
	if got := ScalarTriple(Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}); math.Abs(got-1) > testEps {
		t.Errorf("ScalarTriple basis = %v want 1", got)
	}
	// Coplanar vectors -> 0.
	if got := ScalarTriple(Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{1, 1, 0}); math.Abs(got) > testEps {
		t.Errorf("ScalarTriple coplanar = %v want 0", got)
	}
}
