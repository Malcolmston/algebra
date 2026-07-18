package geom3d

import (
	"math"
	"testing"
)

func TestMat3MulVecIdentity(t *testing.T) {
	v := NewVec3(1, 2, 3)
	if got := Mat3Identity().MulVec(v); !vecClose(got, v) {
		t.Errorf("I*v = %v", got)
	}
}

func TestMat3Mul(t *testing.T) {
	a := Mat3{{1, 2, 3}, {4, 5, 6}, {7, 8, 10}}
	inv, ok := a.Inverse()
	if !ok {
		t.Fatal("expected invertible")
	}
	if got := a.Mul(inv); !got.Equal(Mat3Identity(), 1e-9) {
		t.Errorf("A*A^-1 = %v", got)
	}
	if got := inv.Mul(a); !got.Equal(Mat3Identity(), 1e-9) {
		t.Errorf("A^-1*A = %v", got)
	}
}

func TestMat3Determinant(t *testing.T) {
	if got := Mat3Identity().Determinant(); math.Abs(got-1) > testEps {
		t.Errorf("det I = %v", got)
	}
	a := Mat3{{2, 0, 0}, {0, 3, 0}, {0, 0, 4}}
	if got := a.Determinant(); math.Abs(got-24) > testEps {
		t.Errorf("det diag = %v want 24", got)
	}
	// Singular matrix.
	s := Mat3{{1, 2, 3}, {2, 4, 6}, {7, 8, 9}}
	if _, ok := s.Inverse(); ok {
		t.Errorf("expected singular")
	}
}

func TestMat3TransposeTrace(t *testing.T) {
	a := Mat3{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}
	tr := a.Transpose()
	if tr[0][1] != a[1][0] || tr[2][0] != a[0][2] {
		t.Errorf("transpose wrong: %v", tr)
	}
	if got := a.Trace(); math.Abs(got-15) > testEps {
		t.Errorf("trace = %v want 15", got)
	}
}

func TestMat3FromRowsCols(t *testing.T) {
	r0, r1, r2 := Vec3{1, 2, 3}, Vec3{4, 5, 6}, Vec3{7, 8, 9}
	m := Mat3FromRows(r0, r1, r2)
	if !vecClose(m.Row(0), r0) || !vecClose(m.Row(2), r2) {
		t.Errorf("FromRows/Row wrong")
	}
	c := Mat3FromColumns(r0, r1, r2)
	if !vecClose(c.Col(0), r0) || !vecClose(c.Col(1), r1) {
		t.Errorf("FromColumns/Col wrong")
	}
}

func TestRotationAxes(t *testing.T) {
	// Rz(90deg) maps x -> y.
	if got := RotationZ(math.Pi / 2).MulVec(Vec3{1, 0, 0}); !got.Equal(Vec3{0, 1, 0}, 1e-9) {
		t.Errorf("Rz90 x = %v", got)
	}
	// Rx(90deg) maps y -> z.
	if got := RotationX(math.Pi / 2).MulVec(Vec3{0, 1, 0}); !got.Equal(Vec3{0, 0, 1}, 1e-9) {
		t.Errorf("Rx90 y = %v", got)
	}
	// Ry(90deg) maps z -> x.
	if got := RotationY(math.Pi / 2).MulVec(Vec3{0, 0, 1}); !got.Equal(Vec3{1, 0, 0}, 1e-9) {
		t.Errorf("Ry90 z = %v", got)
	}
	// Rotation preserves length.
	m := RotationX(0.7).Mul(RotationY(0.3)).Mul(RotationZ(-1.1))
	v := Vec3{1, 2, -3}
	if math.Abs(m.MulVec(v).Length()-v.Length()) > 1e-9 {
		t.Errorf("rotation changed length")
	}
	// Rotation matrix is orthogonal: R^T R = I.
	if !m.Transpose().Mul(m).Equal(Mat3Identity(), 1e-9) {
		t.Errorf("R not orthogonal")
	}
}

func TestAxisAngleAndRodrigues(t *testing.T) {
	axis := Vec3{0, 0, 1}
	ang := math.Pi / 2
	m := AxisAngleMatrix(axis, ang)
	// Should match RotationZ.
	if !m.Equal(RotationZ(ang), 1e-9) {
		t.Errorf("AxisAngle z != Rz: %v", m)
	}
	v := Vec3{1, 0, 0}
	// Rodrigues direct and matrix must agree.
	rd := Rodrigues(v, axis, ang)
	if !rd.Equal(m.MulVec(v), 1e-9) {
		t.Errorf("Rodrigues %v != matrix %v", rd, m.MulVec(v))
	}
	if !rd.Equal(Vec3{0, 1, 0}, 1e-9) {
		t.Errorf("Rodrigues result = %v want (0,1,0)", rd)
	}
	// Arbitrary axis, rotation preserves component along axis.
	ax := Vec3{1, 1, 1}
	rot := Rodrigues(Vec3{2, -1, 3}, ax, 1.234)
	u := ax.Unit()
	before := Vec3{2, -1, 3}.Dot(u)
	if math.Abs(rot.Dot(u)-before) > 1e-9 {
		t.Errorf("Rodrigues changed axial component")
	}
	// Zero axis returns input unchanged.
	if got := Rodrigues(v, Vec3{}, 1); !got.Equal(v, 1e-12) {
		t.Errorf("zero axis = %v", got)
	}
}

func TestEulerRoundTrip(t *testing.T) {
	cases := []EulerAngles{
		{Roll: 0, Pitch: 0, Yaw: 0},
		{Roll: 0.3, Pitch: -0.4, Yaw: 1.2},
		{Roll: -1.0, Pitch: 0.5, Yaw: -2.0},
		{Roll: 0.1, Pitch: 1.0, Yaw: 0.2},
	}
	for i, e := range cases {
		m := EulerToMatrix(e)
		got := MatrixToEuler(m)
		m2 := EulerToMatrix(got)
		// Compare via the matrices to sidestep angle wrapping.
		if !m.Equal(m2, 1e-9) {
			t.Errorf("case %d: euler roundtrip matrix mismatch\n%v\n%v", i, m, m2)
		}
	}
}

func TestEulerGimbalLock(t *testing.T) {
	e := EulerAngles{Roll: 0.5, Pitch: math.Pi / 2, Yaw: 0.3}
	m := EulerToMatrix(e)
	got := MatrixToEuler(m)
	if math.Abs(got.Pitch-math.Pi/2) > 1e-7 {
		t.Errorf("gimbal pitch = %v", got.Pitch)
	}
	// Reconstructed matrix must still match.
	if !EulerToMatrix(got).Equal(m, 1e-7) {
		t.Errorf("gimbal reconstruction failed")
	}
}
