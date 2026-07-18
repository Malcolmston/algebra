package hypercomplex

import (
	"math"
	"testing"
)

const qTol = 1e-9

func TestQuaternionMulKnownAnswers(t *testing.T) {
	i := Quat(0, 1, 0, 0)
	j := Quat(0, 0, 1, 0)
	k := Quat(0, 0, 0, 1)
	// i*j = k, j*k = i, k*i = j.
	if got := i.Mul(j); !got.Equal(k, qTol) {
		t.Errorf("i*j = %+v want %+v", got, k)
	}
	if got := j.Mul(k); !got.Equal(i, qTol) {
		t.Errorf("j*k = %+v want %+v", got, i)
	}
	if got := k.Mul(i); !got.Equal(j, qTol) {
		t.Errorf("k*i = %+v want %+v", got, j)
	}
	// i*i = j*j = k*k = -1.
	neg1 := Quat(-1, 0, 0, 0)
	if got := i.Mul(i); !got.Equal(neg1, qTol) {
		t.Errorf("i*i = %+v want -1", got)
	}
	// Non-commutative: j*i = -k.
	if got := j.Mul(i); !got.Equal(k.Neg(), qTol) {
		t.Errorf("j*i = %+v want -k", got)
	}
}

func TestQuaternionConjNormInverse(t *testing.T) {
	q := Quat(1, 2, 3, 4)
	if got := q.NormSq(); math.Abs(got-30) > qTol {
		t.Errorf("normSq = %v want 30", got)
	}
	if got := q.Norm(); math.Abs(got-math.Sqrt(30)) > qTol {
		t.Errorf("norm = %v", got)
	}
	// q * q^-1 = 1.
	if got := q.Mul(q.Inverse()); !got.Equal(IdentityQuat(), qTol) {
		t.Errorf("q*q^-1 = %+v want identity", got)
	}
	// conj(q)*conj(q) relation: (q*)* = q.
	if got := q.Conj().Conj(); !got.Equal(q, qTol) {
		t.Errorf("(q*)* = %+v want %+v", got, q)
	}
}

func TestQuaternionRotateVector(t *testing.T) {
	// 90 degrees about Z sends x-axis to y-axis.
	q := QuatFromAxisAngle(V3(0, 0, 1), math.Pi/2)
	got := q.RotateVector(V3(1, 0, 0))
	if !got.Equal(V3(0, 1, 0), 1e-9) {
		t.Errorf("rotate x by 90 about z = %+v want (0,1,0)", got)
	}
	// 120 degrees about (1,1,1) cyclically permutes the axes: x -> y.
	q2 := QuatFromAxisAngle(V3(1, 1, 1), 2*math.Pi/3)
	got2 := q2.RotateVector(V3(1, 0, 0))
	if !got2.Equal(V3(0, 1, 0), 1e-9) {
		t.Errorf("rotate x by 120 about diagonal = %+v want (0,1,0)", got2)
	}
}

func TestQuaternionAxisAngleRoundTrip(t *testing.T) {
	axis := V3(1, 2, 3).Normalize()
	angle := 1.234
	q := QuatFromAxisAngle(axis, angle)
	gotAxis, gotAngle := q.ToAxisAngle()
	if math.Abs(gotAngle-angle) > 1e-9 {
		t.Errorf("angle = %v want %v", gotAngle, angle)
	}
	if !gotAxis.Equal(axis, 1e-9) {
		t.Errorf("axis = %+v want %+v", gotAxis, axis)
	}
}

func TestQuaternionEulerRoundTrip(t *testing.T) {
	cases := []struct{ roll, pitch, yaw float64 }{
		{0.1, 0.2, 0.3},
		{-0.5, 0.4, 1.1},
		{0.9, -0.7, -0.2},
	}
	for _, c := range cases {
		q := QuatFromEuler(c.roll, c.pitch, c.yaw)
		r, p, y := q.ToEuler()
		if math.Abs(r-c.roll) > 1e-9 || math.Abs(p-c.pitch) > 1e-9 || math.Abs(y-c.yaw) > 1e-9 {
			t.Errorf("euler round trip: got (%v,%v,%v) want (%v,%v,%v)", r, p, y, c.roll, c.pitch, c.yaw)
		}
	}
}

func TestQuaternionRotationMatrixRoundTrip(t *testing.T) {
	q := QuatFromEuler(0.3, -0.6, 1.2).Normalize()
	m := q.ToRotationMatrix()
	back := QuatFromRotationMatrix(m)
	// q and -q represent the same rotation.
	if !back.Equal(q, 1e-9) && !back.Equal(q.Neg(), 1e-9) {
		t.Errorf("matrix round trip = %+v want %+v (or negation)", back, q)
	}
	// Matrix must be orthonormal: applying to basis then comparing lengths.
	col0 := q.RotateVector(V3(1, 0, 0))
	if math.Abs(col0.Norm()-1) > 1e-9 {
		t.Errorf("rotation not length preserving: %v", col0.Norm())
	}
}

func TestQuaternionExpLog(t *testing.T) {
	// exp(0 + (pi/2)k) = k.
	q := Quat(0, 0, 0, math.Pi/2)
	if got := q.Exp(); !got.Equal(Quat(0, 0, 0, 1), 1e-12) {
		t.Errorf("exp((pi/2)k) = %+v want k", got)
	}
	// log(exp(q)) round trip for a unit rotation quaternion.
	u := QuatFromAxisAngle(V3(0, 1, 0), 0.7)
	if got := u.Log().Exp(); !got.Equal(u, 1e-9) {
		t.Errorf("exp(log(u)) = %+v want %+v", got, u)
	}
}

func TestQuaternionPow(t *testing.T) {
	// (45 deg about Z)^2 = 90 deg about Z.
	q := QuatFromAxisAngle(V3(0, 0, 1), math.Pi/4)
	want := QuatFromAxisAngle(V3(0, 0, 1), math.Pi/2)
	if got := q.Pow(2); !got.Equal(want, 1e-9) {
		t.Errorf("q^2 = %+v want %+v", got, want)
	}
	// q^0 = identity.
	if got := q.Pow(0); !got.Equal(IdentityQuat(), 1e-9) {
		t.Errorf("q^0 = %+v want identity", got)
	}
}

func TestSlerp(t *testing.T) {
	q0 := IdentityQuat()
	q1 := QuatFromAxisAngle(V3(0, 0, 1), math.Pi/2)
	// Endpoints.
	if got := Slerp(q0, q1, 0); !got.Equal(q0, 1e-9) {
		t.Errorf("slerp t=0 = %+v want %+v", got, q0)
	}
	if got := Slerp(q0, q1, 1); !got.Equal(q1, 1e-9) {
		t.Errorf("slerp t=1 = %+v want %+v", got, q1)
	}
	// Midpoint is 45 deg about Z.
	mid := QuatFromAxisAngle(V3(0, 0, 1), math.Pi/4)
	if got := Slerp(q0, q1, 0.5); !got.Equal(mid, 1e-9) {
		t.Errorf("slerp t=0.5 = %+v want %+v", got, mid)
	}
}

func TestQuatIsUnit(t *testing.T) {
	if !QuatFromAxisAngle(V3(1, 0, 0), 0.5).IsUnit(1e-12) {
		t.Errorf("axis-angle quaternion should be unit")
	}
	if Quat(2, 0, 0, 0).IsUnit(1e-9) {
		t.Errorf("2 should not be unit")
	}
}

func BenchmarkSlerp(b *testing.B) {
	q0 := QuatFromEuler(0.1, 0.2, 0.3)
	q1 := QuatFromEuler(1.0, -0.5, 0.7)
	var sink Quaternion
	for i := 0; i < b.N; i++ {
		sink = Slerp(q0, q1, float64(i%100)/100)
	}
	_ = sink
}
