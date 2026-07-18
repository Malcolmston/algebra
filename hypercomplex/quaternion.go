package hypercomplex

import "math"

// Quaternion is a Hamilton quaternion q = W + X*i + Y*j + Z*k, where i, j and k
// satisfy i*i = j*j = k*k = i*j*k = -1. W is the scalar (real) part and
// (X, Y, Z) is the vector (imaginary) part.
type Quaternion struct {
	W, X, Y, Z float64
}

// Quat constructs a quaternion from its scalar part w and vector part
// (x, y, z).
func Quat(w, x, y, z float64) Quaternion {
	return Quaternion{W: w, X: x, Y: y, Z: z}
}

// IdentityQuat returns the multiplicative identity quaternion 1 + 0i + 0j + 0k.
func IdentityQuat() Quaternion {
	return Quaternion{W: 1}
}

// QuatFromScalar returns the quaternion whose scalar part is s and whose vector
// part is zero.
func QuatFromScalar(s float64) Quaternion {
	return Quaternion{W: s}
}

// QuatFromVector returns the pure quaternion 0 + v (a quaternion with zero
// scalar part), used to embed a 3-vector into the quaternions.
func QuatFromVector(v Vec3) Quaternion {
	return Quaternion{X: v.X, Y: v.Y, Z: v.Z}
}

// QuatFromAxisAngle returns the unit quaternion representing a rotation of angle
// radians about axis. The axis is normalized internally; a zero axis yields the
// identity quaternion.
func QuatFromAxisAngle(axis Vec3, angle float64) Quaternion {
	u := axis.Normalize()
	if u == (Vec3{}) {
		return IdentityQuat()
	}
	half := angle / 2
	s := math.Sin(half)
	return Quaternion{
		W: math.Cos(half),
		X: u.X * s,
		Y: u.Y * s,
		Z: u.Z * s,
	}
}

// QuatFromEuler returns the unit quaternion for the intrinsic Tait-Bryan
// rotation with the given roll (about X), pitch (about Y) and yaw (about Z),
// applied in the aerospace Z-Y-X order (yaw, then pitch, then roll). Angles are
// in radians.
func QuatFromEuler(roll, pitch, yaw float64) Quaternion {
	cr, sr := math.Cos(roll/2), math.Sin(roll/2)
	cp, sp := math.Cos(pitch/2), math.Sin(pitch/2)
	cy, sy := math.Cos(yaw/2), math.Sin(yaw/2)
	return Quaternion{
		W: cr*cp*cy + sr*sp*sy,
		X: sr*cp*cy - cr*sp*sy,
		Y: cr*sp*cy + sr*cp*sy,
		Z: cr*cp*sy - sr*sp*cy,
	}
}

// QuatFromRotationMatrix returns the unit quaternion corresponding to the
// row-major 3x3 rotation matrix m (nine elements m00, m01, m02, m10, ...). The
// implementation uses the numerically stable trace method and normalizes the
// result.
func QuatFromRotationMatrix(m [9]float64) Quaternion {
	m00, m01, m02 := m[0], m[1], m[2]
	m10, m11, m12 := m[3], m[4], m[5]
	m20, m21, m22 := m[6], m[7], m[8]
	trace := m00 + m11 + m22
	var q Quaternion
	switch {
	case trace > 0:
		s := math.Sqrt(trace+1) * 2 // s = 4*W
		q.W = 0.25 * s
		q.X = (m21 - m12) / s
		q.Y = (m02 - m20) / s
		q.Z = (m10 - m01) / s
	case m00 > m11 && m00 > m22:
		s := math.Sqrt(1+m00-m11-m22) * 2 // s = 4*X
		q.W = (m21 - m12) / s
		q.X = 0.25 * s
		q.Y = (m01 + m10) / s
		q.Z = (m02 + m20) / s
	case m11 > m22:
		s := math.Sqrt(1+m11-m00-m22) * 2 // s = 4*Y
		q.W = (m02 - m20) / s
		q.X = (m01 + m10) / s
		q.Y = 0.25 * s
		q.Z = (m12 + m21) / s
	default:
		s := math.Sqrt(1+m22-m00-m11) * 2 // s = 4*Z
		q.W = (m10 - m01) / s
		q.X = (m02 + m20) / s
		q.Y = (m12 + m21) / s
		q.Z = 0.25 * s
	}
	return q.Normalize()
}

// Scalar returns the scalar (real) part W of q.
func (q Quaternion) Scalar() float64 {
	return q.W
}

// Vector returns the vector (imaginary) part (X, Y, Z) of q as a [Vec3].
func (q Quaternion) Vector() Vec3 {
	return Vec3{q.X, q.Y, q.Z}
}

// Add returns the component-wise sum q + r.
func (q Quaternion) Add(r Quaternion) Quaternion {
	return Quaternion{q.W + r.W, q.X + r.X, q.Y + r.Y, q.Z + r.Z}
}

// Sub returns the component-wise difference q - r.
func (q Quaternion) Sub(r Quaternion) Quaternion {
	return Quaternion{q.W - r.W, q.X - r.X, q.Y - r.Y, q.Z - r.Z}
}

// Scale returns q with every component multiplied by the real factor s.
func (q Quaternion) Scale(s float64) Quaternion {
	return Quaternion{q.W * s, q.X * s, q.Y * s, q.Z * s}
}

// Neg returns the additive inverse -q.
func (q Quaternion) Neg() Quaternion {
	return Quaternion{-q.W, -q.X, -q.Y, -q.Z}
}

// Mul returns the Hamilton product q*r. Quaternion multiplication is
// associative but not commutative.
func (q Quaternion) Mul(r Quaternion) Quaternion {
	return Quaternion{
		W: q.W*r.W - q.X*r.X - q.Y*r.Y - q.Z*r.Z,
		X: q.W*r.X + q.X*r.W + q.Y*r.Z - q.Z*r.Y,
		Y: q.W*r.Y - q.X*r.Z + q.Y*r.W + q.Z*r.X,
		Z: q.W*r.Z + q.X*r.Y - q.Y*r.X + q.Z*r.W,
	}
}

// Conj returns the conjugate q* = W - X*i - Y*j - Z*k.
func (q Quaternion) Conj() Quaternion {
	return Quaternion{q.W, -q.X, -q.Y, -q.Z}
}

// Dot returns the Euclidean dot product of q and r treated as 4-vectors.
func (q Quaternion) Dot(r Quaternion) float64 {
	return q.W*r.W + q.X*r.X + q.Y*r.Y + q.Z*r.Z
}

// NormSq returns the squared norm W^2 + X^2 + Y^2 + Z^2, equal to the scalar
// part of q * q*.
func (q Quaternion) NormSq() float64 {
	return q.W*q.W + q.X*q.X + q.Y*q.Y + q.Z*q.Z
}

// Norm returns the Euclidean norm |q| = sqrt(W^2 + X^2 + Y^2 + Z^2).
func (q Quaternion) Norm() float64 {
	return math.Sqrt(q.NormSq())
}

// Normalize returns q scaled to unit norm. A zero quaternion is returned
// unchanged.
func (q Quaternion) Normalize() Quaternion {
	n := q.Norm()
	if n == 0 {
		return q
	}
	return q.Scale(1 / n)
}

// Inverse returns the multiplicative inverse q^-1 = q* / |q|^2. A zero
// quaternion is returned unchanged.
func (q Quaternion) Inverse() Quaternion {
	n2 := q.NormSq()
	if n2 == 0 {
		return q
	}
	return q.Conj().Scale(1 / n2)
}

// IsUnit reports whether q has unit norm to within the absolute tolerance tol.
func (q Quaternion) IsUnit(tol float64) bool {
	return math.Abs(q.Norm()-1) <= tol
}

// Equal reports whether q and r agree component-wise to within the absolute
// tolerance tol.
func (q Quaternion) Equal(r Quaternion, tol float64) bool {
	return math.Abs(q.W-r.W) <= tol &&
		math.Abs(q.X-r.X) <= tol &&
		math.Abs(q.Y-r.Y) <= tol &&
		math.Abs(q.Z-r.Z) <= tol
}

// Exp returns the quaternion exponential exp(q). Writing q = w + v with v the
// vector part and theta = |v|, exp(q) = e^w (cos theta + (v/theta) sin theta).
func (q Quaternion) Exp() Quaternion {
	ew := math.Exp(q.W)
	v := q.Vector()
	theta := v.Norm()
	if theta == 0 {
		return Quaternion{W: ew}
	}
	s := ew * math.Sin(theta) / theta
	return Quaternion{
		W: ew * math.Cos(theta),
		X: v.X * s,
		Y: v.Y * s,
		Z: v.Z * s,
	}
}

// Log returns the principal quaternion logarithm log(q). Writing q = |q|
// (cos phi + n sin phi) with n the unit vector part, log(q) = ln|q| + n*phi,
// where phi = atan2(|v|, w). A zero quaternion is returned unchanged.
func (q Quaternion) Log() Quaternion {
	n := q.Norm()
	if n == 0 {
		return q
	}
	v := q.Vector()
	vn := v.Norm()
	if vn == 0 {
		// Purely real quaternion; principal log has zero vector part for w>0.
		return Quaternion{W: math.Log(n)}
	}
	phi := math.Atan2(vn, q.W)
	s := phi / vn
	return Quaternion{
		W: math.Log(n),
		X: v.X * s,
		Y: v.Y * s,
		Z: v.Z * s,
	}
}

// Pow returns q raised to the real power t, computed as exp(t*log(q)). For a
// unit quaternion this scales the rotation angle by t.
func (q Quaternion) Pow(t float64) Quaternion {
	return q.Log().Scale(t).Exp()
}

// RotateVector returns the image of v under the rotation represented by q,
// computed as the vector part of q * (0,v) * q^-1. The quaternion need not be a
// unit quaternion; any non-zero scaling cancels.
func (q Quaternion) RotateVector(v Vec3) Vec3 {
	p := QuatFromVector(v)
	r := q.Mul(p).Mul(q.Inverse())
	return r.Vector()
}

// ToAxisAngle returns the rotation axis (a unit vector) and angle in radians
// (in the range [0, pi]) represented by q. The quaternion is normalized
// internally. For the identity rotation the axis defaults to (1, 0, 0).
func (q Quaternion) ToAxisAngle() (axis Vec3, angle float64) {
	u := q.Normalize()
	// Fold to the hemisphere with non-negative scalar part for angle in [0,pi].
	if u.W < 0 {
		u = u.Neg()
	}
	w := math.Min(1, math.Max(-1, u.W))
	angle = 2 * math.Acos(w)
	v := u.Vector()
	vn := v.Norm()
	if vn < 1e-12 {
		return Vec3{1, 0, 0}, 0
	}
	return v.Scale(1 / vn), angle
}

// ToEuler returns the roll (about X), pitch (about Y) and yaw (about Z) angles
// in radians that reproduce the rotation of q under the aerospace Z-Y-X
// convention (the inverse of [QuatFromEuler]). The pitch is clamped to
// [-pi/2, pi/2].
func (q Quaternion) ToEuler() (roll, pitch, yaw float64) {
	u := q.Normalize()
	// Roll (x-axis rotation).
	sinrCosp := 2 * (u.W*u.X + u.Y*u.Z)
	cosrCosp := 1 - 2*(u.X*u.X+u.Y*u.Y)
	roll = math.Atan2(sinrCosp, cosrCosp)
	// Pitch (y-axis rotation), clamped for gimbal-lock safety.
	sinp := 2 * (u.W*u.Y - u.Z*u.X)
	if sinp >= 1 {
		pitch = math.Pi / 2
	} else if sinp <= -1 {
		pitch = -math.Pi / 2
	} else {
		pitch = math.Asin(sinp)
	}
	// Yaw (z-axis rotation).
	sinyCosp := 2 * (u.W*u.Z + u.X*u.Y)
	cosyCosp := 1 - 2*(u.Y*u.Y+u.Z*u.Z)
	yaw = math.Atan2(sinyCosp, cosyCosp)
	return roll, pitch, yaw
}

// ToRotationMatrix returns the row-major 3x3 rotation matrix (nine elements
// m00, m01, m02, m10, ...) equivalent to the rotation represented by q. The
// quaternion is normalized internally.
func (q Quaternion) ToRotationMatrix() [9]float64 {
	u := q.Normalize()
	w, x, y, z := u.W, u.X, u.Y, u.Z
	return [9]float64{
		1 - 2*(y*y+z*z), 2 * (x*y - w*z), 2 * (x*z + w*y),
		2 * (x*y + w*z), 1 - 2*(x*x+z*z), 2 * (y*z - w*x),
		2 * (x*z - w*y), 2 * (y*z + w*x), 1 - 2*(x*x+y*y),
	}
}

// Slerp returns the spherical linear interpolation between the unit quaternions
// q0 and q1 at parameter t in [0, 1]. The interpolation follows the shortest
// arc (choosing the sign of q1 that minimizes the angle) and falls back to
// normalized linear interpolation when the endpoints are nearly parallel.
func Slerp(q0, q1 Quaternion, t float64) Quaternion {
	a := q0.Normalize()
	b := q1.Normalize()
	dot := a.Dot(b)
	if dot < 0 {
		b = b.Neg()
		dot = -dot
	}
	const threshold = 0.9995
	if dot > threshold {
		// Nearly parallel: linear interpolate and renormalize.
		return a.Add(b.Sub(a).Scale(t)).Normalize()
	}
	dot = math.Min(1, math.Max(-1, dot))
	theta0 := math.Acos(dot)
	sinTheta0 := math.Sin(theta0)
	theta := theta0 * t
	s0 := math.Cos(theta) - dot*math.Sin(theta)/sinTheta0
	s1 := math.Sin(theta) / sinTheta0
	return a.Scale(s0).Add(b.Scale(s1))
}
