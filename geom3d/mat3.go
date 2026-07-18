package geom3d

import (
	"fmt"
	"math"
)

// Mat3 is a 3x3 matrix of float64 stored in row-major order: element m[r][c]
// lives in row r and column c. Vectors are treated as column vectors and
// transformed on the left by [Mat3.MulVec].
type Mat3 [3][3]float64

// Mat3Identity returns the 3x3 identity matrix.
func Mat3Identity() Mat3 {
	return Mat3{
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1},
	}
}

// Mat3FromRows builds a matrix whose rows are the given vectors r0, r1 and r2.
func Mat3FromRows(r0, r1, r2 Vec3) Mat3 {
	return Mat3{
		{r0.X, r0.Y, r0.Z},
		{r1.X, r1.Y, r1.Z},
		{r2.X, r2.Y, r2.Z},
	}
}

// Mat3FromColumns builds a matrix whose columns are the given vectors c0, c1
// and c2.
func Mat3FromColumns(c0, c1, c2 Vec3) Mat3 {
	return Mat3{
		{c0.X, c1.X, c2.X},
		{c0.Y, c1.Y, c2.Y},
		{c0.Z, c1.Z, c2.Z},
	}
}

// Row returns row i (0, 1 or 2) of m as a vector. It panics if i is out of
// range.
func (m Mat3) Row(i int) Vec3 {
	return Vec3{m[i][0], m[i][1], m[i][2]}
}

// Col returns column j (0, 1 or 2) of m as a vector. It panics if j is out of
// range.
func (m Mat3) Col(j int) Vec3 {
	return Vec3{m[0][j], m[1][j], m[2][j]}
}

// Add returns the component-wise sum m+n.
func (m Mat3) Add(n Mat3) Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r[i][j] = m[i][j] + n[i][j]
		}
	}
	return r
}

// Sub returns the component-wise difference m-n.
func (m Mat3) Sub(n Mat3) Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r[i][j] = m[i][j] - n[i][j]
		}
	}
	return r
}

// Scale returns m with every element multiplied by the scalar s.
func (m Mat3) Scale(s float64) Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r[i][j] = m[i][j] * s
		}
	}
	return r
}

// Mul returns the matrix product m*n.
func (m Mat3) Mul(n Mat3) Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r[i][j] = m[i][0]*n[0][j] + m[i][1]*n[1][j] + m[i][2]*n[2][j]
		}
	}
	return r
}

// MulVec returns the matrix-vector product m*v, treating v as a column vector.
func (m Mat3) MulVec(v Vec3) Vec3 {
	return Vec3{
		m[0][0]*v.X + m[0][1]*v.Y + m[0][2]*v.Z,
		m[1][0]*v.X + m[1][1]*v.Y + m[1][2]*v.Z,
		m[2][0]*v.X + m[2][1]*v.Y + m[2][2]*v.Z,
	}
}

// Transpose returns the transpose of m.
func (m Mat3) Transpose() Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r[i][j] = m[j][i]
		}
	}
	return r
}

// Trace returns the sum of the diagonal elements of m.
func (m Mat3) Trace() float64 {
	return m[0][0] + m[1][1] + m[2][2]
}

// Determinant returns the determinant of m.
func (m Mat3) Determinant() float64 {
	return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) -
		m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) +
		m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
}

// Inverse returns the matrix inverse of m and true, or the zero matrix and
// false if m is singular (its determinant is zero within the package
// tolerance).
func (m Mat3) Inverse() (Mat3, bool) {
	det := m.Determinant()
	if math.Abs(det) <= geom3dEps {
		return Mat3{}, false
	}
	inv := 1 / det
	var r Mat3
	r[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1]) * inv
	r[0][1] = (m[0][2]*m[2][1] - m[0][1]*m[2][2]) * inv
	r[0][2] = (m[0][1]*m[1][2] - m[0][2]*m[1][1]) * inv
	r[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2]) * inv
	r[1][1] = (m[0][0]*m[2][2] - m[0][2]*m[2][0]) * inv
	r[1][2] = (m[0][2]*m[1][0] - m[0][0]*m[1][2]) * inv
	r[2][0] = (m[1][0]*m[2][1] - m[1][1]*m[2][0]) * inv
	r[2][1] = (m[0][1]*m[2][0] - m[0][0]*m[2][1]) * inv
	r[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0]) * inv
	return r, true
}

// Equal reports whether m and n are equal within the absolute tolerance eps,
// compared element-wise.
func (m Mat3) Equal(n Mat3, eps float64) bool {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if math.Abs(m[i][j]-n[i][j]) > eps {
				return false
			}
		}
	}
	return true
}

// String returns a human-readable, row-by-row representation of m.
func (m Mat3) String() string {
	return fmt.Sprintf("[[%g %g %g] [%g %g %g] [%g %g %g]]",
		m[0][0], m[0][1], m[0][2],
		m[1][0], m[1][1], m[1][2],
		m[2][0], m[2][1], m[2][2])
}

// RotationX returns the matrix that rotates a vector by angle radians about the
// x-axis, following the right-hand rule.
func RotationX(angle float64) Mat3 {
	c, s := math.Cos(angle), math.Sin(angle)
	return Mat3{
		{1, 0, 0},
		{0, c, -s},
		{0, s, c},
	}
}

// RotationY returns the matrix that rotates a vector by angle radians about the
// y-axis, following the right-hand rule.
func RotationY(angle float64) Mat3 {
	c, s := math.Cos(angle), math.Sin(angle)
	return Mat3{
		{c, 0, s},
		{0, 1, 0},
		{-s, 0, c},
	}
}

// RotationZ returns the matrix that rotates a vector by angle radians about the
// z-axis, following the right-hand rule.
func RotationZ(angle float64) Mat3 {
	c, s := math.Cos(angle), math.Sin(angle)
	return Mat3{
		{c, -s, 0},
		{s, c, 0},
		{0, 0, 1},
	}
}

// AxisAngleMatrix returns the rotation matrix for a rotation by angle radians
// about the given axis, using Rodrigues' rotation formula. The axis is
// normalized internally; if it has zero length the identity matrix is
// returned.
func AxisAngleMatrix(axis Vec3, angle float64) Mat3 {
	k, n := axis.Normalize()
	if n == 0 {
		return Mat3Identity()
	}
	c, s := math.Cos(angle), math.Sin(angle)
	t := 1 - c
	x, y, z := k.X, k.Y, k.Z
	return Mat3{
		{t*x*x + c, t*x*y - s*z, t*x*z + s*y},
		{t*x*y + s*z, t*y*y + c, t*y*z - s*x},
		{t*x*z - s*y, t*y*z + s*x, t*z*z + c},
	}
}

// Rodrigues rotates the vector v by angle radians about the given axis using
// Rodrigues' rotation formula, without forming a matrix. The axis is
// normalized internally; if it has zero length v is returned unchanged.
func Rodrigues(v, axis Vec3, angle float64) Vec3 {
	k, n := axis.Normalize()
	if n == 0 {
		return v
	}
	c, s := math.Cos(angle), math.Sin(angle)
	// v_rot = v*cos + (k x v)*sin + k*(k·v)*(1-cos)
	return v.Scale(c).
		Add(k.Cross(v).Scale(s)).
		Add(k.Scale(k.Dot(v) * (1 - c)))
}

// EulerAngles holds a set of Tait-Bryan angles in radians using the aerospace
// convention: Roll about the x-axis, Pitch about the y-axis and Yaw about the
// z-axis.
type EulerAngles struct {
	Roll, Pitch, Yaw float64
}

// EulerToMatrix returns the rotation matrix for the given Tait-Bryan angles
// using the intrinsic Z-Y-X convention, i.e. R = Rz(yaw)·Ry(pitch)·Rx(roll).
// Applying the result to a vector first rolls, then pitches, then yaws it.
func EulerToMatrix(e EulerAngles) Mat3 {
	return RotationZ(e.Yaw).Mul(RotationY(e.Pitch)).Mul(RotationX(e.Roll))
}

// MatrixToEuler extracts the Tait-Bryan angles (intrinsic Z-Y-X convention)
// that produce the rotation matrix m, inverting [EulerToMatrix]. Near the
// gimbal-lock poles (pitch = ±pi/2) roll is set to zero and yaw absorbs the
// combined rotation. The matrix is assumed to be a proper rotation.
func MatrixToEuler(m Mat3) EulerAngles {
	// m[2][0] = -sin(pitch)
	sp := -m[2][0]
	sp = geom3dclamp(sp, -1, 1)
	pitch := math.Asin(sp)
	var roll, yaw float64
	if math.Abs(sp) < 1-geom3dEps {
		roll = math.Atan2(m[2][1], m[2][2])
		yaw = math.Atan2(m[1][0], m[0][0])
	} else {
		// Gimbal lock: pitch is ±pi/2, only roll+yaw are determined.
		roll = 0
		yaw = math.Atan2(-m[0][1], m[1][1])
	}
	return EulerAngles{Roll: roll, Pitch: pitch, Yaw: yaw}
}
