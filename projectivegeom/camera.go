package projectivegeom

import "math"

// Mat3x4 is a 3x4 real matrix in row-major order, the shape of a pinhole camera
// projection matrix mapping homogeneous points of RP^3 to homogeneous points of
// RP^2.
type Mat3x4 [3][4]float64

// MulVec4 returns the 3-vector product m*v.
func (m Mat3x4) MulVec4(v Vec4) Vec3 {
	c := [4]float64{v.X, v.Y, v.Z, v.W}
	var o [3]float64
	for i := 0; i < 3; i++ {
		var s float64
		for k := 0; k < 4; k++ {
			s += m[i][k] * c[k]
		}
		o[i] = s
	}
	return Vec3{o[0], o[1], o[2]}
}

// Camera is a finite pinhole camera given by a 3x4 projection matrix P = K[R|t],
// mapping scene points of RP^3 to image points of RP^2.
type Camera struct {
	P Mat3x4
}

// NewCamera wraps a 3x4 matrix as a camera.
func NewCamera(p Mat3x4) Camera { return Camera{p} }

// CanonicalCamera returns the camera [I | 0] whose center is the origin and
// whose axes coincide with the world axes.
func CanonicalCamera() Camera {
	return Camera{Mat3x4{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}}}
}

// Intrinsics returns the upper-triangular camera calibration matrix K with focal
// lengths fx, fy, principal point (cx, cy) and axis skew s.
func Intrinsics(fx, fy, s, cx, cy float64) Mat3 {
	return NewMat3(fx, s, cx, 0, fy, cy, 0, 0, 1)
}

// CameraFromKRt builds the camera P = K[R | t] from calibration matrix K,
// rotation R and translation t.
func CameraFromKRt(k, r Mat3, t Vec3) Camera {
	var rt Mat3x4
	for i := 0; i < 3; i++ {
		rt[i][0] = r[i][0]
		rt[i][1] = r[i][1]
		rt[i][2] = r[i][2]
	}
	rt[0][3] = t.X
	rt[1][3] = t.Y
	rt[2][3] = t.Z
	var p Mat3x4
	for i := 0; i < 3; i++ {
		for j := 0; j < 4; j++ {
			var sum float64
			for l := 0; l < 3; l++ {
				sum += k[i][l] * rt[l][j]
			}
			p[i][j] = sum
		}
	}
	return Camera{p}
}

// Project returns the image of a scene point under the camera, as a homogeneous
// point of RP^2.
func (c Camera) Project(x SPoint) Point { return Point{c.P.MulVec4(x.V)} }

// ProjectAffine returns the pixel coordinates of a scene point, or ErrAtInfinity
// when the projected point lies on the image line at infinity (the point is on
// the camera's principal plane).
func (c Camera) ProjectAffine(x SPoint) (u, v float64, err error) {
	return c.Project(x).Affine()
}

// Center returns the camera center, the scene point C with P*C = 0, computed as
// the null vector of P via its four 3x3 minors.
func (c Camera) Center() SPoint {
	col := func(j0, j1, j2 int) float64 {
		return det3(
			c.P[0][j0], c.P[0][j1], c.P[0][j2],
			c.P[1][j0], c.P[1][j1], c.P[1][j2],
			c.P[2][j0], c.P[2][j1], c.P[2][j2],
		)
	}
	return SPoint{Vec4{
		X: +col(1, 2, 3),
		Y: -col(0, 2, 3),
		Z: +col(0, 1, 3),
		W: -col(0, 1, 2),
	}}
}

// PrincipalPlane returns the plane of scene points that project onto the image
// line at infinity, namely the third row of P.
func (c Camera) PrincipalPlane() SPlane {
	return SPlane{Vec4{c.P[2][0], c.P[2][1], c.P[2][2], c.P[2][3]}}
}

// LookAt builds a camera positioned at eye, viewing target, with the given world
// up direction and calibration K. It constructs a right-handed camera frame
// (right, down, forward) and returns P = K[R | -R*eye]. It returns
// ErrConfiguration when eye and target coincide or up is parallel to the view
// direction.
func LookAt(k Mat3, eye, target, up Vec3) (Camera, error) {
	fwd, ok := target.Sub(eye).Normalized()
	if !ok {
		return Camera{}, ErrConfiguration
	}
	right, ok := fwd.Cross(up).Normalized()
	if !ok {
		return Camera{}, ErrConfiguration
	}
	down := fwd.Cross(right)
	// Rows of R are the camera axes expressed in world coordinates.
	r := NewMat3(
		right.X, right.Y, right.Z,
		down.X, down.Y, down.Z,
		fwd.X, fwd.Y, fwd.Z,
	)
	t := r.MulVec(eye).Neg()
	return CameraFromKRt(k, r, t), nil
}

// Depth returns the signed depth (cheirality value) of a scene point: it is
// positive when the finite point lies in front of the camera for a camera with
// positive-determinant rotation block. It returns 0 for points at infinity.
func (c Camera) Depth(x SPoint) float64 {
	img := c.P.MulVec4(x.V)
	m := NewMat3(
		c.P[0][0], c.P[0][1], c.P[0][2],
		c.P[1][0], c.P[1][1], c.P[1][2],
		c.P[2][0], c.P[2][1], c.P[2][2],
	)
	det := m.Det()
	w := x.V.W
	if math.Abs(w) < Eps {
		return 0
	}
	sign := 1.0
	if det < 0 {
		sign = -1
	}
	m3 := Vec3{m[2][0], m[2][1], m[2][2]}.Norm()
	if m3 < Eps {
		return 0
	}
	return sign * img.Z / (w * m3)
}
