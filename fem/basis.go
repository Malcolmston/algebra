package fem

// ShapeP1Interval returns the two linear Lagrange basis functions evaluated at
// the reference coordinate xi in [0,1]. Index 0 corresponds to xi=0 and index 1
// to xi=1.
func ShapeP1Interval(xi float64) [2]float64 {
	return [2]float64{1 - xi, xi}
}

// GradShapeP1Interval returns the derivatives d/dxi of the linear basis
// functions on the reference interval. They are constant.
func GradShapeP1Interval() [2]float64 {
	return [2]float64{-1, 1}
}

// ShapeP2Interval returns the three quadratic Lagrange basis functions on the
// reference interval [0,1]. Index 0 is the node at xi=0, index 1 the node at
// xi=1 and index 2 the midpoint node at xi=0.5.
func ShapeP2Interval(xi float64) [3]float64 {
	return [3]float64{
		1 - 3*xi + 2*xi*xi,
		2*xi*xi - xi,
		4 * xi * (1 - xi),
	}
}

// GradShapeP2Interval returns the derivatives d/dxi of the quadratic basis
// functions on the reference interval at xi.
func GradShapeP2Interval(xi float64) [3]float64 {
	return [3]float64{
		-3 + 4*xi,
		4*xi - 1,
		4 - 8*xi,
	}
}

// Barycentric returns the barycentric coordinates (L1,L2,L3) of the point p
// with respect to the triangle with vertices v1,v2,v3.
func Barycentric(v1, v2, v3, p [2]float64) [3]float64 {
	det := (v2[1]-v3[1])*(v1[0]-v3[0]) + (v3[0]-v2[0])*(v1[1]-v3[1])
	l1 := ((v2[1]-v3[1])*(p[0]-v3[0]) + (v3[0]-v2[0])*(p[1]-v3[1])) / det
	l2 := ((v3[1]-v1[1])*(p[0]-v3[0]) + (v1[0]-v3[0])*(p[1]-v3[1])) / det
	l3 := 1 - l1 - l2
	return [3]float64{l1, l2, l3}
}

// ShapeP1Triangle returns the three linear basis functions given the
// barycentric coordinates l = (L1,L2,L3). For P1 the basis functions equal the
// barycentric coordinates.
func ShapeP1Triangle(l [3]float64) [3]float64 {
	return l
}

// GradShapeP1TriangleRef returns the reference gradients (d/dxi, d/deta) of the
// linear triangle basis functions on the reference triangle. They are constant.
func GradShapeP1TriangleRef() [3][2]float64 {
	return [3][2]float64{
		{-1, -1},
		{1, 0},
		{0, 1},
	}
}

// ShapeP2Triangle returns the six quadratic basis functions for barycentric
// coordinates l = (L1,L2,L3). Indices 0..2 are the vertex nodes and 3..5 the
// edge-midpoint nodes opposite v1, v2 and v3 respectively (edges (v2,v3),
// (v1,v3), (v1,v2)).
func ShapeP2Triangle(l [3]float64) [6]float64 {
	l1, l2, l3 := l[0], l[1], l[2]
	return [6]float64{
		l1 * (2*l1 - 1),
		l2 * (2*l2 - 1),
		l3 * (2*l3 - 1),
		4 * l2 * l3,
		4 * l1 * l3,
		4 * l1 * l2,
	}
}

// GradShapeP2TriangleRef returns the reference gradients (d/dxi, d/deta) of the
// six quadratic basis functions, evaluated at reference coordinates (xi,eta).
func GradShapeP2TriangleRef(xi, eta float64) [6][2]float64 {
	l1 := 1 - xi - eta
	l2 := xi
	l3 := eta
	return [6][2]float64{
		{(4*l1 - 1) * -1, (4*l1 - 1) * -1},
		{(4*l2 - 1) * 1, 0},
		{0, (4*l3 - 1) * 1},
		{4 * l3, 4 * l2},
		{-4 * l3, 4 * (l1 - l3)},
		{4 * (l1 - l2), -4 * l2},
	}
}

// TriangleGradients returns the constant physical gradients of the three
// barycentric coordinates (equivalently the P1 basis functions) on the triangle
// with the given vertices, together with the signed area.
func TriangleGradients(v1, v2, v3 [2]float64) (grads [3][2]float64, area float64) {
	area = TriangleSignedArea(v1, v2, v3)
	inv := 1.0 / (2 * area)
	grads[0] = [2]float64{(v2[1] - v3[1]) * inv, (v3[0] - v2[0]) * inv}
	grads[1] = [2]float64{(v3[1] - v1[1]) * inv, (v1[0] - v3[0]) * inv}
	grads[2] = [2]float64{(v1[1] - v2[1]) * inv, (v2[0] - v1[0]) * inv}
	return grads, area
}

// TriangleJacobian returns the 2×2 Jacobian of the affine map from the
// reference triangle to the physical triangle, its determinant and its inverse
// transpose (used to map reference gradients to physical gradients).
func TriangleJacobian(v1, v2, v3 [2]float64) (jac [2][2]float64, det float64, invT [2][2]float64) {
	jac = [2][2]float64{
		{v2[0] - v1[0], v3[0] - v1[0]},
		{v2[1] - v1[1], v3[1] - v1[1]},
	}
	det = jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0]
	invT = [2][2]float64{
		{jac[1][1] / det, -jac[1][0] / det},
		{-jac[0][1] / det, jac[0][0] / det},
	}
	return jac, det, invT
}

// PhysicalGradP2Triangle maps the reference gradients of the P2 basis functions
// at reference coordinates (xi,eta) to physical gradients on the triangle with
// the given vertices.
func PhysicalGradP2Triangle(v1, v2, v3 [2]float64, xi, eta float64) [6][2]float64 {
	_, _, invT := TriangleJacobian(v1, v2, v3)
	ref := GradShapeP2TriangleRef(xi, eta)
	var out [6][2]float64
	for i := 0; i < 6; i++ {
		gx := invT[0][0]*ref[i][0] + invT[0][1]*ref[i][1]
		gy := invT[1][0]*ref[i][0] + invT[1][1]*ref[i][1]
		out[i] = [2]float64{gx, gy}
	}
	return out
}
