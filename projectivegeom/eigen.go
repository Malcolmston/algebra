package projectivegeom

import "math"

// realCubicRoots returns the real roots of x^3 + a*x^2 + b*x + c = 0. It uses
// the depressed-cubic trigonometric and Cardano formulas and returns one or
// three roots.
func realCubicRoots(a, b, c float64) []float64 {
	// Depress: x = t - a/3, giving t^3 + p t + q = 0.
	p := b - a*a/3
	q := 2*a*a*a/27 - a*b/3 + c
	shift := -a / 3
	const eps = 1e-14
	if math.Abs(p) < eps && math.Abs(q) < eps {
		return []float64{shift}
	}
	disc := q*q/4 + p*p*p/27
	if disc > eps {
		// One real root.
		s := math.Sqrt(disc)
		u := math.Cbrt(-q/2 + s)
		v := math.Cbrt(-q/2 - s)
		return []float64{u + v + shift}
	}
	if disc < -eps {
		// Three distinct real roots (trigonometric form).
		m := 2 * math.Sqrt(-p/3)
		theta := math.Acos(3*q/(p*m)) / 3
		r0 := m*math.Cos(theta) + shift
		r1 := m*math.Cos(theta-2*math.Pi/3) + shift
		r2 := m*math.Cos(theta-4*math.Pi/3) + shift
		return []float64{r0, r1, r2}
	}
	// disc ~ 0: multiple real roots.
	if math.Abs(p) < eps {
		return []float64{shift}
	}
	r0 := 3 * q / p
	r1 := -3 * q / (2 * p)
	return []float64{r0 + shift, r1 + shift, r1 + shift}
}

// realEigenvalues3 returns the real eigenvalues of the 3x3 matrix m as roots of
// its characteristic polynomial.
func realEigenvalues3(m Mat3) []float64 {
	tr := m.Trace()
	// Sum of the three principal 2x2 minors.
	m00 := m[1][1]*m[2][2] - m[1][2]*m[2][1]
	m11 := m[0][0]*m[2][2] - m[0][2]*m[2][0]
	m22 := m[0][0]*m[1][1] - m[0][1]*m[1][0]
	minors := m00 + m11 + m22
	det := m.Det()
	// Characteristic poly: λ^3 - tr λ^2 + minors λ - det = 0.
	return realCubicRoots(-tr, minors, -det)
}

// eigenvector3 returns a non-zero eigenvector of m for the eigenvalue lambda by
// solving (m - lambda*I) v = 0. The second result is false when no numerically
// valid eigenvector is found.
func eigenvector3(m Mat3, lambda, tol float64) (Vec3, bool) {
	a := m
	a[0][0] -= lambda
	a[1][1] -= lambda
	a[2][2] -= lambda
	rows := [][]float64{
		{a[0][0], a[0][1], a[0][2]},
		{a[1][0], a[1][1], a[1][2]},
		{a[2][0], a[2][1], a[2][2]},
	}
	x, ok := nullVector(rows)
	if !ok {
		return Vec3{}, false
	}
	v := Vec3{x[0], x[1], x[2]}
	if v.IsZero(Eps) {
		return Vec3{}, false
	}
	// Verify residual.
	res := a.MulVec(v)
	if res.Norm() > tol*(1+v.Norm())*(1+m.MaxAbs()) {
		return Vec3{}, false
	}
	return v, true
}
