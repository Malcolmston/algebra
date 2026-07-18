package diffgeo

import "math"

// Curve is a parametric space curve: a function mapping a scalar parameter t to
// a point r(t) in three-dimensional space. Derivatives are taken with respect
// to t by central finite differences, so any smooth parametrization works.
type Curve func(t float64) Vec3

// Finite-difference step sizes tuned to balance truncation error against
// floating-point cancellation for each derivative order.
const (
	diffgeoCurveH1 = 1e-5
	diffgeoCurveH2 = 1e-4
	diffgeoCurveH3 = 1e-3
)

// Velocity returns the first derivative r'(t), the velocity vector, computed by
// a central difference. Its direction is tangent to the curve and its magnitude
// is the parametric [Speed].
func Velocity(c Curve, t float64) Vec3 {
	const h = diffgeoCurveH1
	return c(t + h).Sub(c(t - h)).Scale(1 / (2 * h))
}

// Acceleration returns the second derivative r”(t) by a central difference.
func Acceleration(c Curve, t float64) Vec3 {
	const h = diffgeoCurveH2
	return c(t + h).Sub(c(t).Scale(2)).Add(c(t - h)).Scale(1 / (h * h))
}

// Jerk returns the third derivative r”'(t), used by [Torsion], via a five-point
// central difference.
func Jerk(c Curve, t float64) Vec3 {
	const h = diffgeoCurveH3
	// (r(t+2h) - 2r(t+h) + 2r(t-h) - r(t-2h)) / (2h³)
	return c(t + 2*h).
		Sub(c(t + h).Scale(2)).
		Add(c(t - h).Scale(2)).
		Sub(c(t - 2*h)).
		Scale(1 / (2 * h * h * h))
}

// Speed returns the parametric speed |r'(t)|, the rate at which arc length
// accumulates per unit parameter.
func Speed(c Curve, t float64) float64 {
	return Velocity(c, t).Norm()
}

// UnitTangent returns the unit tangent vector T(t) = r'(t)/|r'(t)|. At a point
// where the velocity vanishes it returns the zero vector.
func UnitTangent(c Curve, t float64) Vec3 {
	return Velocity(c, t).Normalize()
}

// CurvatureVector returns dT/ds, the derivative of the unit tangent with
// respect to arc length. Its magnitude is the [Curvature] and its direction is
// the principal normal. It is computed as (a − (a·T)T)/|v|² with v = r'(t),
// a = r”(t) and T = v/|v|.
func CurvatureVector(c Curve, t float64) Vec3 {
	v := Velocity(c, t)
	a := Acceleration(c, t)
	s2 := v.Norm2()
	if s2 < Eps {
		return Vec3{}
	}
	tHat := v.Scale(1 / math.Sqrt(s2))
	// component of acceleration orthogonal to the tangent, divided by speed²
	ortho := a.Sub(tHat.Scale(a.Dot(tHat)))
	return ortho.Scale(1 / s2)
}

// Curvature returns the (unsigned) curvature κ(t) = |r'×r”| / |r'|³, the
// reciprocal of the radius of the osculating circle. For a straight line it is
// zero; for a circle of radius R it is 1/R.
func Curvature(c Curve, t float64) float64 {
	v := Velocity(c, t)
	a := Acceleration(c, t)
	sp := v.Norm()
	if sp < Eps {
		return 0
	}
	return v.Cross(a).Norm() / (sp * sp * sp)
}

// RadiusOfCurvature returns 1/κ(t), the radius of the osculating circle. It
// returns +Inf where the curvature vanishes.
func RadiusOfCurvature(c Curve, t float64) float64 {
	k := Curvature(c, t)
	if k < Eps {
		return math.Inf(1)
	}
	return 1 / k
}

// Torsion returns the torsion τ(t) = ((r'×r”)·r”') / |r'×r”|², which measures
// how sharply the curve twists out of its osculating plane. It is zero for any
// planar curve. Where the curvature vanishes the osculating plane is undefined
// and Torsion returns 0.
func Torsion(c Curve, t float64) float64 {
	v := Velocity(c, t)
	a := Acceleration(c, t)
	j := Jerk(c, t)
	cross := v.Cross(a)
	den := cross.Norm2()
	if den < Eps {
		return 0
	}
	return cross.Dot(j) / den
}

// PrincipalNormal returns the unit principal normal N(t), the direction in
// which the curve is turning. It is the normalized [CurvatureVector].
func PrincipalNormal(c Curve, t float64) Vec3 {
	return CurvatureVector(c, t).Normalize()
}

// Binormal returns the unit binormal B(t) = T×N = (r'×r”)/|r'×r”|, the normal
// to the osculating plane.
func Binormal(c Curve, t float64) Vec3 {
	return Velocity(c, t).Cross(Acceleration(c, t)).Normalize()
}

// Frame is an orthonormal Frenet frame at a point of a curve: the unit tangent
// T, principal normal N and binormal B, forming a right-handed basis with
// T×N = B.
type Frame struct {
	T, N, B Vec3
}

// FrenetFrame returns the [Frame] {T, N, B} at parameter t. The tangent is
// r'/|r'|, the binormal is (r'×r”)/|r'×r”|, and the normal is B×T. Where the
// curvature vanishes the normal and binormal are not well defined and are
// returned as zero vectors.
func FrenetFrame(c Curve, t float64) Frame {
	tHat := UnitTangent(c, t)
	b := Binormal(c, t)
	if b.IsZero(Eps) {
		return Frame{T: tHat}
	}
	return Frame{T: tHat, N: b.Cross(tHat), B: b}
}

// ArcLength returns the arc length of the curve between parameters a and b,
// computed as the integral of [Speed] using composite Simpson's rule with n
// subintervals. n is rounded up to the next even number and forced to at least
// 2. The result is the exact geometric length in the limit of large n.
func ArcLength(c Curve, a, b float64, n int) float64 {
	if n < 2 {
		n = 2
	}
	if n%2 == 1 {
		n++
	}
	h := (b - a) / float64(n)
	sum := Speed(c, a) + Speed(c, b)
	for i := 1; i < n; i++ {
		x := a + float64(i)*h
		if i%2 == 1 {
			sum += 4 * Speed(c, x)
		} else {
			sum += 2 * Speed(c, x)
		}
	}
	return math.Abs(h) / 3 * sum
}

// Line returns the straight-line curve r(t) = p + t·d through point p with
// direction d. Its curvature and torsion are identically zero.
func Line(p, d Vec3) Curve {
	return func(t float64) Vec3 {
		return p.Add(d.Scale(t))
	}
}

// Circle returns the planar circle of the given radius in the xy-plane,
// r(t) = (radius·cos t, radius·sin t, 0). Its curvature is 1/radius and its
// torsion is zero.
func Circle(radius float64) Curve {
	return func(t float64) Vec3 {
		return Vec3{radius * math.Cos(t), radius * math.Sin(t), 0}
	}
}

// Ellipse returns the planar ellipse r(t) = (a·cos t, b·sin t, 0) with the
// given semi-axes a and b.
func Ellipse(a, b float64) Curve {
	return func(t float64) Vec3 {
		return Vec3{a * math.Cos(t), b * math.Sin(t), 0}
	}
}

// Helix returns the circular helix r(t) = (radius·cos t, radius·sin t, rise·t).
// It has constant curvature radius/(radius²+rise²) and constant torsion
// rise/(radius²+rise²).
func Helix(radius, rise float64) Curve {
	return func(t float64) Vec3 {
		return Vec3{radius * math.Cos(t), radius * math.Sin(t), rise * t}
	}
}
