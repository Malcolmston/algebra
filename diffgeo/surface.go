package diffgeo

import "math"

// Surface is a parametric surface patch: a function mapping parameters (u, v)
// to a point r(u, v) in three-dimensional space. Partial derivatives are taken
// by central finite differences, so any smooth parametrization works.
type Surface func(u, v float64) Vec3

// diffgeoSurfH is the finite-difference step used for surface partials. A
// single moderate step keeps the second-order mixed and pure partials
// numerically stable while leaving the first partials accurate.
const diffgeoSurfH = 1e-4

// PartialU returns the first partial derivative rᵤ = ∂r/∂u at (u, v).
func PartialU(s Surface, u, v float64) Vec3 {
	const h = diffgeoSurfH
	return s(u+h, v).Sub(s(u-h, v)).Scale(1 / (2 * h))
}

// PartialV returns the first partial derivative rᵥ = ∂r/∂v at (u, v).
func PartialV(s Surface, u, v float64) Vec3 {
	const h = diffgeoSurfH
	return s(u, v+h).Sub(s(u, v-h)).Scale(1 / (2 * h))
}

// PartialUU returns the second partial derivative rᵤᵤ = ∂²r/∂u² at (u, v).
func PartialUU(s Surface, u, v float64) Vec3 {
	const h = diffgeoSurfH
	return s(u+h, v).Sub(s(u, v).Scale(2)).Add(s(u-h, v)).Scale(1 / (h * h))
}

// PartialVV returns the second partial derivative rᵥᵥ = ∂²r/∂v² at (u, v).
func PartialVV(s Surface, u, v float64) Vec3 {
	const h = diffgeoSurfH
	return s(u, v+h).Sub(s(u, v).Scale(2)).Add(s(u, v-h)).Scale(1 / (h * h))
}

// PartialUV returns the mixed second partial derivative rᵤᵥ = ∂²r/∂u∂v at
// (u, v), computed by the four-point central difference.
func PartialUV(s Surface, u, v float64) Vec3 {
	const h = diffgeoSurfH
	return s(u+h, v+h).
		Sub(s(u+h, v-h)).
		Sub(s(u-h, v+h)).
		Add(s(u-h, v-h)).
		Scale(1 / (4 * h * h))
}

// SurfaceNormal returns the unit normal vector n = (rᵤ×rᵥ)/|rᵤ×rᵥ| at (u, v).
// The orientation follows the right-hand rule from the parameter ordering
// (u, v). At a degenerate point where rᵤ and rᵥ are parallel it returns the
// zero vector.
func SurfaceNormal(s Surface, u, v float64) Vec3 {
	return PartialU(s, u, v).Cross(PartialV(s, u, v)).Normalize()
}

// FirstForm is the first fundamental form (metric) of a surface at a point,
// I = E du² + 2F du dv + G dv², with coefficients E = rᵤ·rᵤ, F = rᵤ·rᵥ and
// G = rᵥ·rᵥ. It encodes lengths, angles and areas intrinsic to the surface.
type FirstForm struct {
	E, F, G float64
}

// Determinant returns EG − F², the squared area of the parallelogram spanned by
// rᵤ and rᵥ. It is strictly positive at a regular point.
func (f FirstForm) Determinant() float64 {
	return f.E*f.G - f.F*f.F
}

// FirstFundamental returns the [FirstForm] coefficients at (u, v).
func FirstFundamental(s Surface, u, v float64) FirstForm {
	ru := PartialU(s, u, v)
	rv := PartialV(s, u, v)
	return FirstForm{
		E: ru.Dot(ru),
		F: ru.Dot(rv),
		G: rv.Dot(rv),
	}
}

// SecondForm is the second fundamental form of a surface at a point,
// II = L du² + 2M du dv + N dv², with coefficients L = rᵤᵤ·n, M = rᵤᵥ·n and
// N = rᵥᵥ·n, where n is the unit [SurfaceNormal]. It measures how the surface
// curves away from its tangent plane.
type SecondForm struct {
	L, M, N float64
}

// SecondFundamental returns the [SecondForm] coefficients at (u, v). Its sign
// convention matches the orientation of [SurfaceNormal].
func SecondFundamental(s Surface, u, v float64) SecondForm {
	n := SurfaceNormal(s, u, v)
	return SecondForm{
		L: PartialUU(s, u, v).Dot(n),
		M: PartialUV(s, u, v).Dot(n),
		N: PartialVV(s, u, v).Dot(n),
	}
}

// GaussianCurvature returns the Gaussian curvature K = (LN − M²)/(EG − F²) at
// (u, v). It is an intrinsic invariant (Gauss's Theorema Egregium): positive on
// dome-like regions, negative on saddles, zero on developable surfaces such as
// planes and cylinders.
func GaussianCurvature(s Surface, u, v float64) float64 {
	I := FirstFundamental(s, u, v)
	II := SecondFundamental(s, u, v)
	den := I.Determinant()
	if math.Abs(den) < Eps {
		return 0
	}
	return (II.L*II.N - II.M*II.M) / den
}

// MeanCurvature returns the mean curvature H = (EN + GL − 2FM)/(2(EG − F²)) at
// (u, v), the average of the two [PrincipalCurvatures]. Its sign depends on the
// orientation of [SurfaceNormal]; a minimal surface has H = 0 everywhere.
func MeanCurvature(s Surface, u, v float64) float64 {
	I := FirstFundamental(s, u, v)
	II := SecondFundamental(s, u, v)
	den := I.Determinant()
	if math.Abs(den) < Eps {
		return 0
	}
	return (I.E*II.N + I.G*II.L - 2*I.F*II.M) / (2 * den)
}

// PrincipalCurvatures returns the principal curvatures k1 ≥ k2 at (u, v), the
// extreme values of the [NormalCurvature] over all tangent directions. They are
// the eigenvalues of the shape operator, k = H ± √(H² − K), where H is the
// [MeanCurvature] and K the [GaussianCurvature].
func PrincipalCurvatures(s Surface, u, v float64) (k1, k2 float64) {
	h := MeanCurvature(s, u, v)
	k := GaussianCurvature(s, u, v)
	disc := h*h - k
	if disc < 0 {
		disc = 0 // clamp tiny negative values from rounding
	}
	root := math.Sqrt(disc)
	return h + root, h - root
}

// NormalCurvature returns the normal curvature at (u, v) in the tangent
// direction (du, dv), that is II(du,dv)/I(du,dv) =
// (L du² + 2M du dv + N dv²)/(E du² + 2F du dv + G dv²). The direction need not
// be a unit vector; only its ratio matters.
func NormalCurvature(s Surface, u, v, du, dv float64) float64 {
	I := FirstFundamental(s, u, v)
	II := SecondFundamental(s, u, v)
	num := II.L*du*du + 2*II.M*du*dv + II.N*dv*dv
	den := I.E*du*du + 2*I.F*du*dv + I.G*dv*dv
	if math.Abs(den) < Eps {
		return 0
	}
	return num / den
}

// AreaElement returns the scalar surface area element √(EG − F²) = |rᵤ×rᵥ| at
// (u, v), the local factor relating parameter area to surface area.
func AreaElement(s Surface, u, v float64) float64 {
	d := FirstFundamental(s, u, v).Determinant()
	if d < 0 {
		d = 0
	}
	return math.Sqrt(d)
}

// SurfaceArea returns the area of the patch over the parameter rectangle
// [u0,u1]×[v0,v1], the double integral of the [AreaElement] evaluated by a
// tensor-product composite Simpson's rule with n subintervals in each
// direction. n is rounded up to the next even number and forced to at least 2.
func SurfaceArea(s Surface, u0, u1, v0, v1 float64, n int) float64 {
	if n < 2 {
		n = 2
	}
	if n%2 == 1 {
		n++
	}
	hu := (u1 - u0) / float64(n)
	hv := (v1 - v0) / float64(n)
	simpsonWeight := func(i int) float64 {
		if i == 0 || i == n {
			return 1
		}
		if i%2 == 1 {
			return 4
		}
		return 2
	}
	var sum float64
	for i := 0; i <= n; i++ {
		u := u0 + float64(i)*hu
		wu := simpsonWeight(i)
		for j := 0; j <= n; j++ {
			v := v0 + float64(j)*hv
			wv := simpsonWeight(j)
			sum += wu * wv * AreaElement(s, u, v)
		}
	}
	return math.Abs(hu*hv) / 9 * sum
}

// Christoffel holds the Christoffel symbols of the second kind Γᵏᵢⱼ at a point,
// indexed as Gamma[k][i][j] with parameter indices 0 = u and 1 = v. They are
// symmetric in the lower indices (Γᵏᵢⱼ = Γᵏⱼᵢ) and are the intrinsic connection
// coefficients governing geodesics and parallel transport.
type Christoffel struct {
	Gamma [2][2][2]float64
}

// At returns the symbol Γᵏᵢⱼ = Gamma[k][i][j].
func (c Christoffel) At(k, i, j int) float64 {
	return c.Gamma[k][i][j]
}

// ChristoffelSymbols returns the [Christoffel] symbols of the second kind at
// (u, v). They are derived purely from the first fundamental form and its first
// partial derivatives via Γᵏᵢⱼ = ½ gᵏˡ(∂ᵢgⱼˡ + ∂ⱼgᵢˡ − ∂ˡgᵢⱼ), where g is the
// metric [FirstForm] written as [[E,F],[F,G]] and gᵏˡ its inverse. The metric
// derivatives are taken by central differences of [FirstFundamental].
func ChristoffelSymbols(s Surface, u, v float64) Christoffel {
	const h = diffgeoSurfH

	// metric components g[i][j] as a function of parameters
	metric := func(uu, vv float64) [2][2]float64 {
		I := FirstFundamental(s, uu, vv)
		return [2][2]float64{{I.E, I.F}, {I.F, I.G}}
	}

	g := metric(u, v)
	det := g[0][0]*g[1][1] - g[0][1]*g[1][0]
	var c Christoffel
	if math.Abs(det) < Eps {
		return c
	}
	// inverse metric gInv[k][l]
	gInv := [2][2]float64{
		{g[1][1] / det, -g[0][1] / det},
		{-g[1][0] / det, g[0][0] / det},
	}

	// dg[m] = ∂g/∂x^m via central differences, m = 0 (u), 1 (v)
	gU := metric(u+h, v)
	gUm := metric(u-h, v)
	gV := metric(u, v+h)
	gVm := metric(u, v-h)
	var dg [2][2][2]float64 // dg[m][i][j]
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			dg[0][i][j] = (gU[i][j] - gUm[i][j]) / (2 * h)
			dg[1][i][j] = (gV[i][j] - gVm[i][j]) / (2 * h)
		}
	}

	for k := 0; k < 2; k++ {
		for i := 0; i < 2; i++ {
			for j := 0; j < 2; j++ {
				var sum float64
				for l := 0; l < 2; l++ {
					sum += gInv[k][l] * (dg[i][j][l] + dg[j][i][l] - dg[l][i][j])
				}
				c.Gamma[k][i][j] = 0.5 * sum
			}
		}
	}
	return c
}

// Sphere returns the parametrization of a sphere of the given radius by
// longitude u and latitude v: r(u,v) = radius·(cos v cos u, cos v sin u, sin v),
// with u in [0, 2π) and v in (−π/2, π/2). Its Gaussian curvature is
// 1/radius² everywhere.
func Sphere(radius float64) Surface {
	return func(u, v float64) Vec3 {
		cu, su := math.Cos(u), math.Sin(u)
		cv, sv := math.Cos(v), math.Sin(v)
		return Vec3{radius * cv * cu, radius * cv * su, radius * sv}
	}
}

// Cylinder returns the parametrization of a circular cylinder of the given
// radius aligned with the z-axis: r(u,v) = (radius·cos u, radius·sin u, v). It
// is a developable surface with zero Gaussian curvature.
func Cylinder(radius float64) Surface {
	return func(u, v float64) Vec3 {
		return Vec3{radius * math.Cos(u), radius * math.Sin(u), v}
	}
}

// Torus returns the parametrization of a torus with center-circle radius major
// and tube radius minor: r(u,v) = ((major+minor·cos v)cos u,
// (major+minor·cos v)sin u, minor·sin v). Its Gaussian curvature is
// cos v /(minor·(major+minor·cos v)).
func Torus(major, minor float64) Surface {
	return func(u, v float64) Vec3 {
		w := major + minor*math.Cos(v)
		return Vec3{w * math.Cos(u), w * math.Sin(u), minor * math.Sin(v)}
	}
}

// Graph returns the Monge patch r(u,v) = (u, v, f(u,v)) of the height function
// f, representing the surface z = f(x, y).
func Graph(f func(u, v float64) float64) Surface {
	return func(u, v float64) Vec3 {
		return Vec3{u, v, f(u, v)}
	}
}

// SurfacePlane returns the affine plane r(u,v) = origin + u·du + v·dv spanned by
// the direction vectors du and dv through origin. Provided the spanning vectors
// are independent, its Gaussian and mean curvatures are zero everywhere.
func SurfacePlane(origin, du, dv Vec3) Surface {
	return func(u, v float64) Vec3 {
		return origin.Add(du.Scale(u)).Add(dv.Scale(v))
	}
}
