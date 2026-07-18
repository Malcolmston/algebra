package complexanalysis

import (
	"math"
	"math/cmplx"
)

// IntegrateCircle approximates the contour integral of f around the positively
// oriented circle of the given center and radius, using the composite
// trapezoidal rule with n sample points. Because the integrand is periodic,
// this rule converges geometrically when f is analytic in an annulus around the
// circle. It panics if n <= 0 or radius <= 0.
func IntegrateCircle(f Function, center complex128, radius float64, n int) complex128 {
	if n <= 0 {
		panic("complexanalysis: IntegrateCircle requires n > 0")
	}
	if radius <= 0 {
		panic("complexanalysis: IntegrateCircle requires radius > 0")
	}
	var sum complex128
	for j := 0; j < n; j++ {
		t := 2 * math.Pi * float64(j) / float64(n)
		e := cmplx.Rect(1, t)
		z := center + complex(radius, 0)*e
		// dz = radius * i * e^{i t} dt, dt = 2*pi/n.
		dz := complex(radius, 0) * complex(0, 1) * e * complex(2*math.Pi/float64(n), 0)
		sum += f(z) * dz
	}
	return sum
}

// IntegrateSegment approximates the integral of f along the straight line
// segment from a to b using the composite Simpson rule with n subintervals
// (n is rounded up to the next even number). It panics if n <= 0.
func IntegrateSegment(f Function, a, b complex128, n int) complex128 {
	if n <= 0 {
		panic("complexanalysis: IntegrateSegment requires n > 0")
	}
	if n%2 == 1 {
		n++
	}
	h := (b - a) / complex(float64(n), 0)
	sum := f(a) + f(b)
	for j := 1; j < n; j++ {
		z := a + complex(float64(j), 0)*h
		if j%2 == 1 {
			sum += 4 * f(z)
		} else {
			sum += 2 * f(z)
		}
	}
	return sum * h / 3
}

// IntegratePolygon approximates the integral of f along the closed polygonal
// contour through the given vertices (the last vertex is joined back to the
// first), applying IntegrateSegment with nPerEdge points on each edge. It
// returns 0 for fewer than two vertices.
func IntegratePolygon(f Function, vertices []complex128, nPerEdge int) complex128 {
	if len(vertices) < 2 {
		return 0
	}
	var sum complex128
	for i := 0; i < len(vertices); i++ {
		a := vertices[i]
		b := vertices[(i+1)%len(vertices)]
		sum += IntegrateSegment(f, a, b, nPerEdge)
	}
	return sum
}

// ContourIntegral approximates the integral of f along a general contour given
// by the parametrization gamma over t in [0, 1]. The derivative gamma'(t) is
// estimated by central finite differences, and the composite trapezoidal rule
// with n subintervals is used. If the contour is closed (gamma(0) == gamma(1)),
// the periodic trapezoidal rule is highly accurate. It panics if n <= 0.
func ContourIntegral(f Function, gamma func(float64) complex128, n int) complex128 {
	if n <= 0 {
		panic("complexanalysis: ContourIntegral requires n > 0")
	}
	const h = 1e-6
	dt := 1.0 / float64(n)
	integrand := func(t float64) complex128 {
		dgamma := (gamma(t+h) - gamma(t-h)) / complex(2*h, 0)
		return f(gamma(t)) * dgamma
	}
	// Trapezoidal rule on [0,1].
	sum := (integrand(0) + integrand(1)) / 2
	for j := 1; j < n; j++ {
		sum += integrand(float64(j) * dt)
	}
	return sum * complex(dt, 0)
}

// Residue approximates the residue of f at the isolated singularity z0 by
// integrating f around a small circle of the given radius and dividing by
// 2*pi*i. The radius must enclose no singularity other than z0.
func Residue(f Function, z0 complex128, radius float64, n int) complex128 {
	return IntegrateCircle(f, z0, radius, n) / complexanalysisTwoPiI
}

// ResidueSimplePole returns the residue of f at a simple pole z0, computed as
// the limit of (z-z0)f(z) as z approaches z0 via a small symmetric average of
// radius eps around z0.
func ResidueSimplePole(f Function, z0 complex128, eps float64) complex128 {
	g := func(z complex128) complex128 { return (z - z0) * f(z) }
	// Average four symmetric points to cancel the leading error term.
	pts := []complex128{
		z0 + complex(eps, 0),
		z0 - complex(eps, 0),
		z0 + complex(0, eps),
		z0 - complex(0, eps),
	}
	var sum complex128
	for _, p := range pts {
		sum += g(p)
	}
	return sum / 4
}

// ResidueOrderM returns the residue of f at a pole z0 of order m using the
// derivative formula Res = lim_{z->z0} 1/(m-1)! d^{m-1}/dz^{m-1} [(z-z0)^m f(z)].
// The (m-1)-th derivative is obtained from a Cauchy integral around a circle of
// the given radius, so f must be analytic on and inside that circle except at
// z0. It panics if m < 1.
func ResidueOrderM(f Function, z0 complex128, m int, radius float64, n int) complex128 {
	if m < 1 {
		panic("complexanalysis: ResidueOrderM requires m >= 1")
	}
	if m == 1 {
		return Residue(f, z0, radius, n)
	}
	g := func(z complex128) complex128 { return cmplx.Pow(z-z0, complex(float64(m), 0)) * f(z) }
	return CauchyDerivative(g, z0, z0, m-1, radius, n)
}

// CauchyIntegralValue evaluates f at the point z0 inside a circle of the given
// center and radius using the Cauchy integral formula,
// f(z0) = 1/(2*pi*i) * integral of f(z)/(z-z0) around the circle. The point z0
// must lie strictly inside the circle and f must be analytic there.
func CauchyIntegralValue(f Function, z0, center complex128, radius float64, n int) complex128 {
	g := func(z complex128) complex128 { return f(z) / (z - z0) }
	return IntegrateCircle(g, center, radius, n) / complexanalysisTwoPiI
}

// CauchyDerivative evaluates the k-th derivative of f at z0 using the Cauchy
// integral formula for derivatives,
// f^(k)(z0) = k!/(2*pi*i) * integral of f(z)/(z-z0)^(k+1) around the circle of
// the given center and radius. The point z0 must lie strictly inside the circle
// and f must be analytic there. k == 0 returns f(z0). It panics if k < 0.
func CauchyDerivative(f Function, z0, center complex128, k int, radius float64, n int) complex128 {
	if k < 0 {
		panic("complexanalysis: CauchyDerivative requires k >= 0")
	}
	g := func(z complex128) complex128 {
		return f(z) / cmplx.Pow(z-z0, complex(float64(k+1), 0))
	}
	integral := IntegrateCircle(g, center, radius, n) / complexanalysisTwoPiI
	return integral * complex(complexanalysisFactorialFloat(k), 0)
}

// WindingNumber returns the winding number (index) of the closed curve gamma
// about the point z0, computed by accumulating the change in the argument of
// gamma(t)-z0 over t in [0, 1] with n samples and dividing by 2*pi. The result
// is rounded to the nearest integer. z0 must not lie on the curve.
func WindingNumber(gamma func(float64) complex128, z0 complex128, n int) int {
	if n < 2 {
		n = 2
	}
	prev := gamma(0) - z0
	var total float64
	for j := 1; j <= n; j++ {
		cur := gamma(float64(j)/float64(n)) - z0
		total += cmplx.Phase(cur / prev) // change in argument in (-pi, pi]
		prev = cur
	}
	return int(math.Round(total / (2 * math.Pi)))
}

// ArgumentPrinciple returns the value of 1/(2*pi*i) * integral of f'(z)/f(z)
// around the positively oriented circle of the given center and radius. By the
// argument principle this equals the number of zeros minus the number of poles
// of f inside the circle, counted with multiplicity. The derivative f' is
// approximated by central finite differences. The exact count is returned by
// CountZeros after rounding.
func ArgumentPrinciple(f Function, center complex128, radius float64, n int) complex128 {
	const h = 1e-6
	logDeriv := func(z complex128) complex128 {
		fp := (f(z+complex(h, 0)) - f(z-complex(h, 0))) / complex(2*h, 0)
		return fp / f(z)
	}
	return IntegrateCircle(logDeriv, center, radius, n) / complexanalysisTwoPiI
}

// CountZeros returns the number of zeros of an analytic function f inside the
// circle of the given center and radius, counted with multiplicity, by rounding
// the argument-principle integral. It assumes f has no poles inside the circle
// and no zeros on it.
func CountZeros(f Function, center complex128, radius float64, n int) int {
	return int(math.Round(real(ArgumentPrinciple(f, center, radius, n))))
}
