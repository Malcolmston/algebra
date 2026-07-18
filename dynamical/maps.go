package dynamical

import "math"

// Map1D is a real function of one real variable, used as the iteration rule of
// a one-dimensional discrete dynamical system x_{n+1} = f(x_n).
type Map1D func(float64) float64

// Point2D is a point in the plane, used as the state of a two-dimensional
// discrete dynamical system.
type Point2D struct {
	X, Y float64
}

// Add returns the vector sum p + q.
func (p Point2D) Add(q Point2D) Point2D { return Point2D{p.X + q.X, p.Y + q.Y} }

// Sub returns the vector difference p - q.
func (p Point2D) Sub(q Point2D) Point2D { return Point2D{p.X - q.X, p.Y - q.Y} }

// Scale returns the point p with both coordinates multiplied by s.
func (p Point2D) Scale(s float64) Point2D { return Point2D{p.X * s, p.Y * s} }

// Norm returns the Euclidean norm (length) of p treated as a vector.
func (p Point2D) Norm() float64 { return math.Hypot(p.X, p.Y) }

// Map2D is a map of the plane to itself, used as the iteration rule of a
// two-dimensional discrete dynamical system p_{n+1} = f(p_n).
type Map2D func(Point2D) Point2D

// Logistic evaluates the logistic map f(x) = r*x*(1-x) at x with parameter r.
// For r in [0,4] and x in [0,1] the unit interval is mapped into itself.
func Logistic(r, x float64) float64 { return r * x * (1 - x) }

// LogisticMap returns the logistic map x -> r*x*(1-x) as a [Map1D] closure for
// the fixed parameter r.
func LogisticMap(r float64) Map1D {
	return func(x float64) float64 { return r * x * (1 - x) }
}

// LogisticDeriv evaluates the derivative f'(x) = r*(1-2x) of the logistic map.
func LogisticDeriv(r, x float64) float64 { return r * (1 - 2*x) }

// Tent evaluates the tent map at x with slope parameter mu: the value is mu*x
// for x < 1/2 and mu*(1-x) for x >= 1/2.
func Tent(mu, x float64) float64 {
	if x < 0.5 {
		return mu * x
	}
	return mu * (1 - x)
}

// TentMap returns the tent map with slope parameter mu as a [Map1D] closure.
func TentMap(mu float64) Map1D {
	return func(x float64) float64 { return Tent(mu, x) }
}

// Doubling evaluates the dyadic doubling (Bernoulli) map f(x) = 2x mod 1.
func Doubling(x float64) float64 {
	v := 2 * x
	return v - math.Floor(v)
}

// DoublingMap returns the doubling map x -> 2x mod 1 as a [Map1D] closure.
func DoublingMap() Map1D { return func(x float64) float64 { return Doubling(x) } }

// SineMap evaluates the sine map f(x) = r*sin(pi*x) at x with parameter r.
func SineMap(r, x float64) float64 { return r * math.Sin(math.Pi*x) }

// GaussMap evaluates the Gauss (mouse) map f(x) = exp(-alpha*x^2) + beta.
func GaussMap(alpha, beta, x float64) float64 { return math.Exp(-alpha*x*x) + beta }

// CubicMap evaluates the cubic map f(x) = a*x - x^3 at x with parameter a.
func CubicMap(a, x float64) float64 { return a*x - x*x*x }

// QuadraticMap evaluates the quadratic map f(x) = x^2 + c, the canonical form
// to which every real quadratic map is conjugate.
func QuadraticMap(c, x float64) float64 { return x*x + c }

// CircleMap evaluates the standard circle map
// f(theta) = theta + omega - (k/2pi)*sin(2pi*theta), reduced modulo 1.
// With k = 0 it reduces to rigid rotation by omega.
func CircleMap(omega, k, theta float64) float64 {
	v := theta + omega - (k/(2*math.Pi))*math.Sin(2*math.Pi*theta)
	return v - math.Floor(v)
}

// CircleRotation returns the rigid rotation theta -> theta + omega mod 1 as a
// [Map1D] closure.
func CircleRotation(omega float64) Map1D {
	return func(theta float64) float64 { return CircleMap(omega, 0, theta) }
}

// Henon evaluates one step of the Henon map
// (x, y) -> (1 - a*x^2 + y, b*x) with parameters a and b.
func Henon(a, b, x, y float64) (float64, float64) {
	return 1 - a*x*x + y, b * x
}

// HenonMap returns the Henon map with parameters a and b as a [Map2D] closure.
func HenonMap(a, b float64) Map2D {
	return func(p Point2D) Point2D {
		nx, ny := Henon(a, b, p.X, p.Y)
		return Point2D{nx, ny}
	}
}

// HenonJacobian returns the 2x2 Jacobian matrix of the Henon map at point p,
// stored in row-major order: {{-2*a*x, 1}, {b, 0}}.
func HenonJacobian(a, b float64, p Point2D) [2][2]float64 {
	return [2][2]float64{{-2 * a * p.X, 1}, {b, 0}}
}
