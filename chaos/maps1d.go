package chaos

import "math"

// Map1D is a one-dimensional real map x_{n+1} = f(x_n).
type Map1D func(float64) float64

// Logistic returns the logistic map x -> r*x*(1-x).
func Logistic(r float64) Map1D {
	return func(x float64) float64 { return r * x * (1 - x) }
}

// LogisticDeriv returns the derivative of the logistic map, x -> r*(1-2x).
func LogisticDeriv(r float64) Map1D {
	return func(x float64) float64 { return r * (1 - 2*x) }
}

// Tent returns the tent map with slope mu on [0,1]: 2mu*x for x<1/2 and
// 2mu*(1-x) otherwise.
func Tent(mu float64) Map1D {
	return func(x float64) float64 {
		if x < 0.5 {
			return 2 * mu * x
		}
		return 2 * mu * (1 - x)
	}
}

// TentDeriv returns the (piecewise constant) derivative of the tent map.
func TentDeriv(mu float64) Map1D {
	return func(x float64) float64 {
		if x < 0.5 {
			return 2 * mu
		}
		return -2 * mu
	}
}

// SineMap returns the sine map x -> a*sin(pi*x), a common smooth analogue of
// the logistic map.
func SineMap(a float64) Map1D {
	return func(x float64) float64 { return a * math.Sin(math.Pi*x) }
}

// SineMapDeriv returns the derivative of the sine map.
func SineMapDeriv(a float64) Map1D {
	return func(x float64) float64 { return a * math.Pi * math.Cos(math.Pi*x) }
}

// GaussMap returns the Gaussian (mouse) map x -> exp(-alpha*x^2) + beta.
func GaussMap(alpha, beta float64) Map1D {
	return func(x float64) float64 { return math.Exp(-alpha*x*x) + beta }
}

// GaussMapDeriv returns the derivative of the Gaussian map.
func GaussMapDeriv(alpha, beta float64) Map1D {
	return func(x float64) float64 { return -2 * alpha * x * math.Exp(-alpha*x*x) }
}

// CubicMap returns the cubic map x -> a*x^3 + (1-a)*x, a symmetric map with a
// period-doubling route to chaos.
func CubicMap(a float64) Map1D {
	return func(x float64) float64 { return a*x*x*x + (1-a)*x }
}

// CubicMapDeriv returns the derivative of the cubic map.
func CubicMapDeriv(a float64) Map1D {
	return func(x float64) float64 { return 3*a*x*x + (1 - a) }
}

// CircleMap returns the (uncoupled) circle map theta -> theta + Omega -
// (K/2pi) sin(2pi theta), reduced modulo 1.
func CircleMap(omega, k float64) Map1D {
	return func(x float64) float64 {
		y := x + omega - k/(2*math.Pi)*math.Sin(2*math.Pi*x)
		return y - math.Floor(y)
	}
}

// BernoulliMap returns the doubling (Bernoulli) map x -> (2x) mod 1.
func BernoulliMap() Map1D {
	return func(x float64) float64 {
		y := 2 * x
		return y - math.Floor(y)
	}
}

// Iterate applies f to x n times and returns the final value.
func Iterate(f Map1D, x float64, n int) float64 {
	for i := 0; i < n; i++ {
		x = f(x)
	}
	return x
}

// Orbit returns the sequence x, f(x), ..., f^n(x) of length n+1.
func Orbit(f Map1D, x float64, n int) []float64 {
	out := make([]float64, n+1)
	out[0] = x
	for i := 1; i <= n; i++ {
		x = f(x)
		out[i] = x
	}
	return out
}

// OrbitAfterTransient discards the first transient iterates and then returns
// the next n points of the orbit as a slice of length n.
func OrbitAfterTransient(f Map1D, x float64, transient, n int) []float64 {
	for i := 0; i < transient; i++ {
		x = f(x)
	}
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = x
		x = f(x)
	}
	return out
}

// Compose returns the n-fold composition f^n as a single Map1D.
func Compose(f Map1D, n int) Map1D {
	return func(x float64) float64 { return Iterate(f, x, n) }
}

// Cobweb returns the vertices of the cobweb (staircase) diagram of f started
// at x0 for n steps, as an alternating sequence of points on the diagonal and
// on the graph of f. The returned slices xs and ys have length 2n+1.
func Cobweb(f Map1D, x0 float64, n int) (xs, ys []float64) {
	xs = make([]float64, 0, 2*n+1)
	ys = make([]float64, 0, 2*n+1)
	x := x0
	xs = append(xs, x)
	ys = append(ys, 0)
	for i := 0; i < n; i++ {
		y := f(x)
		xs = append(xs, x)
		ys = append(ys, y)
		xs = append(xs, y)
		ys = append(ys, y)
		x = y
	}
	return xs, ys
}
