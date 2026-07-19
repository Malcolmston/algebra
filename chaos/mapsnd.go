package chaos

import "math"

// MapN is an n-dimensional discrete map x_{k+1} = F(x_k) acting on Vec values.
// The returned vector must have the same dimension as the argument.
type MapN func(Vec) Vec

// Map2D is a two-dimensional map returning the next (x, y) from the current
// pair. It is a convenience form used by the planar canonical systems.
type Map2D func(x, y float64) (float64, float64)

// Lift2D turns a Map2D into a MapN acting on 2-vectors.
func Lift2D(f Map2D) MapN {
	return func(v Vec) Vec {
		x, y := f(v[0], v[1])
		return Vec{x, y}
	}
}

// HenonMap returns the Henon map (x,y) -> (1 - a*x^2 + y, b*x).
func HenonMap(a, b float64) Map2D {
	return func(x, y float64) (float64, float64) {
		return 1 - a*x*x + y, b * x
	}
}

// HenonJacobian returns the Jacobian matrix of the Henon map at (x, y).
func HenonJacobian(a, b, x, y float64) Mat {
	return Mat{
		{-2 * a * x, 1},
		{b, 0},
	}
}

// StandardMap returns the Chirikov standard map on the cylinder,
// p' = p + K sin(theta), theta' = theta + p', with theta reduced modulo 2pi.
func StandardMap(k float64) Map2D {
	return func(theta, p float64) (float64, float64) {
		pn := p + k*math.Sin(theta)
		tn := theta + pn
		tn = math.Mod(tn, 2*math.Pi)
		if tn < 0 {
			tn += 2 * math.Pi
		}
		return tn, pn
	}
}

// StandardMapJacobian returns the Jacobian of the standard map at (theta, p).
func StandardMapJacobian(k, theta, _ float64) Mat {
	c := k * math.Cos(theta)
	return Mat{
		{1 + c, 1},
		{c, 1},
	}
}

// IkedaMap returns the Ikeda map, a dissipative planar map with a spiral
// attractor, using the standard parameterisation with parameter u.
func IkedaMap(u float64) Map2D {
	return func(x, y float64) (float64, float64) {
		t := 0.4 - 6/(1+x*x+y*y)
		ct, st := math.Cos(t), math.Sin(t)
		return 1 + u*(x*ct-y*st), u * (x*st + y*ct)
	}
}

// TinkerbellMap returns the Tinkerbell map with parameters a, b, c, d.
func TinkerbellMap(a, b, c, d float64) Map2D {
	return func(x, y float64) (float64, float64) {
		return x*x - y*y + a*x + b*y, 2*x*y + c*x + d*y
	}
}

// GingerbreadmanMap returns the area-preserving Gingerbreadman map
// (x,y) -> (1 - y + |x|, x).
func GingerbreadmanMap() Map2D {
	return func(x, y float64) (float64, float64) {
		return 1 - y + math.Abs(x), x
	}
}

// IterateN applies the map F to v n times, returning the final vector.
func IterateN(F MapN, v Vec, n int) Vec {
	x := v.Clone()
	for i := 0; i < n; i++ {
		x = F(x)
	}
	return x
}

// OrbitN returns the sequence v, F(v), ..., F^n(v) as n+1 vectors.
func OrbitN(F MapN, v Vec, n int) []Vec {
	out := make([]Vec, n+1)
	out[0] = v.Clone()
	x := v.Clone()
	for i := 1; i <= n; i++ {
		x = F(x)
		out[i] = x.Clone()
	}
	return out
}

// OrbitNAfterTransient discards the first transient iterates and returns the
// next n points of the orbit.
func OrbitNAfterTransient(F MapN, v Vec, transient, n int) []Vec {
	x := v.Clone()
	for i := 0; i < transient; i++ {
		x = F(x)
	}
	out := make([]Vec, n)
	for i := 0; i < n; i++ {
		out[i] = x.Clone()
		x = F(x)
	}
	return out
}
