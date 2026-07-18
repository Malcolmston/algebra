package dynamical

import "math"

// Vec3 is a point or vector in three-dimensional space, used as the state of a
// continuous dynamical system.
type Vec3 struct {
	X, Y, Z float64
}

// Add returns the vector sum v + w.
func (v Vec3) Add(w Vec3) Vec3 { return Vec3{v.X + w.X, v.Y + w.Y, v.Z + w.Z} }

// Sub returns the vector difference v - w.
func (v Vec3) Sub(w Vec3) Vec3 { return Vec3{v.X - w.X, v.Y - w.Y, v.Z - w.Z} }

// Scale returns v with every component multiplied by s.
func (v Vec3) Scale(s float64) Vec3 { return Vec3{v.X * s, v.Y * s, v.Z * s} }

// Norm returns the Euclidean norm (length) of v.
func (v Vec3) Norm() float64 { return math.Sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z) }

// Dot returns the dot product v . w.
func (v Vec3) Dot(w Vec3) float64 { return v.X*w.X + v.Y*w.Y + v.Z*w.Z }

// Flow3D is the vector field of a continuous three-dimensional dynamical
// system: it maps a state s to its time derivative ds/dt = f(s).
type Flow3D func(Vec3) Vec3

// RK4Step3D advances the state s by one step of size dt under the vector field
// f using the classical fourth-order Runge-Kutta method, returning the new
// state. The local truncation error is O(dt^5) per step.
func RK4Step3D(f Flow3D, s Vec3, dt float64) Vec3 {
	k1 := f(s)
	k2 := f(s.Add(k1.Scale(dt / 2)))
	k3 := f(s.Add(k2.Scale(dt / 2)))
	k4 := f(s.Add(k3.Scale(dt)))
	incr := k1.Add(k2.Scale(2)).Add(k3.Scale(2)).Add(k4).Scale(dt / 6)
	return s.Add(incr)
}

// Integrate3D integrates the vector field f from the initial state s0 for n
// steps of size dt using [RK4Step3D], returning the trajectory as a slice of
// length n+1 whose first element is s0.
func Integrate3D(f Flow3D, s0 Vec3, dt float64, n int) []Vec3 {
	if n < 0 {
		n = 0
	}
	out := make([]Vec3, n+1)
	out[0] = s0
	s := s0
	for i := 1; i <= n; i++ {
		s = RK4Step3D(f, s, dt)
		out[i] = s
	}
	return out
}

// LorenzParams holds the three parameters of the Lorenz system.
type LorenzParams struct {
	Sigma, Rho, Beta float64
}

// DefaultLorenz returns the classical chaotic Lorenz parameters
// sigma = 10, rho = 28, beta = 8/3.
func DefaultLorenz() LorenzParams { return LorenzParams{Sigma: 10, Rho: 28, Beta: 8.0 / 3.0} }

// Field returns the Lorenz vector field
//
//	dx/dt = sigma*(y - x)
//	dy/dt = x*(rho - z) - y
//	dz/dt = x*y - beta*z
//
// as a [Flow3D] closure for these parameters.
func (p LorenzParams) Field() Flow3D {
	return func(s Vec3) Vec3 {
		return Vec3{
			X: p.Sigma * (s.Y - s.X),
			Y: s.X*(p.Rho-s.Z) - s.Y,
			Z: s.X*s.Y - p.Beta*s.Z,
		}
	}
}

// LorenzFixedPoints returns the equilibria of the Lorenz system with parameters
// p. The origin is always an equilibrium; for rho > 1 the two symmetric points
// C+ and C- at (+-sqrt(beta*(rho-1)), +-sqrt(beta*(rho-1)), rho-1) are included.
func (p LorenzParams) LorenzFixedPoints() []Vec3 {
	if p.Rho <= 1 {
		return []Vec3{{0, 0, 0}}
	}
	c := math.Sqrt(p.Beta * (p.Rho - 1))
	return []Vec3{
		{0, 0, 0},
		{c, c, p.Rho - 1},
		{-c, -c, p.Rho - 1},
	}
}

// LorenzTrajectory integrates the Lorenz system with parameters p from state s0
// for n RK4 steps of size dt, returning the trajectory (length n+1).
func LorenzTrajectory(p LorenzParams, s0 Vec3, dt float64, n int) []Vec3 {
	return Integrate3D(p.Field(), s0, dt, n)
}

// RosslerParams holds the three parameters of the Rossler system.
type RosslerParams struct {
	A, B, C float64
}

// DefaultRossler returns commonly used chaotic Rossler parameters
// a = 0.2, b = 0.2, c = 5.7.
func DefaultRossler() RosslerParams { return RosslerParams{A: 0.2, B: 0.2, C: 5.7} }

// Field returns the Rossler vector field
//
//	dx/dt = -y - z
//	dy/dt = x + a*y
//	dz/dt = b + z*(x - c)
//
// as a [Flow3D] closure for these parameters.
func (p RosslerParams) Field() Flow3D {
	return func(s Vec3) Vec3 {
		return Vec3{
			X: -s.Y - s.Z,
			Y: s.X + p.A*s.Y,
			Z: p.B + s.Z*(s.X-p.C),
		}
	}
}

// RosslerTrajectory integrates the Rossler system with parameters p from state
// s0 for n RK4 steps of size dt, returning the trajectory (length n+1).
func RosslerTrajectory(p RosslerParams, s0 Vec3, dt float64, n int) []Vec3 {
	return Integrate3D(p.Field(), s0, dt, n)
}
