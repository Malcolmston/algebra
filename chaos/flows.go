package chaos

import "math"

// Field is the right-hand side of an autonomous ODE system dx/dt = f(x).
// It must return a vector of the same dimension as its argument.
type Field func(Vec) Vec

// Flow bundles a vector field with the phase-space dimension it acts on.
type Flow struct {
	// F is the vector field dx/dt = F(x).
	F Field
	// Dim is the dimension of the phase space.
	Dim int
}

// NewFlow constructs a Flow of the given dimension from a vector field.
func NewFlow(dim int, f Field) *Flow {
	return &Flow{F: f, Dim: dim}
}

// Eval returns the field value F(x).
func (fl *Flow) Eval(x Vec) Vec { return fl.F(x) }

// StepEuler advances the state x by one explicit-Euler step of size h.
func StepEuler(f Field, x Vec, h float64) Vec {
	return x.AddScaled(h, f(x))
}

// StepRK4 advances the state x by one classical fourth-order Runge-Kutta step
// of size h.
func StepRK4(f Field, x Vec, h float64) Vec {
	k1 := f(x)
	k2 := f(x.AddScaled(h/2, k1))
	k3 := f(x.AddScaled(h/2, k2))
	k4 := f(x.AddScaled(h, k3))
	next := x.Clone()
	for i := range next {
		next[i] += h / 6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
	}
	return next
}

// StepRK2 advances the state x by one midpoint (second-order Runge-Kutta)
// step of size h.
func StepRK2(f Field, x Vec, h float64) Vec {
	k1 := f(x)
	k2 := f(x.AddScaled(h/2, k1))
	return x.AddScaled(h, k2)
}

// Integrate advances x for n steps of size h with RK4 and returns the final
// state.
func Integrate(f Field, x Vec, h float64, n int) Vec {
	s := x.Clone()
	for i := 0; i < n; i++ {
		s = StepRK4(f, s, h)
	}
	return s
}

// Trajectory integrates x for n steps of size h with RK4 and returns the
// n+1 states x(0), x(h), ..., x(nh).
func Trajectory(f Field, x Vec, h float64, n int) []Vec {
	out := make([]Vec, n+1)
	out[0] = x.Clone()
	s := x.Clone()
	for i := 1; i <= n; i++ {
		s = StepRK4(f, s, h)
		out[i] = s.Clone()
	}
	return out
}

// TrajectoryAfterTransient integrates away the first transient steps and then
// returns the next n+1 recorded states.
func TrajectoryAfterTransient(f Field, x Vec, h float64, transient, n int) []Vec {
	s := x.Clone()
	for i := 0; i < transient; i++ {
		s = StepRK4(f, s, h)
	}
	return Trajectory(f, s, h, n)
}

// Times returns the sample times 0, h, 2h, ..., n*h.
func Times(h float64, n int) []float64 {
	t := make([]float64, n+1)
	for i := 0; i <= n; i++ {
		t[i] = float64(i) * h
	}
	return t
}

// Lorenz returns the Lorenz vector field with parameters sigma, rho, beta.
func Lorenz(sigma, rho, beta float64) Field {
	return func(v Vec) Vec {
		x, y, z := v[0], v[1], v[2]
		return Vec{
			sigma * (y - x),
			x*(rho-z) - y,
			x*y - beta*z,
		}
	}
}

// LorenzStandard returns the Lorenz field with the classic chaotic parameters
// sigma=10, rho=28, beta=8/3.
func LorenzStandard() Field {
	return Lorenz(10, 28, 8.0/3.0)
}

// LorenzJacobian returns the Jacobian of the Lorenz field at state v.
func LorenzJacobian(sigma, rho, beta float64, v Vec) Mat {
	x, y, z := v[0], v[1], v[2]
	return Mat{
		{-sigma, sigma, 0},
		{rho - z, -1, -x},
		{y, x, -beta},
	}
}

// Rossler returns the Rossler vector field with parameters a, b, c.
func Rossler(a, b, c float64) Field {
	return func(v Vec) Vec {
		x, y, z := v[0], v[1], v[2]
		return Vec{
			-y - z,
			x + a*y,
			b + z*(x-c),
		}
	}
}

// RosslerStandard returns the Rossler field with the common chaotic
// parameters a=0.2, b=0.2, c=5.7.
func RosslerStandard() Field {
	return Rossler(0.2, 0.2, 5.7)
}

// RosslerJacobian returns the Jacobian of the Rossler field at state v.
func RosslerJacobian(a, _, c float64, v Vec) Mat {
	x, _, z := v[0], v[1], v[2]
	return Mat{
		{0, -1, -1},
		{1, a, 0},
		{z, 0, x - c},
	}
}

// DampedPendulum returns the field of a periodically forced, damped pendulum
// written as an autonomous 3-D system (theta, omega, phase) with damping q,
// forcing amplitude A and drive frequency wd:
//
//	dtheta/dt = omega
//	domega/dt = -sin(theta) - q*omega + A*cos(phi)
//	dphi/dt   = wd
func DampedPendulum(q, a, wd float64) Field {
	return func(v Vec) Vec {
		theta, omega, phi := v[0], v[1], v[2]
		return Vec{
			omega,
			-math.Sin(theta) - q*omega + a*math.Cos(phi),
			wd,
		}
	}
}

// DuffingField returns the field of the forced Duffing oscillator
// x” + delta x' - x + x^3 = gamma cos(omega t) as an autonomous 3-D system
// in (x, v, phase).
func DuffingField(delta, gamma, omega float64) Field {
	return func(u Vec) Vec {
		x, v, phi := u[0], u[1], u[2]
		return Vec{
			v,
			-delta*v + x - x*x*x + gamma*math.Cos(phi),
			omega,
		}
	}
}

// VanDerPolField returns the field of the (unforced) Van der Pol oscillator
// x” - mu(1-x^2)x' + x = 0 as a planar system in (x, v).
func VanDerPolField(mu float64) Field {
	return func(u Vec) Vec {
		x, v := u[0], u[1]
		return Vec{v, mu*(1-x*x)*v - x}
	}
}
