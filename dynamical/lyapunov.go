package dynamical

import "math"

// Lyapunov estimates the Lyapunov exponent of the one-dimensional map f from
// the analytic derivative df. Starting at x0 it discards the first transient
// iterates, then averages log|df(x)| over the next n iterates. The result is
// the mean exponential rate of separation of nearby trajectories; a positive
// value indicates sensitive dependence on initial conditions (chaos).
//
// For the tent map with slope mu the derivative has constant magnitude mu, so
// the estimate equals log(mu) exactly; for the logistic map at r = 4 it
// approaches log(2).
func Lyapunov(f, df Map1D, x0 float64, transient, n int) float64 {
	x := x0
	for i := 0; i < transient; i++ {
		x = f(x)
	}
	if n <= 0 {
		return 0
	}
	var sum float64
	for i := 0; i < n; i++ {
		d := math.Abs(df(x))
		if d == 0 {
			sum += math.Log(math.SmallestNonzeroFloat64)
		} else {
			sum += math.Log(d)
		}
		x = f(x)
	}
	return sum / float64(n)
}

// LyapunovSeparation estimates the Lyapunov exponent of the map f without an
// analytic derivative, using the trajectory-separation (renormalization)
// method. Two orbits are launched an initial distance d0 apart; after each
// step the separation is measured, its logarithmic growth accumulated, and the
// perturbed orbit rescaled back to distance d0. The first transient iterates
// are discarded before accumulation begins.
func LyapunovSeparation(f Map1D, x0, d0 float64, transient, n int) float64 {
	if d0 == 0 {
		d0 = 1e-8
	}
	x := x0
	for i := 0; i < transient; i++ {
		x = f(x)
	}
	if n <= 0 {
		return 0
	}
	xp := x + d0
	var sum float64
	for i := 0; i < n; i++ {
		x = f(x)
		xp = f(xp)
		d := math.Abs(xp - x)
		if d == 0 {
			d = d0 * math.SmallestNonzeroFloat64
		}
		sum += math.Log(d / d0)
		// Renormalize the perturbed orbit back to distance d0 from x.
		xp = x + (xp-x)*(d0/d)
	}
	return sum / float64(n)
}

// Lyapunov2D estimates the largest Lyapunov exponent of a two-dimensional map f
// whose Jacobian at any point is supplied by jac (row-major 2x2). A tangent
// vector is transported along the orbit by repeated multiplication with the
// Jacobian and renormalized to unit length each step; the accumulated
// logarithmic growth divided by n is returned. The first transient iterates of
// the base orbit are discarded first.
func Lyapunov2D(f Map2D, jac func(Point2D) [2][2]float64, p0 Point2D, transient, n int) float64 {
	p := p0
	for i := 0; i < transient; i++ {
		p = f(p)
	}
	if n <= 0 {
		return 0
	}
	vx, vy := 1.0, 0.0
	var sum float64
	for i := 0; i < n; i++ {
		j := jac(p)
		nx := j[0][0]*vx + j[0][1]*vy
		ny := j[1][0]*vx + j[1][1]*vy
		norm := math.Hypot(nx, ny)
		if norm == 0 {
			return math.Inf(-1)
		}
		sum += math.Log(norm)
		vx, vy = nx/norm, ny/norm
		p = f(p)
	}
	return sum / float64(n)
}

// LyapunovHenon estimates the largest Lyapunov exponent of the Henon map with
// parameters a and b by [Lyapunov2D], starting from the origin. For the
// classical parameters a = 1.4, b = 0.3 the exponent is approximately 0.419.
func LyapunovHenon(a, b float64, transient, n int) float64 {
	return Lyapunov2D(HenonMap(a, b), func(p Point2D) [2][2]float64 {
		return HenonJacobian(a, b, p)
	}, Point2D{0, 0}, transient, n)
}
