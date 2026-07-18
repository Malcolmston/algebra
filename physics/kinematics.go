package physics

import "math"

// --- Projectile motion (constant standard gravity, no air resistance) ---

// ProjectileMaxHeight returns the maximum height reached by a projectile
// launched from ground level with initial speed v (m/s) at launch angle theta
// (radians), under standard gravity: H = (v·sinθ)²/(2·g). The result is in
// metres (m).
func ProjectileMaxHeight(v, theta float64) float64 {
	vy := v * math.Sin(theta)
	return vy * vy / (2 * StandardGravity)
}

// ProjectileTimeOfFlight returns the total time a projectile launched from and
// landing at the same height stays aloft, with initial speed v (m/s) at launch
// angle theta (radians), under standard gravity: t = 2·v·sinθ/g. The result is
// in seconds (s).
func ProjectileTimeOfFlight(v, theta float64) float64 {
	return 2 * v * math.Sin(theta) / StandardGravity
}

// ProjectilePosition returns the horizontal and vertical position (x, y), in
// metres, of a projectile at time t (s), launched from the origin with initial
// speed v (m/s) at launch angle theta (radians) under standard gravity:
// x = v·cosθ·t, y = v·sinθ·t − ½·g·t².
func ProjectilePosition(v, theta, t float64) (x, y float64) {
	x = v * math.Cos(theta) * t
	y = v*math.Sin(theta)*t - 0.5*StandardGravity*t*t
	return x, y
}

// --- Uniformly accelerated linear motion ---

// VelocityFinal returns the final velocity v = v0 + a·t of a body with initial
// velocity v0 (m/s) under constant acceleration a (m/s²) after time t (s). The
// result is in metres per second (m/s).
func VelocityFinal(v0, a, t float64) float64 { return v0 + a*t }

// Displacement returns the displacement s = v0·t + ½·a·t² of a body with
// initial velocity v0 (m/s) under constant acceleration a (m/s²) after time t
// (s). The result is in metres (m).
func Displacement(v0, a, t float64) float64 { return v0*t + 0.5*a*t*t }

// CentripetalAcceleration returns the centripetal acceleration a = v²/r
// required to keep a body moving at speed v (m/s) on a circle of radius r (m).
// The result is in metres per second squared (m/s²) and points toward the
// centre. r must be non-zero.
func CentripetalAcceleration(v, r float64) float64 { return v * v / r }

// --- Orbital mechanics (two-body, point masses) ---

// OrbitalVelocity returns the speed of a body in a circular orbit of radius r
// (m) about a central mass m (kg): v = √(G·m/r), where G is the Newtonian
// gravitational constant. The result is in metres per second (m/s). r must be
// non-zero.
func OrbitalVelocity(m, r float64) float64 {
	return math.Sqrt(GravitationalConstant * m / r)
}

// EscapeVelocity returns the escape speed from the surface of a body of mass m
// (kg) and radius r (m): v = √(2·G·m/r), where G is the Newtonian
// gravitational constant. The result is in metres per second (m/s). r must be
// non-zero.
func EscapeVelocity(m, r float64) float64 {
	return math.Sqrt(2 * GravitationalConstant * m / r)
}

// OrbitalPeriod returns the period of a circular orbit of radius r (m) about a
// central mass m (kg), from Kepler's third law: T = 2π·√(r³/(G·m)), where G is
// the Newtonian gravitational constant. The result is in seconds (s). m must be
// non-zero.
func OrbitalPeriod(m, r float64) float64 {
	return 2 * math.Pi * math.Sqrt(r*r*r/(GravitationalConstant*m))
}

// GravitationalForce returns the magnitude of the gravitational attraction
// F = G·m1·m2/r² between two point masses m1 and m2 (kg) separated by distance
// r (m), where G is the Newtonian gravitational constant. The result is in
// newtons (N). r must be non-zero.
func GravitationalForce(m1, m2, r float64) float64 {
	return GravitationalConstant * m1 * m2 / (r * r)
}

// --- Momentum and energy conservation in 1D collisions ---

// ElasticCollision1D returns the final velocities v1f and v2f (m/s) of two
// bodies after a one-dimensional perfectly elastic head-on collision, given
// masses m1 and m2 (kg) and initial velocities v1 and v2 (m/s). Both linear
// momentum and kinetic energy are conserved. The sum m1+m2 must be non-zero.
func ElasticCollision1D(m1, v1, m2, v2 float64) (v1f, v2f float64) {
	sum := m1 + m2
	v1f = ((m1-m2)*v1 + 2*m2*v2) / sum
	v2f = ((m2-m1)*v2 + 2*m1*v1) / sum
	return v1f, v2f
}

// InelasticCollision1D returns the common final velocity of two bodies that
// stick together after a one-dimensional perfectly inelastic collision, given
// masses m1 and m2 (kg) and initial velocities v1 and v2 (m/s):
// vf = (m1·v1 + m2·v2)/(m1 + m2). Momentum is conserved; kinetic energy is not.
// The result is in metres per second (m/s). The sum m1+m2 must be non-zero.
func InelasticCollision1D(m1, v1, m2, v2 float64) float64 {
	return (m1*v1 + m2*v2) / (m1 + m2)
}
