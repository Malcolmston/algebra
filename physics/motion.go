package physics

import (
	"errors"
	"math"
)

// This file deepens the mechanics helpers that begin with the single-scalar
// formulas in formulas.go. Every routine is deterministic closed-form
// arithmetic: no iteration, no randomness, and no heap allocation apart from
// the two slice-consuming helpers, each of which makes a single pass over its
// input. The exported names here are additive and distinct from those in
// formulas.go.

// ErrCenterOfMass is returned by [CenterOfMass1D] when its inputs are unusable:
// the masses and positions slices have different lengths, the slices are empty,
// or the total mass is zero (so the centre of mass is undefined).
var ErrCenterOfMass = errors.New("physics: invalid masses/positions for centre of mass")

// --- Kinematics under constant acceleration ---

// Velocity returns the velocity v = v0 + a·t of a body with initial velocity
// v0 (m/s) under constant acceleration a (m/s²) after time t (s). The result is
// in metres per second (m/s). It is pure arithmetic and always well defined.
func Velocity(v0, a, t float64) float64 { return v0 + a*t }

// Displacement is deliberately not defined here; see the constant-acceleration
// displacement helper elsewhere in the package. This comment documents the
// intentional omission so the file's kinematics story reads completely.

// VelocityFromDistance returns the speed v = √(v0² + 2·a·d) reached by a body
// with initial velocity v0 (m/s) under constant acceleration a (m/s²) after
// travelling distance d (m), from the time-independent kinematic relation
// v² = v0² + 2·a·d. The result is in metres per second (m/s). The caller must
// ensure the radicand v0² + 2·a·d is non-negative; a negative radicand (for
// example strong deceleration over a distance the body cannot reach) yields
// NaN, signalling that the requested state is not physically attained.
func VelocityFromDistance(v0, a, d float64) float64 {
	return math.Sqrt(v0*v0 + 2*a*d)
}

// TimeToStop returns the time t = −v0/a for a body with initial velocity v0
// (m/s) under constant acceleration a (m/s²) to reach zero velocity. The result
// is in seconds (s). The acceleration a must be non-zero; with a = 0 the body
// never stops and the result is an infinity whose sign reflects v0. When a
// opposes v0 (genuine deceleration) the result is positive.
func TimeToStop(v0, a float64) float64 { return -v0 / a }

// --- Energy and momentum ---

// ImpulseMomentum returns the impulse J = force·dt delivered by a constant force
// (N) acting over a time interval dt (s). By the impulse–momentum theorem the
// impulse equals the change in linear momentum it produces, so the result is in
// kilogram-metres per second (kg·m/s), equivalently newton-seconds (N·s). It is
// pure arithmetic and always well defined; dt is normally positive.
func ImpulseMomentum(force, dt float64) float64 { return force * dt }

// CenterOfMass1D returns the mass-weighted mean position (Σ mᵢ·xᵢ)/(Σ mᵢ) of a
// system of point masses arranged on a line, where masses[i] (kg) is the mass
// at positions[i] (m). The result is in metres (m).
//
// It returns [ErrCenterOfMass] if masses and positions have different lengths,
// if they are empty, or if the total mass is zero (the centre of mass is then
// undefined). The two slices are read in a single pass and no memory is
// allocated.
func CenterOfMass1D(masses, positions []float64) (float64, error) {
	if len(masses) != len(positions) || len(masses) == 0 {
		return 0, ErrCenterOfMass
	}
	var totalMass, weighted float64
	for i, m := range masses {
		totalMass += m
		weighted += m * positions[i]
	}
	if totalMass == 0 {
		return 0, ErrCenterOfMass
	}
	return weighted / totalMass, nil
}

// --- Gravitation and orbital mechanics ---

// GravitationalPotentialEnergy returns the gravitational potential energy
// U = −G·m1·m2/r of two point masses m1 and m2 (kg) separated by distance r
// (m), where G is the Newtonian gravitational constant. The energy is measured
// relative to the conventional zero at infinite separation, so it is negative
// for the attractive gravitational interaction. The result is in joules (J).
// The separation r must be non-zero.
func GravitationalPotentialEnergy(m1, m2, r float64) float64 {
	return -GravitationalConstant * m1 * m2 / r
}

// VisViva returns the orbital speed v = √(G·M·(2/r − 1/a)) of a body at
// distance r (m) from a central mass mCentral (kg) on an orbit of semi-major
// axis a (m), from the vis-viva equation, where G is the Newtonian
// gravitational constant. The result is in metres per second (m/s). Both r and
// a must be non-zero; for a bound orbit r lies within the orbit so the radicand
// is non-negative, while a radicand driven negative by inconsistent inputs
// yields NaN.
func VisViva(mCentral, r, a float64) float64 {
	return math.Sqrt(GravitationalConstant * mCentral * (2/r - 1/a))
}
