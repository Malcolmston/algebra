package physics

import (
	"math"
	"testing"
)

func TestProjectileMotion(t *testing.T) {
	// 45° launch at 10 m/s.
	assertApprox(t, "MaxHeight", ProjectileMaxHeight(10, math.Pi/4), 50/(2*StandardGravity), 1e-12)
	assertApprox(t, "TimeOfFlight", ProjectileTimeOfFlight(10, math.Pi/4),
		2*10*math.Sin(math.Pi/4)/StandardGravity, 1e-12)

	// At the apex time, vertical position equals the max height.
	tf := ProjectileTimeOfFlight(10, math.Pi/4)
	_, y := ProjectilePosition(10, math.Pi/4, tf/2)
	assertApprox(t, "apex y", y, ProjectileMaxHeight(10, math.Pi/4), 1e-9)
	// At landing time, y returns to zero.
	_, yEnd := ProjectilePosition(10, math.Pi/4, tf)
	assertApprox(t, "landing y", yEnd, 0, 1e-9)
}

func TestUniformAcceleration(t *testing.T) {
	assertApprox(t, "VelocityFinal", VelocityFinal(2, 3, 4), 14, 1e-12)
	assertApprox(t, "Displacement", Displacement(0, StandardGravity, 2), FreeFallDistance(2), 1e-12)
	assertApprox(t, "Centripetal", CentripetalAcceleration(10, 5), 20, 1e-12)
}

func TestOrbitalMechanics(t *testing.T) {
	const earthMass = 5.972e24 // kg
	const earthRadius = 6.371e6

	// Escape velocity from Earth's surface ≈ 11.19 km/s.
	assertApprox(t, "EscapeVelocity Earth", EscapeVelocity(earthMass, earthRadius), 11186, 2e-3)

	// Low Earth orbit at ~400 km altitude: v ≈ 7.67 km/s, T ≈ 92.5 min.
	rLEO := earthRadius + 400e3
	assertApprox(t, "OrbitalVelocity LEO", OrbitalVelocity(earthMass, rLEO), 7673, 2e-3)
	assertApprox(t, "OrbitalPeriod LEO", OrbitalPeriod(earthMass, rLEO)/60, 92.4, 5e-3)

	// Circular orbit relation: v = 2πr/T.
	v := OrbitalVelocity(earthMass, rLEO)
	tp := OrbitalPeriod(earthMass, rLEO)
	assertApprox(t, "v = 2πr/T", v, 2*math.Pi*rLEO/tp, 1e-9)

	// Earth's surface gravity from G·M/r² ≈ 9.8 m/s².
	assertApprox(t, "surface g", GravitationalForce(earthMass, 1, earthRadius), 9.82, 5e-3)
}

func TestCollisions1D(t *testing.T) {
	// Equal masses, elastic: velocities are exchanged.
	v1f, v2f := ElasticCollision1D(1, 1, 1, 0)
	assertApprox(t, "elastic v1f", v1f, 0, 1e-12)
	assertApprox(t, "elastic v2f", v2f, 1, 1e-12)

	// Momentum and kinetic energy conserved for arbitrary masses.
	v1f, v2f = ElasticCollision1D(2, 3, 5, -1)
	pBefore := Momentum(2, 3) + Momentum(5, -1)
	pAfter := Momentum(2, v1f) + Momentum(5, v2f)
	assertApprox(t, "elastic momentum", pAfter, pBefore, 1e-12)
	keBefore := KineticEnergy(2, 3) + KineticEnergy(5, -1)
	keAfter := KineticEnergy(2, v1f) + KineticEnergy(5, v2f)
	assertApprox(t, "elastic KE", keAfter, keBefore, 1e-12)

	// Perfectly inelastic: stuck bodies, momentum conserved.
	vf := InelasticCollision1D(1, 1, 1, 0)
	assertApprox(t, "inelastic vf", vf, 0.5, 1e-12)
	assertApprox(t, "inelastic momentum", Momentum(2, vf), Momentum(1, 1), 1e-12)
}

func BenchmarkOrbitalPeriod(b *testing.B) {
	var acc float64
	for i := 0; i < b.N; i++ {
		acc += OrbitalPeriod(5.972e24, 6.771e6)
	}
	_ = acc
}

func BenchmarkElasticCollision1D(b *testing.B) {
	var acc float64
	for i := 0; i < b.N; i++ {
		v1f, _ := ElasticCollision1D(2, 3, 5, -1)
		acc += v1f
	}
	_ = acc
}
