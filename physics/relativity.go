package physics

import "math"

// physicsC2 is the square of the speed of light, precomputed once so the
// relativistic energy and momentum helpers avoid repeating the multiplication.
var physicsC2 = SpeedOfLight * SpeedOfLight

// Beta returns the dimensionless velocity ratio β = v/c for a speed v (m/s),
// where c is the speed of light. Values approach ±1 as |v| approaches c.
func Beta(v float64) float64 { return v / SpeedOfLight }

// RelativisticMomentum returns the relativistic momentum p = γ·m·v of a body of
// rest mass m (kg) moving at speed v (m/s), where γ is the Lorentz factor. The
// result is in kilogram-metres per second (kg·m/s). |v| must be less than c.
func RelativisticMomentum(m, v float64) float64 {
	return LorentzFactor(v) * m * v
}

// RelativisticEnergy returns the total relativistic energy E = γ·m·c² of a body
// of rest mass m (kg) moving at speed v (m/s), where γ is the Lorentz factor
// and c the speed of light. The result is in joules (J) and reduces to the rest
// energy m·c² at v = 0. |v| must be less than c.
func RelativisticEnergy(m, v float64) float64 {
	return LorentzFactor(v) * m * physicsC2
}

// RelativisticKineticEnergy returns the relativistic kinetic energy
// K = (γ − 1)·m·c² of a body of rest mass m (kg) moving at speed v (m/s), where
// γ is the Lorentz factor and c the speed of light. The result is in joules (J)
// and approaches the classical ½·m·v² for v ≪ c. |v| must be less than c.
func RelativisticKineticEnergy(m, v float64) float64 {
	return (LorentzFactor(v) - 1) * m * physicsC2
}

// EnergyMomentumRelation returns the total energy E of a body from the
// relativistic energy–momentum relation E = √((p·c)² + (m·c²)²), for rest mass
// m (kg) and momentum p (kg·m/s), where c is the speed of light. The result is
// in joules (J).
func EnergyMomentumRelation(m, p float64) float64 {
	pc := p * SpeedOfLight
	mc2 := m * physicsC2
	return math.Sqrt(pc*pc + mc2*mc2)
}

// RelativisticVelocityAddition returns the combined velocity of two collinear
// velocities u and v (m/s) under special relativity: w = (u + v)/(1 + u·v/c²),
// where c is the speed of light. The result never exceeds c in magnitude when
// the inputs do not. The result is in metres per second (m/s).
func RelativisticVelocityAddition(u, v float64) float64 {
	return (u + v) / (1 + u*v/physicsC2)
}

// TimeDilation returns the dilated time interval Δt = γ·properTime measured in a
// frame moving at speed v (m/s) relative to the clock, where properTime (s) is
// the interval in the clock's rest frame and γ the Lorentz factor. The result
// is in seconds (s). |v| must be less than c.
func TimeDilation(properTime, v float64) float64 {
	return LorentzFactor(v) * properTime
}

// LengthContraction returns the contracted length L = properLength/γ of an
// object measured in a frame in which it moves at speed v (m/s), where
// properLength (m) is its length in its rest frame and γ the Lorentz factor.
// The result is in metres (m). |v| must be less than c.
func LengthContraction(properLength, v float64) float64 {
	return properLength / LorentzFactor(v)
}

// RelativisticDopplerFactor returns the relativistic Doppler frequency ratio
// f_observed/f_source for motion directly along the line of sight at speed v
// (m/s). A positive v denotes the source and observer separating (redshift,
// factor < 1); a negative v denotes approach (blueshift, factor > 1). The
// factor is √((1 − β)/(1 + β)) with β = v/c. |v| must be less than c.
func RelativisticDopplerFactor(v float64) float64 {
	beta := v / SpeedOfLight
	return math.Sqrt((1 - beta) / (1 + beta))
}
