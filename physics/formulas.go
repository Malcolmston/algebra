package physics

import "math"

// --- Kinematics and mechanics ---

// KineticEnergy returns the translational kinetic energy ½·m·v² of a body of
// mass m (kg) moving at speed v (m/s). The result is in joules (J).
func KineticEnergy(m, v float64) float64 { return 0.5 * m * v * v }

// PotentialEnergyGravity returns the gravitational potential energy m·g·h of a
// body of mass m (kg) at height h (m) near Earth's surface, using the standard
// gravity g. The result is in joules (J).
func PotentialEnergyGravity(m, h float64) float64 { return m * StandardGravity * h }

// Force returns the net force m·a on a body of mass m (kg) undergoing
// acceleration a (m/s²), per Newton's second law. The result is in newtons (N).
func Force(m, a float64) float64 { return m * a }

// Momentum returns the linear momentum m·v of a body of mass m (kg) moving at
// velocity v (m/s). The result is in kilogram-metres per second (kg·m/s).
func Momentum(m, v float64) float64 { return m * v }

// Work returns the work done by a constant force f (N) acting over a
// displacement d (m) parallel to the force: W = f·d. The result is in joules
// (J). For a force at angle θ to the displacement, pass f·cos(θ).
func Work(f, d float64) float64 { return f * d }

// Power returns the average power work/time for work (J) done in time t (s).
// The result is in watts (W). t must be non-zero.
func Power(work, t float64) float64 { return work / t }

// FreeFallDistance returns the distance ½·g·t² fallen from rest under standard
// gravity after time t (s), ignoring air resistance. The result is in metres
// (m).
func FreeFallDistance(t float64) float64 { return 0.5 * StandardGravity * t * t }

// ProjectileRange returns the horizontal range of a projectile launched from
// and landing at the same height, with initial speed v (m/s) at launch angle
// theta (radians), under standard gravity and ignoring air resistance:
// R = v²·sin(2θ)/g. The result is in metres (m).
func ProjectileRange(v, theta float64) float64 {
	return v * v * math.Sin(2*theta) / StandardGravity
}

// --- Waves and optics ---

// WaveSpeed returns the propagation speed of a wave, freq·wavelength, for a
// wave of frequency freq (Hz) and wavelength wavelength (m). The result is in
// metres per second (m/s).
func WaveSpeed(freq, wavelength float64) float64 { return freq * wavelength }

// PhotonEnergy returns the energy h·f of a photon of frequency freq (Hz), where
// h is the Planck constant. The result is in joules (J).
func PhotonEnergy(freq float64) float64 { return PlanckConstant * freq }

// DeBroglieWavelength returns the de Broglie wavelength h/(m·v) of a particle of
// mass m (kg) moving at speed v (m/s), where h is the Planck constant. The
// result is in metres (m). m and v must be non-zero.
func DeBroglieWavelength(m, v float64) float64 { return PlanckConstant / (m * v) }

// --- Relativity ---

// LorentzFactor returns the Lorentz factor γ = 1/√(1 − v²/c²) for a speed v
// (m/s), where c is the speed of light. It is dimensionless and equals 1 at
// v = 0 and grows without bound as v approaches c. |v| must be less than c.
func LorentzFactor(v float64) float64 {
	beta := v / SpeedOfLight
	return 1 / math.Sqrt(1-beta*beta)
}

// MassEnergy returns the rest energy m·c² of a mass m (kg), where c is the
// speed of light. The result is in joules (J).
func MassEnergy(m float64) float64 { return m * SpeedOfLight * SpeedOfLight }

// --- Thermodynamics ---

// IdealGasPressure returns the pressure of an ideal gas from the ideal gas law
// P = n·R·T/V, for n moles at temperature T (K) in volume V (m³), where R is
// the molar gas constant. The result is in pascals (Pa). V must be non-zero.
func IdealGasPressure(n, T, V float64) float64 { return n * GasConstant * T / V }

// --- Electromagnetism ---

// CoulombForce returns the magnitude and sign of the electrostatic force
// between two point charges q1 and q2 (C) separated by distance r (m):
// F = q1·q2/(4·π·ε0·r²). A positive result indicates repulsion (like charges),
// a negative result attraction (opposite charges). The result is in newtons
// (N). r must be non-zero.
func CoulombForce(q1, q2, r float64) float64 {
	return q1 * q2 / (4 * math.Pi * VacuumPermittivity * r * r)
}

// VoltageOhm returns the voltage V = I·R across a resistor carrying current i
// (A) with resistance r (Ω), per Ohm's law. The result is in volts (V).
func VoltageOhm(i, r float64) float64 { return i * r }

// CurrentOhm returns the current I = V/R through a resistor of resistance r (Ω)
// with voltage v (V) across it, per Ohm's law. The result is in amperes (A).
// r must be non-zero.
func CurrentOhm(v, r float64) float64 { return v / r }

// ResistanceOhm returns the resistance R = V/I of a component with voltage v
// (V) across it carrying current i (A), per Ohm's law. The result is in ohms
// (Ω). i must be non-zero.
func ResistanceOhm(v, i float64) float64 { return v / i }
