package physics

import (
	"errors"
	"math"
)

// This file broadens the electromagnetism coverage of the package beyond the
// scalar Coulomb/Ohm helpers in formulas.go: point-charge fields and potentials,
// energy densities, the magnetic field of a straight wire, magnetic and Lorentz
// forces, parallel-plate capacitance, stored energies, series/parallel network
// combinators for resistors and capacitors, and RC/LC transient helpers. All
// quantities are in SI units and every function is deterministic.
//
// Performance: the electric-field helpers reuse fourPiEps0, precomputed once at
// package initialisation, instead of recomputing 4·π·ε0 on every call, and the
// straight-wire field reuses mu0Over2Pi. The variadic network combinators
// accumulate their result in a single pass with no intermediate slices.

// fourPiEps0 is the constant 4·π·ε0 (F/m), precomputed once so the point-charge
// field and potential helpers avoid recomputing 4·π·ε0 on every call.
var fourPiEps0 = 4 * math.Pi * VacuumPermittivity

// mu0Over2Pi is the constant μ0/(2·π) (H/m over radians), precomputed once so
// MagneticFieldStraightWire avoids recomputing μ0/(2·π) on every call.
var mu0Over2Pi = VacuumPermeability / (2 * math.Pi)

// errEmptyNetwork is returned by the reciprocal network combinators when
// no element values are supplied.
var errEmptyNetwork = errors.New("physics: no element values provided")

// errZeroElement is returned by the reciprocal network combinators when
// an element value is zero, which would divide by zero.
var errZeroElement = errors.New("physics: zero-valued element in network")

// --- Fields, potentials, and energy density ---

// ElectricField returns the magnitude of the electric field of a point charge,
// E = q/(4·π·ε0·r²), for a charge q (C) at distance r (m). The result is in
// volts per metre (V/m); its sign follows the sign of q, so a positive charge
// yields a positive (radially outward) field. r must be non-zero.
func ElectricField(q, r float64) float64 {
	return q / (fourPiEps0 * r * r)
}

// ElectricPotential returns the electric potential of a point charge,
// V = q/(4·π·ε0·r), for a charge q (C) at distance r (m), taking the potential
// to be zero at infinity. The result is in volts (V); its sign follows the sign
// of q. r must be non-zero.
func ElectricPotential(q, r float64) float64 {
	return q / (fourPiEps0 * r)
}

// FieldEnergyDensity returns the energy density of an electric field in vacuum,
// u = ½·ε0·E², for a field magnitude E (V/m). The result is in joules per cubic
// metre (J/m³) and is always non-negative.
func FieldEnergyDensity(E float64) float64 {
	return 0.5 * VacuumPermittivity * E * E
}

// --- Magnetic fields and forces ---

// MagneticFieldStraightWire returns the magnitude of the magnetic field around
// a long straight wire, B = μ0·I/(2·π·r), for a current current (A) at
// perpendicular distance r (m). The result is in teslas (T). r must be
// non-zero.
func MagneticFieldStraightWire(current, r float64) float64 {
	return mu0Over2Pi * current / r
}

// MagneticForceOnCharge returns the magnitude of the magnetic force on a moving
// charge, F = q·v·B·sin(θ), for a charge q (C) moving at speed v (m/s) through a
// field of magnitude B (T), where theta (radians) is the angle between the
// velocity and the field. The result is in newtons (N). The sign follows the
// signs of q and sin(θ).
func MagneticForceOnCharge(q, v, B, theta float64) float64 {
	return q * v * B * math.Sin(theta)
}

// LorentzForceMag returns the combined electric and magnetic force on a charge
// as the scalar sum F = q·(E + v·B·sin(θ)), for a charge q (C) in an electric
// field of magnitude E (V/m) and a magnetic field of magnitude B (T), moving at
// speed v (m/s) at angle theta (radians) to the magnetic field. This assumes the
// electric force q·E and the magnetic force q·v·B·sin(θ) act along the same line
// (the collinear case); for non-collinear fields combine the vector components
// directly. The result is in newtons (N).
func LorentzForceMag(q, E, v, B, theta float64) float64 {
	return q * (E + v*B*math.Sin(theta))
}

// --- Capacitance, inductance, and stored energy ---

// ParallelPlateCapacitance returns the capacitance of a parallel-plate
// capacitor, C = ε0·εr·A/d, for plate area (m²), plate separation (m), and
// relative permittivity relPermittivity (dimensionless; 1 for vacuum). The
// result is in farads (F). separation must be non-zero.
func ParallelPlateCapacitance(area, separation, relPermittivity float64) float64 {
	return VacuumPermittivity * relPermittivity * area / separation
}

// InductorEnergy returns the energy stored in an inductor, E = ½·L·I², for an
// inductance l (H) carrying current i (A). The result is in joules (J).
func InductorEnergy(l, i float64) float64 {
	return 0.5 * l * i * i
}

// --- Circuit networks ---

// SeriesResistance returns the equivalent resistance of resistors connected in
// series, the sum of the individual resistances r (Ω), accumulated in a single
// pass. With no arguments it returns 0. The result is in ohms (Ω).
func SeriesResistance(r ...float64) float64 {
	var sum float64
	for _, v := range r {
		sum += v
	}
	return sum
}

// ParallelResistance returns the equivalent resistance of resistors connected in
// parallel, 1/Σ(1/rᵢ), for resistances r (Ω), accumulated in a single pass. It
// returns an error if no resistances are supplied or if any resistance is zero
// (a zero resistance short-circuits the network). The result is in ohms (Ω).
func ParallelResistance(r ...float64) (float64, error) {
	if len(r) == 0 {
		return 0, errEmptyNetwork
	}
	var recip float64
	for _, v := range r {
		if v == 0 {
			return 0, errZeroElement
		}
		recip += 1 / v
	}
	return 1 / recip, nil
}

// SeriesCapacitance returns the equivalent capacitance of capacitors connected
// in series, 1/Σ(1/cᵢ), for capacitances c (F), accumulated in a single pass. It
// returns an error if no capacitances are supplied or if any capacitance is
// zero. The result is in farads (F).
func SeriesCapacitance(c ...float64) (float64, error) {
	if len(c) == 0 {
		return 0, errEmptyNetwork
	}
	var recip float64
	for _, v := range c {
		if v == 0 {
			return 0, errZeroElement
		}
		recip += 1 / v
	}
	return 1 / recip, nil
}

// ParallelCapacitance returns the equivalent capacitance of capacitors connected
// in parallel, the sum of the individual capacitances c (F), accumulated in a
// single pass. With no arguments it returns 0. The result is in farads (F).
func ParallelCapacitance(c ...float64) float64 {
	var sum float64
	for _, v := range c {
		sum += v
	}
	return sum
}

// --- RC and LC transients ---

// RCTimeConstant returns the time constant τ = R·C of an RC circuit, for a
// resistance r (Ω) and capacitance c (F). The result is in seconds (s).
func RCTimeConstant(r, c float64) float64 {
	return r * c
}

// RCChargeVoltage returns the capacitor voltage while charging through a
// resistor, V = Vs·(1 − e^(−t/RC)), for a source voltage vSource (V) applied to
// a resistance r (Ω) and capacitance c (F) at time t (s) after the switch
// closes. The result is in volts (V). The product r·c must be non-zero.
func RCChargeVoltage(vSource, r, c, t float64) float64 {
	return vSource * (1 - math.Exp(-t/(r*c)))
}

// RCDischargeVoltage returns the capacitor voltage while discharging through a
// resistor, V = V0·e^(−t/RC), for an initial voltage v0 (V) across a resistance
// r (Ω) and capacitance c (F) at time t (s) after discharge begins. The result
// is in volts (V). The product r·c must be non-zero.
func RCDischargeVoltage(v0, r, c, t float64) float64 {
	return v0 * math.Exp(-t/(r*c))
}

// ResonantFrequency returns the resonant frequency of an LC circuit,
// f = 1/(2·π·√(L·C)), for an inductance l (H) and capacitance c (F). The result
// is in hertz (Hz). The product l·c must be positive.
func ResonantFrequency(l, c float64) float64 {
	return 1 / (2 * math.Pi * math.Sqrt(l*c))
}
