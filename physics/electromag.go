package physics

// CoulombConstant is the Coulomb constant k_e = 1/(4·π·ε0), in newton square
// metres per coulomb squared (N·m²·C⁻²). It is precomputed once from the vacuum
// permittivity so the electrostatic helpers avoid recomputing 4·π·ε0 on every
// call. Its value is approximately 8.9875×10⁹.
var CoulombConstant = physicsCoulombK

// physicsCoulombK is the internal precomputed value backing CoulombConstant.
var physicsCoulombK = 1 / (4 * physicsPi * VacuumPermittivity)

// physicsPi is the value of π used by the electromagnetism helpers, kept local
// to avoid importing math for a single constant.
const physicsPi = 3.14159265358979323846

// ElectricFieldPointCharge returns the magnitude of the electric field
// E = k_e·q/r² produced by a point charge q (C) at distance r (m), where k_e is
// the Coulomb constant. The result is in volts per metre (V/m); its sign
// follows the sign of q (positive q points radially outward). r must be
// non-zero.
func ElectricFieldPointCharge(q, r float64) float64 {
	return physicsCoulombK * q / (r * r)
}

// ElectricPotentialPointCharge returns the electric potential V = k_e·q/r of a
// point charge q (C) at distance r (m), taking the potential to be zero at
// infinity, where k_e is the Coulomb constant. The result is in volts (V). r
// must be non-zero.
func ElectricPotentialPointCharge(q, r float64) float64 {
	return physicsCoulombK * q / r
}

// PowerDissipated returns the electrical power dissipated by a resistor,
// P = I²·R, carrying current i (A) with resistance r (Ω). The result is in
// watts (W).
func PowerDissipated(i, r float64) float64 { return i * i * r }

// ResistorsSeries returns the equivalent resistance of resistors connected in
// series, the sum of the individual resistances r (Ω). With no arguments it
// returns 0. The result is in ohms (Ω).
func ResistorsSeries(r ...float64) float64 {
	var sum float64
	for _, v := range r {
		sum += v
	}
	return sum
}

// ResistorsParallel returns the equivalent resistance of resistors connected in
// parallel, 1/Σ(1/rᵢ), for resistances r (Ω). With no arguments it returns 0. A
// zero resistance short-circuits the network and yields 0. The result is in
// ohms (Ω).
func ResistorsParallel(r ...float64) float64 {
	if len(r) == 0 {
		return 0
	}
	var recip float64
	for _, v := range r {
		if v == 0 {
			return 0
		}
		recip += 1 / v
	}
	return 1 / recip
}

// CapacitorEnergy returns the energy E = ½·C·V² stored in a capacitor of
// capacitance c (F) charged to voltage v (V). The result is in joules (J).
func CapacitorEnergy(c, v float64) float64 { return 0.5 * c * v * v }

// CapacitorCharge returns the charge Q = C·V stored on a capacitor of
// capacitance c (F) at voltage v (V). The result is in coulombs (C).
func CapacitorCharge(c, v float64) float64 { return c * v }

// CapacitorsSeries returns the equivalent capacitance of capacitors connected
// in series, 1/Σ(1/cᵢ), for capacitances c (F). With no arguments it returns 0.
// A zero capacitance yields 0. The result is in farads (F).
func CapacitorsSeries(c ...float64) float64 {
	if len(c) == 0 {
		return 0
	}
	var recip float64
	for _, v := range c {
		if v == 0 {
			return 0
		}
		recip += 1 / v
	}
	return 1 / recip
}

// CapacitorsParallel returns the equivalent capacitance of capacitors connected
// in parallel, the sum of the individual capacitances c (F). With no arguments
// it returns 0. The result is in farads (F).
func CapacitorsParallel(c ...float64) float64 {
	var sum float64
	for _, v := range c {
		sum += v
	}
	return sum
}
