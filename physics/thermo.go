package physics

import "math"

// This file expands the thermodynamics coverage of the physics package beyond
// the lone IdealGasPressure defined in formulas.go. Everything here is SI
// throughout and deterministic. Functions that divide by an argument document
// the non-zero requirement on that argument; passing zero yields ±Inf or NaN
// per IEEE 754 rather than a panic.

// physicsWienX is the dimensionless root of (x-5)·eˣ + 5 = 0, which arises when
// maximising the Planck spectral radiance with respect to wavelength. It is used
// to derive the Wien displacement constant from the base physical constants.
const physicsWienX = 4.965114231744276

// Package-level constants precomputed once from the base physical constants in
// constants.go, so the radiation formulas below never recompute these grouped
// products on every call.
var (
	// wienB is the Wien displacement-law constant b = h·c/(k_B·x) in metre-kelvin
	// (m·K), where x is physicsWienX. Its value is approximately 2.897771955e-3
	// m·K. WienPeakWavelength reuses it directly.
	wienB = PlanckConstant * SpeedOfLight / (BoltzmannConstant * physicsWienX)

	// physicsPlanckC1 is the grouped constant 2·h·c² (W·m²) that forms the
	// numerator prefactor of Planck's spectral radiance law.
	physicsPlanckC1 = 2 * PlanckConstant * SpeedOfLight * SpeedOfLight

	// physicsPlanckC2 is the grouped constant h·c/k_B (m·K) that forms the
	// exponent scale of Planck's spectral radiance law.
	physicsPlanckC2 = PlanckConstant * SpeedOfLight / BoltzmannConstant
)

// --- Ideal gas law, PV = nRT ---

// IdealGasVolume returns the volume of an ideal gas from V = n·R·T/P, for n
// moles at temperature T (K) and pressure P (Pa), where R is the molar gas
// constant. The result is in cubic metres (m³). P must be non-zero.
func IdealGasVolume(n, T, P float64) float64 { return n * GasConstant * T / P }

// IdealGasTemperature returns the temperature of an ideal gas from
// T = P·V/(n·R), for pressure P (Pa) in volume V (m³) with n moles, where R is
// the molar gas constant. The result is in kelvin (K). n must be non-zero.
func IdealGasTemperature(P, V, n float64) float64 { return P * V / (n * GasConstant) }

// IdealGasMoles returns the amount of substance of an ideal gas from
// n = P·V/(R·T), for pressure P (Pa) in volume V (m³) at temperature T (K),
// where R is the molar gas constant. The result is in moles (mol). T must be
// non-zero.
func IdealGasMoles(P, V, T float64) float64 { return P * V / (GasConstant * T) }

// NumberDensity returns the number of gas particles per unit volume from
// n = P/(k_B·T), for pressure P (Pa) at temperature T (K), where k_B is the
// Boltzmann constant. The result is in particles per cubic metre (m⁻³). T must
// be non-zero.
func NumberDensity(P, T float64) float64 { return P / (BoltzmannConstant * T) }

// --- Kinetic theory of gases ---

// RMSSpeed returns the root-mean-square molecular speed of an ideal gas,
// v_rms = √(3·R·T/M), for temperature T (K) and molar mass molarMass (kg/mol),
// where R is the molar gas constant. The result is in metres per second (m/s).
// molarMass must be non-zero.
func RMSSpeed(T, molarMass float64) float64 {
	return math.Sqrt(3 * GasConstant * T / molarMass)
}

// MeanKineticEnergy returns the mean translational kinetic energy of a gas
// particle, E = 3/2·k_B·T, at temperature T (K), where k_B is the Boltzmann
// constant. The result is in joules (J).
func MeanKineticEnergy(T float64) float64 { return 1.5 * BoltzmannConstant * T }

// MostProbableSpeed returns the most probable molecular speed of the
// Maxwell–Boltzmann distribution, v_p = √(2·R·T/M), for temperature T (K) and
// molar mass molarMass (kg/mol), where R is the molar gas constant. The result
// is in metres per second (m/s). molarMass must be non-zero.
func MostProbableSpeed(T, molarMass float64) float64 {
	return math.Sqrt(2 * GasConstant * T / molarMass)
}

// --- Heat transfer ---

// HeatCapacity returns the heat Q = m·c·ΔT required to change the temperature
// of a mass (kg) with specific heat capacity specificHeat (J·kg⁻¹·K⁻¹) by dT
// (K). The result is in joules (J); it is negative when dT is negative (heat
// released).
func HeatCapacity(mass, specificHeat, dT float64) float64 { return mass * specificHeat * dT }

// LatentHeat returns the heat Q = m·L absorbed or released during a phase
// change of a mass (kg) with specific latent heat L (J/kg). The result is in
// joules (J).
func LatentHeat(mass, L float64) float64 { return mass * L }

// ThermalConduction returns the steady-state rate of heat conduction through a
// slab, Q̇ = k·A·ΔT/L (Fourier's law), for thermal conductivity k
// (W·m⁻¹·K⁻¹), cross-sectional area area (m²), temperature difference dT (K)
// across the slab and thickness (m). The result is in watts (W). thickness must
// be non-zero.
func ThermalConduction(k, area, dT, thickness float64) float64 {
	return k * area * dT / thickness
}

// --- Cycles and efficiency ---

// CarnotEfficiency returns the maximum thermal efficiency 1 − Tc/Th of a heat
// engine operating between a cold reservoir at tCold (K) and a hot reservoir at
// tHot (K). The result is dimensionless (0–1 for tHot > tCold > 0). It performs
// no division guard, but tHot must be greater than zero for a physically
// meaningful result.
func CarnotEfficiency(tCold, tHot float64) float64 { return 1 - tCold/tHot }

// COPRefrigerator returns the ideal (Carnot) coefficient of performance of a
// refrigerator, COP = Tc/(Th − Tc), moving heat from a cold reservoir at tCold
// (K) to a hot reservoir at tHot (K). The result is dimensionless. tHot must
// differ from tCold.
func COPRefrigerator(tCold, tHot float64) float64 { return tCold / (tHot - tCold) }

// COPHeatPump returns the ideal (Carnot) coefficient of performance of a heat
// pump, COP = Th/(Th − Tc), delivering heat to a hot reservoir at tHot (K) while
// drawing from a cold reservoir at tCold (K). The result is dimensionless. tHot
// must differ from tCold.
func COPHeatPump(tCold, tHot float64) float64 { return tHot / (tHot - tCold) }

// --- Radiation and blackbody physics ---

// StefanBoltzmannPower returns the thermal power radiated by a grey body from
// the Stefan–Boltzmann law P = ε·σ·A·T⁴, for emissivity (dimensionless, 0–1),
// surface area area (m²) and absolute temperature T (K), where σ is the
// Stefan–Boltzmann constant. The result is in watts (W).
func StefanBoltzmannPower(emissivity, area, T float64) float64 {
	return emissivity * StefanBoltzmann * area * T * T * T * T
}

// WienPeakWavelength returns the wavelength λ_max = b/T at which a blackbody's
// spectral radiance peaks (Wien's displacement law), for absolute temperature T
// (K), where b is the Wien displacement constant. The result is in metres (m).
// T must be non-zero.
func WienPeakWavelength(T float64) float64 { return wienB / T }

// PlanckSpectralRadiance returns the spectral radiance of a blackbody per unit
// wavelength from Planck's law,
//
//	B(λ,T) = (2·h·c²/λ⁵) · 1/(exp(h·c/(λ·k_B·T)) − 1),
//
// for wavelength (m) and absolute temperature T (K), where h is the Planck
// constant, c the speed of light and k_B the Boltzmann constant. The result is
// in W·sr⁻¹·m⁻³ (watts per steradian per square metre per metre of wavelength).
// The grouped constants 2·h·c² and h·c/k_B are precomputed once at package
// initialisation, so this call performs no redundant constant products. Both
// wavelength and T must be non-zero.
func PlanckSpectralRadiance(wavelength, T float64) float64 {
	l2 := wavelength * wavelength
	l5 := l2 * l2 * wavelength
	return physicsPlanckC1 / (l5 * (math.Exp(physicsPlanckC2/(wavelength*T)) - 1))
}
