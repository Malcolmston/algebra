package physics

import "math"

// --- Ideal gas law, PV = nRT ---

// IdealGasVolume returns the volume of an ideal gas from V = n·R·T/P, for n
// moles at temperature T (K) and pressure P (Pa), where R is the molar gas
// constant. The result is in cubic metres (m³). P must be non-zero.
func IdealGasVolume(n, T, P float64) float64 { return n * GasConstant * T / P }

// IdealGasMoles returns the amount of substance of an ideal gas from
// n = P·V/(R·T), for pressure P (Pa) in volume V (m³) at temperature T (K),
// where R is the molar gas constant. The result is in moles (mol). T must be
// non-zero.
func IdealGasMoles(P, V, T float64) float64 { return P * V / (GasConstant * T) }

// IdealGasTemperature returns the temperature of an ideal gas from
// T = P·V/(n·R), for pressure P (Pa) in volume V (m³) with n moles, where R is
// the molar gas constant. The result is in kelvin (K). n must be non-zero.
func IdealGasTemperature(P, V, n float64) float64 { return P * V / (n * GasConstant) }

// RMSSpeed returns the root-mean-square molecular speed of an ideal gas,
// v_rms = √(3·R·T/M), for temperature T (K) and molar mass molarMass (kg/mol),
// where R is the molar gas constant. The result is in metres per second (m/s).
// molarMass must be non-zero.
func RMSSpeed(T, molarMass float64) float64 {
	return math.Sqrt(3 * GasConstant * T / molarMass)
}

// --- Heat transfer ---

// HeatEnergy returns the heat Q = m·c·ΔT required to change the temperature of
// a mass m (kg) with specific heat capacity c (J·kg⁻¹·K⁻¹) by deltaT (K). The
// result is in joules (J); it is negative when deltaT is negative (heat
// released).
func HeatEnergy(m, c, deltaT float64) float64 { return m * c * deltaT }

// ConductionRate returns the steady-state rate of heat conduction through a
// slab, Q/t = k·A·ΔT/d (Fourier's law), for thermal conductivity k
// (W·m⁻¹·K⁻¹), cross-sectional area area (m²), temperature difference deltaT
// (K) across the slab and thickness d (m). The result is in watts (W). d must
// be non-zero.
func ConductionRate(k, area, deltaT, d float64) float64 {
	return k * area * deltaT / d
}

// RadiatedPower returns the thermal power radiated by a grey body from the
// Stefan–Boltzmann law P = ε·σ·A·T⁴, for emissivity emissivity (dimensionless,
// 0–1), surface area area (m²) and absolute temperature T (K), where σ is the
// Stefan–Boltzmann constant. The result is in watts (W).
func RadiatedPower(emissivity, area, T float64) float64 {
	return emissivity * StefanBoltzmann * area * T * T * T * T
}

// ThermalExpansion returns the change in length ΔL = α·L0·ΔT of a linear
// object, for coefficient of linear expansion alpha (1/K), original length L0
// (m) and temperature change deltaT (K). The result is in metres (m).
func ThermalExpansion(alpha, L0, deltaT float64) float64 {
	return alpha * L0 * deltaT
}
