// Package physics is the "science math" layer of the algebra module: a
// standard-library-only collection of physical constants, unit conversions and
// common physics formulas expressed in SI units.
//
// It has three parts:
//
//   - Physical constants. Exported float64 values in SI units, each documented
//     with its name, symbol, value and unit. Constants that are exact by
//     definition of the SI (2019 redefinition) are marked as such. The full
//     set is also available programmatically through [Constants] and [Lookup],
//     each entry being a [Constant] record — this powers generated
//     documentation.
//
//   - Unit conversions. A compact but useful set of conversions grouped by
//     physical dimension (length, mass, time, temperature, energy and angle).
//     Use the generic [Convert] for symbol-driven conversion — it returns an
//     error when the two units belong to different dimensions — or the typed
//     helpers such as [CelsiusToKelvin], [DegToRad] and [EVToJoules].
//
//   - Common formulas. Well-documented functions covering kinematics
//     ([KineticEnergy], [Force], [Momentum], [Work], [Power], [FreeFallDistance],
//     [ProjectileRange]), waves and optics ([WaveSpeed], [PhotonEnergy],
//     [DeBroglieWavelength]), relativity ([LorentzFactor], [MassEnergy]),
//     thermodynamics ([IdealGasPressure]) and electromagnetism ([CoulombForce],
//     [VoltageOhm], [CurrentOhm], [ResistanceOhm]).
//
// All quantities are SI unless a function's documentation states otherwise.
// The package depends only on the math standard-library package.
package physics
