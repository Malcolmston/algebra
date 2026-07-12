package physics

import (
	"fmt"
	"math"
)

// A dimension groups units that measure the same physical quantity. Conversion
// is only defined between units of the same dimension.
type dimension int

const (
	dimLength dimension = iota
	dimMass
	dimTime
	dimTemperature
	dimEnergy
	dimAngle
)

// unitDef describes a linear unit: one unit equals factor base-SI units, so a
// value v in this unit equals v*factor in the base unit of its dimension.
// Temperature units are affine and are handled separately, not by factor.
type unitDef struct {
	dim    dimension
	factor float64
}

// units maps a unit symbol to its definition. Symbols are case-sensitive.
// Temperature symbols carry a factor of 1 and are converted through the affine
// helpers, never by multiplication.
var units = map[string]unitDef{
	// Length, base metre (m).
	"m":  {dimLength, 1},
	"km": {dimLength, 1000},
	"cm": {dimLength, 0.01},
	"mm": {dimLength, 0.001},
	"mi": {dimLength, 1609.344},
	"ft": {dimLength, 0.3048},
	"in": {dimLength, 0.0254},
	"nm": {dimLength, 1e-9},
	"Å":  {dimLength, 1e-10},

	// Mass, base kilogram (kg).
	"kg": {dimMass, 1},
	"g":  {dimMass, 0.001},
	"mg": {dimMass, 1e-6},
	"lb": {dimMass, 0.45359237},
	"oz": {dimMass, 0.028349523125},

	// Time, base second (s). One year is the Julian year of 365.25 days.
	"s":   {dimTime, 1},
	"min": {dimTime, 60},
	"hr":  {dimTime, 3600},
	"day": {dimTime, 86400},
	"yr":  {dimTime, 31557600},

	// Temperature, base kelvin (K); affine, converted via helpers.
	"K": {dimTemperature, 1},
	"C": {dimTemperature, 1},
	"F": {dimTemperature, 1},

	// Energy, base joule (J).
	"J":   {dimEnergy, 1},
	"eV":  {dimEnergy, ElectronVolt},
	"cal": {dimEnergy, 4.184},
	"kWh": {dimEnergy, 3.6e6},

	// Angle, base radian (rad).
	"rad": {dimAngle, 1},
	"deg": {dimAngle, math.Pi / 180},
}

// Convert converts value from the unit fromUnit to the unit toUnit and returns
// the result. Both symbols must be known and belong to the same physical
// dimension (for example "m" and "km", or "C" and "F"); otherwise Convert
// returns an error and a zero value. Temperature conversions are affine and
// handled correctly. Recognised symbols are:
//
//	length:      m, km, cm, mm, mi, ft, in, nm, Å
//	mass:        kg, g, mg, lb, oz
//	time:        s, min, hr, day, yr
//	temperature: K, C, F
//	energy:      J, eV, cal, kWh
//	angle:       rad, deg
func Convert(value float64, fromUnit, toUnit string) (float64, error) {
	from, ok := units[fromUnit]
	if !ok {
		return 0, fmt.Errorf("physics: unknown unit %q", fromUnit)
	}
	to, ok := units[toUnit]
	if !ok {
		return 0, fmt.Errorf("physics: unknown unit %q", toUnit)
	}
	if from.dim != to.dim {
		return 0, fmt.Errorf("physics: incompatible units %q and %q (different dimensions)", fromUnit, toUnit)
	}
	if from.dim == dimTemperature {
		return convertTemperature(value, fromUnit, toUnit), nil
	}
	// Linear: to base then to target.
	return value * from.factor / to.factor, nil
}

// convertTemperature converts an affine temperature between C, K and F. It
// assumes both symbols are valid temperature units.
func convertTemperature(value float64, fromUnit, toUnit string) float64 {
	// Convert the input to kelvin first.
	var k float64
	switch fromUnit {
	case "K":
		k = value
	case "C":
		k = value + 273.15
	case "F":
		k = (value-32)*5.0/9.0 + 273.15
	}
	// Convert kelvin to the target.
	switch toUnit {
	case "K":
		return k
	case "C":
		return k - 273.15
	case "F":
		return (k-273.15)*9.0/5.0 + 32
	}
	return k
}

// CelsiusToKelvin converts a temperature in degrees Celsius to kelvin.
func CelsiusToKelvin(c float64) float64 { return c + 273.15 }

// KelvinToCelsius converts a temperature in kelvin to degrees Celsius.
func KelvinToCelsius(k float64) float64 { return k - 273.15 }

// CelsiusToFahrenheit converts a temperature in degrees Celsius to degrees
// Fahrenheit.
func CelsiusToFahrenheit(c float64) float64 { return c*9.0/5.0 + 32 }

// FahrenheitToCelsius converts a temperature in degrees Fahrenheit to degrees
// Celsius.
func FahrenheitToCelsius(f float64) float64 { return (f - 32) * 5.0 / 9.0 }

// KelvinToFahrenheit converts a temperature in kelvin to degrees Fahrenheit.
func KelvinToFahrenheit(k float64) float64 { return (k-273.15)*9.0/5.0 + 32 }

// FahrenheitToKelvin converts a temperature in degrees Fahrenheit to kelvin.
func FahrenheitToKelvin(f float64) float64 { return (f-32)*5.0/9.0 + 273.15 }

// DegToRad converts an angle in degrees to radians.
func DegToRad(deg float64) float64 { return deg * math.Pi / 180 }

// RadToDeg converts an angle in radians to degrees.
func RadToDeg(rad float64) float64 { return rad * 180 / math.Pi }

// EVToJoules converts an energy in electronvolts to joules.
func EVToJoules(ev float64) float64 { return ev * ElectronVolt }

// JoulesToEV converts an energy in joules to electronvolts.
func JoulesToEV(j float64) float64 { return j / ElectronVolt }

// MilesToMeters converts a length in miles to metres.
func MilesToMeters(mi float64) float64 { return mi * 1609.344 }

// MetersToMiles converts a length in metres to miles.
func MetersToMiles(m float64) float64 { return m / 1609.344 }

// PoundsToKilograms converts a mass in pounds (avoirdupois) to kilograms.
func PoundsToKilograms(lb float64) float64 { return lb * 0.45359237 }

// KilogramsToPounds converts a mass in kilograms to pounds (avoirdupois).
func KilogramsToPounds(kg float64) float64 { return kg / 0.45359237 }
