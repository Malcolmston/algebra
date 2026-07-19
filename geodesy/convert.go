package geodesy

import (
	"math"
	"strings"
)

// Length unit conversion factors (metres per unit).
const (
	// MetersPerNauticalMile is the international nautical mile in metres.
	MetersPerNauticalMile = 1852.0
	// MetersPerStatuteMile is the statute (land) mile in metres.
	MetersPerStatuteMile = 1609.344
	// MetersPerFoot is the international foot in metres.
	MetersPerFoot = 0.3048
	// MetersPerKilometer is a kilometre in metres.
	MetersPerKilometer = 1000.0
)

// MetersToKilometers converts metres to kilometres.
func MetersToKilometers(m float64) float64 { return m / MetersPerKilometer }

// KilometersToMeters converts kilometres to metres.
func KilometersToMeters(km float64) float64 { return km * MetersPerKilometer }

// MetersToNauticalMiles converts metres to nautical miles.
func MetersToNauticalMiles(m float64) float64 { return m / MetersPerNauticalMile }

// NauticalMilesToMeters converts nautical miles to metres.
func NauticalMilesToMeters(nm float64) float64 { return nm * MetersPerNauticalMile }

// MetersToMiles converts metres to statute miles.
func MetersToMiles(m float64) float64 { return m / MetersPerStatuteMile }

// MilesToMeters converts statute miles to metres.
func MilesToMeters(mi float64) float64 { return mi * MetersPerStatuteMile }

// MetersToFeet converts metres to feet.
func MetersToFeet(m float64) float64 { return m / MetersPerFoot }

// FeetToMeters converts feet to metres.
func FeetToMeters(ft float64) float64 { return ft * MetersPerFoot }

// DegToGrad converts degrees to gradians (400 gradians = 360 degrees).
func DegToGrad(deg float64) float64 { return deg * 10.0 / 9.0 }

// GradToDeg converts gradians to degrees.
func GradToDeg(grad float64) float64 { return grad * 9.0 / 10.0 }

// RadToGrad converts radians to gradians.
func RadToGrad(r float64) float64 { return r * 200 / math.Pi }

// GradToRad converts gradians to radians.
func GradToRad(g float64) float64 { return g * math.Pi / 200 }

// BackAzimuth returns the reverse bearing (degrees, [0,360)) of a forward
// bearing.
func BackAzimuth(bearing float64) float64 { return NormalizeDegrees(bearing + 180) }

// RelativeBearing returns the signed angle (degrees, (-180,180]) from bearing
// from to bearing to, that is the turn required to go from one heading to the
// other (positive clockwise).
func RelativeBearing(from, to float64) float64 {
	return WrapPi180(to - from)
}

// compass32 holds the 32-point compass rose abbreviations in clockwise order
// starting at north.
var compass32 = []string{
	"N", "NbE", "NNE", "NEbN", "NE", "NEbE", "ENE", "EbN",
	"E", "EbS", "ESE", "SEbE", "SE", "SEbS", "SSE", "SbE",
	"S", "SbW", "SSW", "SWbS", "SW", "SWbW", "WSW", "WbS",
	"W", "WbN", "WNW", "NWbW", "NW", "NWbN", "NNW", "NbW",
}

// compassPoint returns the compass abbreviation for a bearing using the given
// number of points (must be 4, 8, 16 or 32).
func compassPoint(bearing float64, points int) string {
	step := 32 / points
	b := NormalizeDegrees(bearing)
	idx := int(math.Round(b/(360.0/float64(points)))) % points
	return compass32[idx*step]
}

// BearingToCompass8 returns the 8-point compass abbreviation (N, NE, E, ...)
// nearest to the given bearing (degrees).
func BearingToCompass8(bearing float64) string { return compassPoint(bearing, 8) }

// BearingToCompass16 returns the 16-point compass abbreviation (N, NNE, NE,
// ...) nearest to the given bearing (degrees).
func BearingToCompass16(bearing float64) string { return compassPoint(bearing, 16) }

// BearingToCompass32 returns the 32-point compass abbreviation nearest to the
// given bearing (degrees).
func BearingToCompass32(bearing float64) string { return compassPoint(bearing, 32) }

// CompassToBearing returns the central bearing (degrees) of the named compass
// point (case-insensitive, for example "NE" or "SSW"). It returns
// ErrInvalidCoordinate if the name is not a recognised 32-point abbreviation.
func CompassToBearing(name string) (float64, error) {
	name = strings.TrimSpace(name)
	for i, p := range compass32 {
		if strings.EqualFold(p, name) {
			return float64(i) * (360.0 / 32.0), nil
		}
	}
	return 0, ErrInvalidCoordinate
}
