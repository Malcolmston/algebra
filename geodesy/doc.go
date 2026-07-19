// Package geodesy implements geodetic computations and map projections on
// the WGS-84 reference ellipsoid (and other configurable ellipsoids) using
// only the Go standard library.
//
// The package covers the classical geodesy toolbox:
//
//   - Great-circle (spherical) calculations: Haversine, spherical law of
//     cosines and equirectangular distances, initial and final bearings,
//     midpoints, intermediate/interpolated points, destination points,
//     cross-track and along-track distances, path intersection, maximum
//     latitude of a great circle, and spherical polygon/triangle area.
//   - Rhumb-line (loxodrome) distance, bearing, midpoint and destination.
//   - Ellipsoidal geodesics via Vincenty's direct and inverse formulae,
//     including geodesic interpolation and waypoint generation.
//   - Map projections: spherical and Web (EPSG:3857) Mercator, the
//     Krüger-series transverse Mercator, and the Universal Transverse
//     Mercator (UTM) grid, with forward and inverse transforms.
//   - Military Grid Reference System (MGRS) encoding and decoding.
//   - Earth-Centred Earth-Fixed (ECEF) <-> geodetic conversions and local
//     East-North-Up (ENU) tangent-plane transforms.
//   - Angle helpers: degree/radian conversion, normalisation, and
//     degrees-minutes-seconds (DMS) parsing and formatting.
//
// All computation is performed in float64. Angles supplied to and returned
// from the public API are in degrees unless a function name explicitly says
// otherwise (Rad suffix), and distances are in metres. Functions are
// deterministic: identical inputs yield identical outputs.
package geodesy

import (
	"math"
)

// Mean Earth radii and reference constants (metres).
const (
	// EarthRadiusMean is the WGS-84 mean radius R1 = (2a+b)/3.
	EarthRadiusMean = 6371008.7714150598
	// EarthRadiusAuthalic is the radius of a sphere of equal surface area.
	EarthRadiusAuthalic = 6371007.1809184747
	// EarthRadiusVolumetric is the radius of a sphere of equal volume.
	EarthRadiusVolumetric = 6371000.7900091587
	// EarthRadiusEquatorial is the WGS-84 semi-major axis a.
	EarthRadiusEquatorial = 6378137.0
	// EarthRadiusPolar is the WGS-84 semi-minor axis b.
	EarthRadiusPolar = 6356752.314245179
)

// DegToRad converts an angle from degrees to radians.
func DegToRad(deg float64) float64 { return deg * math.Pi / 180 }

// RadToDeg converts an angle from radians to degrees.
func RadToDeg(rad float64) float64 { return rad * 180 / math.Pi }

// rad is the internal degrees->radians helper.
func rad(d float64) float64 { return d * math.Pi / 180 }

// deg is the internal radians->degrees helper.
func deg(r float64) float64 { return r * 180 / math.Pi }

// NormalizeDegrees reduces an angle to the half-open interval [0, 360).
func NormalizeDegrees(deg float64) float64 {
	d := math.Mod(deg, 360)
	if d < 0 {
		d += 360
	}
	return d
}

// NormalizeBearing reduces a compass bearing to the half-open interval
// [0, 360). It is an alias for NormalizeDegrees expressing intent.
func NormalizeBearing(bearing float64) float64 { return NormalizeDegrees(bearing) }

// NormalizeLongitude wraps a longitude into the half-open interval
// [-180, 180); the antimeridian is represented as -180.
func NormalizeLongitude(lon float64) float64 {
	l := math.Mod(lon+180, 360)
	if l < 0 {
		l += 360
	}
	return l - 180
}

// NormalizeLatitude clamps a latitude to the closed interval [-90, 90].
func NormalizeLatitude(lat float64) float64 {
	if lat > 90 {
		return 90
	}
	if lat < -90 {
		return -90
	}
	return lat
}

// NormalizeRadians reduces an angle in radians to the half-open interval
// [0, 2π).
func NormalizeRadians(r float64) float64 {
	x := math.Mod(r, 2*math.Pi)
	if x < 0 {
		x += 2 * math.Pi
	}
	return x
}

// WrapPi wraps an angle in radians to the half-open interval (-π, π].
func WrapPi(r float64) float64 {
	x := math.Mod(r+math.Pi, 2*math.Pi)
	if x <= 0 {
		x += 2 * math.Pi
	}
	return x - math.Pi
}

// WrapPi180 wraps an angle in degrees to the half-open interval (-180, 180].
func WrapPi180(d float64) float64 { return NormalizeLongitude(d) }

// DMSToDegrees converts a degrees/minutes/seconds triple to decimal degrees.
// The sign of the result follows the sign convention of deg (all components
// are treated with the magnitude of deg's sign); minutes and seconds must be
// non-negative magnitudes.
func DMSToDegrees(d, m, s float64) float64 {
	sign := 1.0
	if d < 0 {
		sign = -1
	}
	return sign * (math.Abs(d) + m/60 + s/3600)
}

// DegreesToDMS decomposes decimal degrees into a signed degree component and
// non-negative minutes and seconds. The sign is carried by the degrees
// return value (and, when degrees is zero, cannot be represented; callers
// needing the sign of a sub-degree value should inspect the input).
func DegreesToDMS(deg float64) (d int, m int, s float64) {
	sign := 1.0
	if deg < 0 {
		sign = -1
		deg = -deg
	}
	di := math.Floor(deg)
	rem := (deg - di) * 60
	mi := math.Floor(rem)
	sec := (rem - mi) * 60
	// Guard against rounding pushing seconds to 60.
	if sec >= 60 {
		sec -= 60
		mi++
	}
	if mi >= 60 {
		mi -= 60
		di++
	}
	return int(sign * di), int(mi), sec
}

// FormatDMS renders decimal degrees as a degrees/minutes/seconds string with
// the given hemisphere characters for positive and negative values (for
// example 'N'/'S' or 'E'/'W'). Seconds are shown with secPrec decimals.
func FormatDMS(deg float64, pos, neg byte, secPrec int) string {
	hemi := pos
	if deg < 0 {
		hemi = neg
	}
	d, m, s := DegreesToDMS(math.Abs(deg))
	// d is non-negative here because we passed the absolute value.
	return sprintfDMS(d, m, s, hemi, secPrec)
}
