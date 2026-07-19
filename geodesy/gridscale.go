package geodesy

import "math"

// GridConvergence returns the meridian convergence angle (degrees) of a
// transverse Mercator projection at the given geodetic point: the angle from
// grid north to true north, such that true azimuth = grid azimuth +
// convergence. It is computed numerically from the projection with the given
// parameters, so it is valid for any transverse Mercator grid.
func GridConvergence(lat, lon, lon0, k0, falseEasting, falseNorthing float64, e Ellipsoid) float64 {
	const d = 1e-6 // degrees
	e0, n0 := TransverseMercatorForward(lat, lon, lon0, k0, falseEasting, falseNorthing, e)
	e1, n1 := TransverseMercatorForward(lat+d, lon, lon0, k0, falseEasting, falseNorthing, e)
	// Grid azimuth of a step due true north (bearing 0).
	gridAz := math.Atan2(e1-e0, n1-n0)
	return -deg(gridAz)
}

// UTMGridConvergence returns the grid convergence angle (degrees) of the UTM
// projection at the given WGS-84 point.
func UTMGridConvergence(p LatLon) float64 {
	zone := UTMZone(p.Lat, p.Lon)
	lon0 := UTMCentralMeridian(zone)
	fn := 0.0
	if p.Lat < 0 {
		fn = UTMFalseNorthingSouth
	}
	return GridConvergence(p.Lat, p.Lon, lon0, UTMScaleFactor, UTMFalseEasting, fn, WGS84)
}

// PointScaleFactor returns the local scale factor of a transverse Mercator
// projection at the given geodetic point (the ratio of grid distance to true
// ellipsoidal distance), computed numerically for the given parameters.
func PointScaleFactor(lat, lon, lon0, k0, falseEasting, falseNorthing float64, e Ellipsoid) float64 {
	const d = 1e-6 // degrees east step
	e0, n0 := TransverseMercatorForward(lat, lon, lon0, k0, falseEasting, falseNorthing, e)
	e1, n1 := TransverseMercatorForward(lat, lon+d, lon0, k0, falseEasting, falseNorthing, e)
	gridDist := math.Hypot(e1-e0, n1-n0)
	trueDist := ParallelChordDistance(lat, d, e)
	if trueDist == 0 {
		return k0
	}
	return gridDist / trueDist
}

// ParallelChordDistance returns the true east-west distance (metres) along the
// parallel at latitude latDeg spanning dLon degrees of longitude on the
// ellipsoid. It is used internally by PointScaleFactor and is a thin alias of
// ParallelArcLength preserving sign-free magnitude.
func ParallelChordDistance(latDeg, dLonDeg float64, e Ellipsoid) float64 {
	return ParallelArcLength(latDeg, dLonDeg, e)
}

// UTMPointScale returns the local scale factor of the UTM projection at the
// given WGS-84 point.
func UTMPointScale(p LatLon) float64 {
	zone := UTMZone(p.Lat, p.Lon)
	lon0 := UTMCentralMeridian(zone)
	fn := 0.0
	if p.Lat < 0 {
		fn = UTMFalseNorthingSouth
	}
	return PointScaleFactor(p.Lat, p.Lon, lon0, UTMScaleFactor, UTMFalseEasting, fn, WGS84)
}
