package geodesy

import "math"

// meridionalArc returns the meridian arc length (metres) from the equator to
// geodetic latitude φ (radians) on the given ellipsoid, using a series in the
// eccentricity accurate to well under a millimetre.
func meridionalArc(φ float64, e Ellipsoid) float64 {
	a := e.A
	e2 := e.FirstEccentricitySquared()
	e4 := e2 * e2
	e6 := e4 * e2
	e8 := e6 * e2
	c0 := 1 - e2/4 - 3*e4/64 - 5*e6/256 - 175*e8/16384
	c2 := 3.0/8*e2 + 3.0/32*e4 + 45.0/1024*e6 + 105.0/4096*e8
	c4 := 15.0/256*e4 + 45.0/1024*e6 + 525.0/16384*e8
	c6 := 35.0/3072*e6 + 175.0/12288*e8
	c8 := 315.0 / 131072 * e8
	return a * (c0*φ - c2*math.Sin(2*φ) + c4*math.Sin(4*φ) -
		c6*math.Sin(6*φ) + c8*math.Sin(8*φ))
}

// MeridianArcLength returns the meridian arc distance (metres) from the equator
// to the given geodetic latitude (degrees) on the given ellipsoid.
func MeridianArcLength(latDeg float64, e Ellipsoid) float64 {
	return meridionalArc(rad(latDeg), e)
}

// MeridianArcBetween returns the meridian arc distance (metres) between two
// latitudes (degrees) along the same meridian on the given ellipsoid.
func MeridianArcBetween(lat1, lat2 float64, e Ellipsoid) float64 {
	return math.Abs(meridionalArc(rad(lat2), e) - meridionalArc(rad(lat1), e))
}

// QuarterMeridian returns the length (metres) of a meridian quadrant (equator
// to pole) on the given ellipsoid.
func QuarterMeridian(e Ellipsoid) float64 { return meridionalArc(math.Pi/2, e) }

// LengthOfDegreeOfLatitude returns the north-south distance (metres) spanned by
// one degree of latitude at the given latitude (degrees) on the ellipsoid.
func LengthOfDegreeOfLatitude(latDeg float64, e Ellipsoid) float64 {
	return e.MeridianRadius(latDeg) * math.Pi / 180
}

// LengthOfDegreeOfLongitude returns the east-west distance (metres) spanned by
// one degree of longitude at the given latitude (degrees) on the ellipsoid.
func LengthOfDegreeOfLongitude(latDeg float64, e Ellipsoid) float64 {
	return e.PrimeVerticalRadius(latDeg) * math.Cos(rad(latDeg)) * math.Pi / 180
}

// RadiusOfParallel returns the radius (metres) of the circle of latitude (small
// circle) at the given latitude (degrees), N·cosφ, on the ellipsoid.
func RadiusOfParallel(latDeg float64, e Ellipsoid) float64 {
	return e.PrimeVerticalRadius(latDeg) * math.Cos(rad(latDeg))
}

// ParallelArcLength returns the distance (metres) along a parallel of latitude
// spanning dLon degrees of longitude at the given latitude on the ellipsoid.
func ParallelArcLength(latDeg, dLonDeg float64, e Ellipsoid) float64 {
	return RadiusOfParallel(latDeg, e) * rad(math.Abs(dLonDeg))
}

// ReducedLatitude returns the reduced (parametric) latitude β in degrees for a
// point on the ellipsoid surface at geodetic latitude latDeg.
func ReducedLatitude(latDeg float64, e Ellipsoid) float64 {
	return deg(math.Atan((1 - e.F) * math.Tan(rad(latDeg))))
}

// GeocentricLatitude returns the geocentric latitude (degrees) of a point on the
// ellipsoid surface at the given geodetic latitude.
func GeocentricLatitude(latDeg float64, e Ellipsoid) float64 {
	e2 := e.FirstEccentricitySquared()
	return deg(math.Atan((1 - e2) * math.Tan(rad(latDeg))))
}

// GeocentricToGeodeticLatitude returns the geodetic latitude (degrees) for a
// point on the ellipsoid surface at the given geocentric latitude.
func GeocentricToGeodeticLatitude(geocentricDeg float64, e Ellipsoid) float64 {
	e2 := e.FirstEccentricitySquared()
	return deg(math.Atan(math.Tan(rad(geocentricDeg)) / (1 - e2)))
}

// IsometricLatitude returns the isometric latitude ψ (dimensionless) at the
// given geodetic latitude (degrees) on the ellipsoid. It is the vertical
// coordinate (scaled by the equatorial radius) of the Mercator projection.
func IsometricLatitude(latDeg float64, e Ellipsoid) float64 {
	φ := rad(latDeg)
	ecc := e.FirstEccentricity()
	return math.Asinh(math.Tan(φ)) - ecc*math.Atanh(ecc*math.Sin(φ))
}

// ConformalLatitude returns the conformal latitude χ (degrees) at the given
// geodetic latitude on the ellipsoid.
func ConformalLatitude(latDeg float64, e Ellipsoid) float64 {
	ψ := IsometricLatitude(latDeg, e)
	return deg(2*math.Atan(math.Exp(ψ)) - math.Pi/2)
}

// RectifyingLatitude returns the rectifying latitude μ (degrees) at the given
// geodetic latitude, the latitude scaled so that meridian arc length is
// proportional to it.
func RectifyingLatitude(latDeg float64, e Ellipsoid) float64 {
	m := meridionalArc(rad(latDeg), e)
	mp := meridionalArc(math.Pi/2, e)
	return deg(math.Pi / 2 * m / mp)
}

// AuthalicLatitude returns the authalic (equal-area) latitude ξ (degrees) at the
// given geodetic latitude on the ellipsoid.
func AuthalicLatitude(latDeg float64, e Ellipsoid) float64 {
	if e.F == 0 {
		return latDeg
	}
	ecc := e.FirstEccentricity()
	q := func(φ float64) float64 {
		s := math.Sin(φ)
		return (1 - ecc*ecc) * (s/(1-ecc*ecc*s*s) - 1/(2*ecc)*math.Log((1-ecc*s)/(1+ecc*s)))
	}
	φ := rad(latDeg)
	qp := q(math.Pi / 2)
	ratio := q(φ) / qp
	ratio = clamp(ratio, -1, 1)
	return deg(math.Asin(ratio))
}
