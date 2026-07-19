package geodesy

import "math"

// ECEF is a point in the Earth-Centred, Earth-Fixed Cartesian frame, with the
// origin at the centre of mass, the Z axis toward the reference pole, the X
// axis toward the intersection of the equator and the prime meridian, and the
// Y axis completing a right-handed system. Units are metres.
type ECEF struct {
	X, Y, Z float64
}

// ENU is a local tangent-plane coordinate: East, North and Up offsets in metres
// relative to a reference geodetic origin.
type ENU struct {
	E, N, U float64
}

// GeodeticToECEF converts geodetic latitude/longitude (degrees) and ellipsoidal
// height (metres) to ECEF coordinates on the given ellipsoid.
func GeodeticToECEF(lat, lon, height float64, e Ellipsoid) ECEF {
	φ := rad(lat)
	λ := rad(lon)
	es := e.FirstEccentricitySquared()
	sinφ := math.Sin(φ)
	cosφ := math.Cos(φ)
	N := e.A / math.Sqrt(1-es*sinφ*sinφ)
	x := (N + height) * cosφ * math.Cos(λ)
	y := (N + height) * cosφ * math.Sin(λ)
	z := (N*(1-es) + height) * sinφ
	return ECEF{X: x, Y: y, Z: z}
}

// ECEFToGeodetic converts ECEF coordinates to geodetic latitude/longitude
// (degrees) and ellipsoidal height (metres) on the given ellipsoid, using
// Bowring's method refined by one iteration for sub-millimetre accuracy.
func ECEFToGeodetic(c ECEF, e Ellipsoid) (lat, lon, height float64) {
	a := e.A
	b := e.SemiMinorAxis()
	es := e.FirstEccentricitySquared()
	eps := e.SecondEccentricitySquared()

	p := math.Hypot(c.X, c.Y)
	λ := math.Atan2(c.Y, c.X)

	if p < 1e-12 {
		// On the polar axis.
		lat = 90
		if c.Z < 0 {
			lat = -90
		}
		return lat, deg(λ), math.Abs(c.Z) - b
	}

	θ := math.Atan2(c.Z*a, p*b)
	sinθ := math.Sin(θ)
	cosθ := math.Cos(θ)
	φ := math.Atan2(c.Z+eps*b*sinθ*sinθ*sinθ, p-es*a*cosθ*cosθ*cosθ)

	// One Newton-style refinement.
	for i := 0; i < 2; i++ {
		sinφ := math.Sin(φ)
		N := a / math.Sqrt(1-es*sinφ*sinφ)
		h := p/math.Cos(φ) - N
		φ = math.Atan2(c.Z, p*(1-es*N/(N+h)))
	}
	sinφ := math.Sin(φ)
	N := a / math.Sqrt(1-es*sinφ*sinφ)
	h := p/math.Cos(φ) - N
	return deg(φ), deg(λ), h
}

// Add returns the vector sum of two ECEF points.
func (c ECEF) Add(o ECEF) ECEF { return ECEF{c.X + o.X, c.Y + o.Y, c.Z + o.Z} }

// Sub returns the vector difference c - o.
func (c ECEF) Sub(o ECEF) ECEF { return ECEF{c.X - o.X, c.Y - o.Y, c.Z - o.Z} }

// Norm returns the Euclidean length of the ECEF vector.
func (c ECEF) Norm() float64 { return math.Sqrt(c.X*c.X + c.Y*c.Y + c.Z*c.Z) }

// Distance returns the straight-line (chord) distance between two ECEF points.
func (c ECEF) Distance(o ECEF) float64 { return c.Sub(o).Norm() }

// ChordDistance returns the straight-line distance in metres through the Earth
// between two geodetic surface points on the given ellipsoid (height 0).
func ChordDistance(p1, p2 LatLon, e Ellipsoid) float64 {
	a := GeodeticToECEF(p1.Lat, p1.Lon, 0, e)
	b := GeodeticToECEF(p2.Lat, p2.Lon, 0, e)
	return a.Distance(b)
}

// ECEFToENU converts an ECEF point to local East-North-Up coordinates relative
// to a geodetic reference origin (degrees, metres) on the given ellipsoid.
func ECEFToENU(c ECEF, refLat, refLon, refHeight float64, e Ellipsoid) ENU {
	origin := GeodeticToECEF(refLat, refLon, refHeight, e)
	d := c.Sub(origin)
	φ := rad(refLat)
	λ := rad(refLon)
	sinφ, cosφ := math.Sin(φ), math.Cos(φ)
	sinλ, cosλ := math.Sin(λ), math.Cos(λ)
	east := -sinλ*d.X + cosλ*d.Y
	north := -sinφ*cosλ*d.X - sinφ*sinλ*d.Y + cosφ*d.Z
	up := cosφ*cosλ*d.X + cosφ*sinλ*d.Y + sinφ*d.Z
	return ENU{E: east, N: north, U: up}
}

// ENUToECEF converts a local East-North-Up offset back to an ECEF point,
// relative to a geodetic reference origin on the given ellipsoid.
func ENUToECEF(l ENU, refLat, refLon, refHeight float64, e Ellipsoid) ECEF {
	origin := GeodeticToECEF(refLat, refLon, refHeight, e)
	φ := rad(refLat)
	λ := rad(refLon)
	sinφ, cosφ := math.Sin(φ), math.Cos(φ)
	sinλ, cosλ := math.Sin(λ), math.Cos(λ)
	dx := -sinλ*l.E - sinφ*cosλ*l.N + cosφ*cosλ*l.U
	dy := cosλ*l.E - sinφ*sinλ*l.N + cosφ*sinλ*l.U
	dz := cosφ*l.N + sinφ*l.U
	return ECEF{X: origin.X + dx, Y: origin.Y + dy, Z: origin.Z + dz}
}

// GeodeticToENU converts a geodetic point (degrees, metres) to local
// East-North-Up coordinates relative to a geodetic reference origin.
func GeodeticToENU(lat, lon, height, refLat, refLon, refHeight float64, e Ellipsoid) ENU {
	return ECEFToENU(GeodeticToECEF(lat, lon, height, e), refLat, refLon, refHeight, e)
}

// ENUToGeodetic converts local East-North-Up coordinates back to a geodetic
// position (degrees latitude/longitude, metres height).
func ENUToGeodetic(l ENU, refLat, refLon, refHeight float64, e Ellipsoid) (lat, lon, height float64) {
	c := ENUToECEF(l, refLat, refLon, refHeight, e)
	return ECEFToGeodetic(c, e)
}
