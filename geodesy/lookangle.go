package geodesy

import "math"

// AER is a local look-angle: azimuth and elevation in degrees and slant range
// in metres, relative to an observer at a geodetic reference point.
type AER struct {
	// Azimuth is the compass bearing to the target in degrees, [0,360).
	Azimuth float64
	// Elevation is the angle above the local horizontal plane in degrees.
	Elevation float64
	// Range is the straight-line (slant) distance in metres.
	Range float64
}

// Dot returns the dot product of two ECEF vectors.
func (c ECEF) Dot(o ECEF) float64 { return c.X*o.X + c.Y*o.Y + c.Z*o.Z }

// Cross returns the cross product of two ECEF vectors.
func (c ECEF) Cross(o ECEF) ECEF {
	return ECEF{
		X: c.Y*o.Z - c.Z*o.Y,
		Y: c.Z*o.X - c.X*o.Z,
		Z: c.X*o.Y - c.Y*o.X,
	}
}

// Scale returns the ECEF vector scaled by s.
func (c ECEF) Scale(s float64) ECEF { return ECEF{c.X * s, c.Y * s, c.Z * s} }

// Normalize returns the unit vector in the direction of c (the zero vector is
// returned unchanged).
func (c ECEF) Normalize() ECEF {
	n := c.Norm()
	if n == 0 {
		return c
	}
	return c.Scale(1 / n)
}

// HorizontalRange returns the horizontal (ground) distance in metres of the ENU
// offset, ignoring the vertical component.
func (l ENU) HorizontalRange() float64 { return math.Hypot(l.E, l.N) }

// SlantRange returns the full 3-D distance in metres of the ENU offset.
func (l ENU) SlantRange() float64 { return math.Sqrt(l.E*l.E + l.N*l.N + l.U*l.U) }

// ENUToAER converts a local East-North-Up offset to azimuth/elevation/range
// look-angles.
func ENUToAER(l ENU) AER {
	az := NormalizeDegrees(deg(math.Atan2(l.E, l.N)))
	r := l.SlantRange()
	el := 0.0
	if r > 0 {
		el = deg(math.Asin(clamp(l.U/r, -1, 1)))
	}
	return AER{Azimuth: az, Elevation: el, Range: r}
}

// AERToENU converts azimuth/elevation/range look-angles to a local
// East-North-Up offset.
func AERToENU(a AER) ENU {
	az := rad(a.Azimuth)
	el := rad(a.Elevation)
	horiz := a.Range * math.Cos(el)
	return ENU{
		E: horiz * math.Sin(az),
		N: horiz * math.Cos(az),
		U: a.Range * math.Sin(el),
	}
}

// LookAngles returns the azimuth/elevation/range from an observer at the given
// geodetic reference point (degrees, metres height) to a target geodetic point,
// on the given ellipsoid.
func LookAngles(obsLat, obsLon, obsHeight, tgtLat, tgtLon, tgtHeight float64, e Ellipsoid) AER {
	enu := GeodeticToENU(tgtLat, tgtLon, tgtHeight, obsLat, obsLon, obsHeight, e)
	return ENUToAER(enu)
}
