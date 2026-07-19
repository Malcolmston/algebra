package geodesy

import "math"

// RhumbDistanceR returns the rhumb-line (loxodrome) distance in metres between
// two points on a sphere of the given radius. A rhumb line crosses every
// meridian at the same angle and appears as a straight line on a Mercator map.
func RhumbDistanceR(p1, p2 LatLon, radius float64) float64 {
	φ1, φ2 := rad(p1.Lat), rad(p2.Lat)
	dφ := φ2 - φ1
	dλ := rad(p2.Lon - p1.Lon)
	if math.Abs(dλ) > math.Pi {
		if dλ > 0 {
			dλ -= 2 * math.Pi
		} else {
			dλ += 2 * math.Pi
		}
	}
	dψ := math.Log(math.Tan(math.Pi/4+φ2/2) / math.Tan(math.Pi/4+φ1/2))
	var q float64
	if math.Abs(dψ) > 1e-12 {
		q = dφ / dψ
	} else {
		q = math.Cos(φ1)
	}
	return math.Hypot(dφ, q*dλ) * radius
}

// RhumbDistance is RhumbDistanceR on the WGS-84 mean sphere.
func RhumbDistance(p1, p2 LatLon) float64 {
	return RhumbDistanceR(p1, p2, EarthRadiusMean)
}

// RhumbBearing returns the constant compass bearing of the rhumb line from p1
// to p2, in degrees clockwise from north in [0,360).
func RhumbBearing(p1, p2 LatLon) float64 {
	φ1, φ2 := rad(p1.Lat), rad(p2.Lat)
	dλ := rad(p2.Lon - p1.Lon)
	if math.Abs(dλ) > math.Pi {
		if dλ > 0 {
			dλ -= 2 * math.Pi
		} else {
			dλ += 2 * math.Pi
		}
	}
	dψ := math.Log(math.Tan(math.Pi/4+φ2/2) / math.Tan(math.Pi/4+φ1/2))
	return NormalizeDegrees(deg(math.Atan2(dλ, dψ)))
}

// RhumbDestinationR returns the point reached by travelling distance metres
// along a rhumb line from start on the given constant bearing (degrees), on a
// sphere of the given radius.
func RhumbDestinationR(start LatLon, bearingDeg, distance, radius float64) LatLon {
	δ := distance / radius
	θ := rad(bearingDeg)
	φ1, λ1 := rad(start.Lat), rad(start.Lon)
	dφ := δ * math.Cos(θ)
	φ2 := φ1 + dφ
	// Clamp latitude to the poles.
	if math.Abs(φ2) > math.Pi/2 {
		if φ2 > 0 {
			φ2 = math.Pi / 2
		} else {
			φ2 = -math.Pi / 2
		}
	}
	dψ := math.Log(math.Tan(math.Pi/4+φ2/2) / math.Tan(math.Pi/4+φ1/2))
	var q float64
	if math.Abs(dψ) > 1e-12 {
		q = dφ / dψ
	} else {
		q = math.Cos(φ1)
	}
	dλ := δ * math.Sin(θ) / q
	λ2 := λ1 + dλ
	return LatLon{Lat: deg(φ2), Lon: NormalizeLongitude(deg(λ2))}
}

// RhumbDestination is RhumbDestinationR on the WGS-84 mean sphere.
func RhumbDestination(start LatLon, bearingDeg, distance float64) LatLon {
	return RhumbDestinationR(start, bearingDeg, distance, EarthRadiusMean)
}

// RhumbMidpoint returns the point half-way along the rhumb line between p1 and
// p2 (loxodromic midpoint).
func RhumbMidpoint(p1, p2 LatLon) LatLon {
	φ1, λ1 := rad(p1.Lat), rad(p1.Lon)
	φ2, λ2 := rad(p2.Lat), rad(p2.Lon)
	if math.Abs(λ2-λ1) > math.Pi {
		λ1 += 2 * math.Pi
	}
	φ3 := (φ1 + φ2) / 2
	f1 := math.Tan(math.Pi/4 + φ1/2)
	f2 := math.Tan(math.Pi/4 + φ2/2)
	f3 := math.Tan(math.Pi/4 + φ3/2)
	var λ3 float64
	den := math.Log(f2 / f1)
	if math.Abs(den) < 1e-12 {
		λ3 = (λ1 + λ2) / 2
	} else {
		λ3 = ((λ2-λ1)*math.Log(f3) + λ1*math.Log(f2) - λ2*math.Log(f1)) / den
	}
	return LatLon{Lat: deg(φ3), Lon: NormalizeLongitude(deg(λ3))}
}
