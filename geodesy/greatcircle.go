package geodesy

import "math"

// AngularDistance returns the great-circle angular separation between two
// points, in radians, using the Haversine formula.
func AngularDistance(p1, p2 LatLon) float64 {
	φ1, φ2 := rad(p1.Lat), rad(p2.Lat)
	dφ := rad(p2.Lat - p1.Lat)
	dλ := rad(p2.Lon - p1.Lon)
	a := math.Sin(dφ/2)*math.Sin(dφ/2) +
		math.Cos(φ1)*math.Cos(φ2)*math.Sin(dλ/2)*math.Sin(dλ/2)
	return 2 * math.Atan2(math.Sqrt(a), math.Sqrt(1-a))
}

// HaversineDistanceR returns the great-circle distance between two points on a
// sphere of the given radius, using the Haversine formula.
func HaversineDistanceR(p1, p2 LatLon, radius float64) float64 {
	return AngularDistance(p1, p2) * radius
}

// HaversineDistance returns the great-circle distance in metres on a sphere of
// the WGS-84 mean radius.
func HaversineDistance(p1, p2 LatLon) float64 {
	return HaversineDistanceR(p1, p2, EarthRadiusMean)
}

// SphericalLawOfCosinesDistance returns the great-circle distance in metres on
// a sphere of the given radius, using the spherical law of cosines. It is
// simpler than the Haversine formula but less accurate for very short
// distances.
func SphericalLawOfCosinesDistance(p1, p2 LatLon, radius float64) float64 {
	φ1, φ2 := rad(p1.Lat), rad(p2.Lat)
	dλ := rad(p2.Lon - p1.Lon)
	c := math.Sin(φ1)*math.Sin(φ2) + math.Cos(φ1)*math.Cos(φ2)*math.Cos(dλ)
	if c > 1 {
		c = 1
	} else if c < -1 {
		c = -1
	}
	return math.Acos(c) * radius
}

// EquirectangularDistance returns an approximate distance in metres using the
// equirectangular (flat-Earth) projection. It is fast and adequate for small
// separations but degrades over long distances and near the poles.
func EquirectangularDistance(p1, p2 LatLon, radius float64) float64 {
	φ1, φ2 := rad(p1.Lat), rad(p2.Lat)
	dλ := rad(p2.Lon - p1.Lon)
	x := dλ * math.Cos((φ1+φ2)/2)
	y := φ2 - φ1
	return math.Hypot(x, y) * radius
}

// InitialBearingRad returns the initial (forward azimuth) bearing from p1 to p2
// along the great circle, in radians measured clockwise from north in [0,2π).
func InitialBearingRad(p1, p2 LatLon) float64 {
	φ1, φ2 := rad(p1.Lat), rad(p2.Lat)
	dλ := rad(p2.Lon - p1.Lon)
	y := math.Sin(dλ) * math.Cos(φ2)
	x := math.Cos(φ1)*math.Sin(φ2) - math.Sin(φ1)*math.Cos(φ2)*math.Cos(dλ)
	return NormalizeRadians(math.Atan2(y, x))
}

// InitialBearing returns the initial bearing from p1 to p2 along the great
// circle, in degrees clockwise from north in [0,360).
func InitialBearing(p1, p2 LatLon) float64 {
	return NormalizeDegrees(deg(InitialBearingRad(p1, p2)))
}

// FinalBearing returns the bearing at which the great circle from p1 arrives at
// p2, in degrees clockwise from north in [0,360).
func FinalBearing(p1, p2 LatLon) float64 {
	return NormalizeDegrees(InitialBearing(p2, p1) + 180)
}

// Midpoint returns the half-way point along the great circle between p1 and p2.
func Midpoint(p1, p2 LatLon) LatLon {
	φ1, λ1 := rad(p1.Lat), rad(p1.Lon)
	φ2 := rad(p2.Lat)
	dλ := rad(p2.Lon - p1.Lon)
	bx := math.Cos(φ2) * math.Cos(dλ)
	by := math.Cos(φ2) * math.Sin(dλ)
	φm := math.Atan2(math.Sin(φ1)+math.Sin(φ2),
		math.Sqrt((math.Cos(φ1)+bx)*(math.Cos(φ1)+bx)+by*by))
	λm := λ1 + math.Atan2(by, math.Cos(φ1)+bx)
	return LatLon{Lat: deg(φm), Lon: NormalizeLongitude(deg(λm))}
}

// IntermediatePoint returns the point at the given fraction of the great-circle
// distance from p1 (fraction 0) to p2 (fraction 1).
func IntermediatePoint(p1, p2 LatLon, fraction float64) LatLon {
	δ := AngularDistance(p1, p2)
	φ1, λ1 := rad(p1.Lat), rad(p1.Lon)
	φ2, λ2 := rad(p2.Lat), rad(p2.Lon)
	if δ == 0 {
		return p1
	}
	sinδ := math.Sin(δ)
	a := math.Sin((1-fraction)*δ) / sinδ
	b := math.Sin(fraction*δ) / sinδ
	x := a*math.Cos(φ1)*math.Cos(λ1) + b*math.Cos(φ2)*math.Cos(λ2)
	y := a*math.Cos(φ1)*math.Sin(λ1) + b*math.Cos(φ2)*math.Sin(λ2)
	z := a*math.Sin(φ1) + b*math.Sin(φ2)
	φi := math.Atan2(z, math.Hypot(x, y))
	λi := math.Atan2(y, x)
	return LatLon{Lat: deg(φi), Lon: NormalizeLongitude(deg(λi))}
}

// InterpolatePoints returns n+1 points evenly spaced (by fraction of the total
// great-circle distance) from p1 to p2 inclusive. n must be >= 1.
func InterpolatePoints(p1, p2 LatLon, n int) []LatLon {
	if n < 1 {
		n = 1
	}
	pts := make([]LatLon, 0, n+1)
	for i := 0; i <= n; i++ {
		pts = append(pts, IntermediatePoint(p1, p2, float64(i)/float64(n)))
	}
	return pts
}

// DestinationR returns the point reached by travelling distance metres along a
// great circle from the start point on the given initial bearing (degrees), on
// a sphere of the given radius.
func DestinationR(start LatLon, bearingDeg, distance, radius float64) LatLon {
	δ := distance / radius
	θ := rad(bearingDeg)
	φ1, λ1 := rad(start.Lat), rad(start.Lon)
	sinφ2 := math.Sin(φ1)*math.Cos(δ) + math.Cos(φ1)*math.Sin(δ)*math.Cos(θ)
	φ2 := math.Asin(sinφ2)
	y := math.Sin(θ) * math.Sin(δ) * math.Cos(φ1)
	x := math.Cos(δ) - math.Sin(φ1)*sinφ2
	λ2 := λ1 + math.Atan2(y, x)
	return LatLon{Lat: deg(φ2), Lon: NormalizeLongitude(deg(λ2))}
}

// Destination returns the point reached from start along a great circle on the
// given initial bearing (degrees) for distance metres, on the WGS-84 mean
// sphere.
func Destination(start LatLon, bearingDeg, distance float64) LatLon {
	return DestinationR(start, bearingDeg, distance, EarthRadiusMean)
}

// CrossTrackDistanceR returns the signed distance in metres of point p from the
// great circle defined by the start and end points, on a sphere of the given
// radius. The result is positive if p lies to the left of the start->end path
// and negative if to the right.
func CrossTrackDistanceR(p, start, end LatLon, radius float64) float64 {
	δ13 := AngularDistance(start, p)
	θ13 := InitialBearingRad(start, p)
	θ12 := InitialBearingRad(start, end)
	return math.Asin(math.Sin(δ13)*math.Sin(θ13-θ12)) * radius
}

// CrossTrackDistance is CrossTrackDistanceR on the WGS-84 mean sphere.
func CrossTrackDistance(p, start, end LatLon) float64 {
	return CrossTrackDistanceR(p, start, end, EarthRadiusMean)
}

// AlongTrackDistanceR returns the distance in metres from the start point to the
// closest point on the great circle start->end to the foot of the
// perpendicular from p, on a sphere of the given radius.
func AlongTrackDistanceR(p, start, end LatLon, radius float64) float64 {
	δ13 := AngularDistance(start, p)
	θ13 := InitialBearingRad(start, p)
	θ12 := InitialBearingRad(start, end)
	δxt := math.Asin(math.Sin(δ13) * math.Sin(θ13-θ12))
	c := math.Cos(δ13) / math.Cos(δxt)
	if c > 1 {
		c = 1
	} else if c < -1 {
		c = -1
	}
	δat := math.Acos(c)
	if math.Cos(θ12-θ13) < 0 {
		δat = -δat
	}
	return δat * radius
}

// AlongTrackDistance is AlongTrackDistanceR on the WGS-84 mean sphere.
func AlongTrackDistance(p, start, end LatLon) float64 {
	return AlongTrackDistanceR(p, start, end, EarthRadiusMean)
}

// Intersection returns the point of intersection of two great-circle paths,
// each defined by a point and an initial bearing (degrees). It returns
// ErrNoConvergence if the paths are (anti)parallel and do not intersect
// uniquely.
func Intersection(p1 LatLon, bearing1 float64, p2 LatLon, bearing2 float64) (LatLon, error) {
	φ1, λ1 := rad(p1.Lat), rad(p1.Lon)
	φ2, λ2 := rad(p2.Lat), rad(p2.Lon)
	θ13 := rad(bearing1)
	θ23 := rad(bearing2)
	dφ := φ2 - φ1
	dλ := λ2 - λ1

	δ12 := 2 * math.Asin(math.Sqrt(math.Sin(dφ/2)*math.Sin(dφ/2)+
		math.Cos(φ1)*math.Cos(φ2)*math.Sin(dλ/2)*math.Sin(dλ/2)))
	if math.Abs(δ12) < 1e-15 {
		return p1, nil // coincident points
	}

	cosθa := (math.Sin(φ2) - math.Sin(φ1)*math.Cos(δ12)) / (math.Sin(δ12) * math.Cos(φ1))
	cosθb := (math.Sin(φ1) - math.Sin(φ2)*math.Cos(δ12)) / (math.Sin(δ12) * math.Cos(φ2))
	cosθa = clamp(cosθa, -1, 1)
	cosθb = clamp(cosθb, -1, 1)
	θa := math.Acos(cosθa)
	θb := math.Acos(cosθb)

	var θ12, θ21 float64
	if math.Sin(dλ) > 0 {
		θ12 = θa
		θ21 = 2*math.Pi - θb
	} else {
		θ12 = 2*math.Pi - θa
		θ21 = θb
	}

	α1 := θ13 - θ12
	α2 := θ21 - θ23
	if math.Sin(α1) == 0 && math.Sin(α2) == 0 {
		return LatLon{}, ErrNoConvergence // infinite intersections
	}
	if math.Sin(α1)*math.Sin(α2) < 0 {
		return LatLon{}, ErrNoConvergence // ambiguous intersection
	}

	cosα3 := -math.Cos(α1)*math.Cos(α2) + math.Sin(α1)*math.Sin(α2)*math.Cos(δ12)
	δ13 := math.Atan2(math.Sin(δ12)*math.Sin(α1)*math.Sin(α2),
		math.Cos(α2)+math.Cos(α1)*cosα3)
	φ3 := math.Asin(clamp(math.Sin(φ1)*math.Cos(δ13)+
		math.Cos(φ1)*math.Sin(δ13)*math.Cos(θ13), -1, 1))
	dλ13 := math.Atan2(math.Sin(θ13)*math.Sin(δ13)*math.Cos(φ1),
		math.Cos(δ13)-math.Sin(φ1)*math.Sin(φ3))
	λ3 := λ1 + dλ13
	return LatLon{Lat: deg(φ3), Lon: NormalizeLongitude(deg(λ3))}, nil
}

// MaxLatitude returns the maximum latitude reached by a great circle that
// crosses the given point at the given bearing (degrees), via Clairaut's
// relation. The minimum latitude on the same circle is the negation.
func MaxLatitude(p LatLon, bearingDeg float64) float64 {
	θ := rad(bearingDeg)
	φ := rad(p.Lat)
	φmax := math.Acos(math.Abs(math.Sin(θ) * math.Cos(φ)))
	return deg(φmax)
}

// FinalBearingRad returns the bearing (radians, [0,2π)) at which the great
// circle from p1 arrives at p2.
func FinalBearingRad(p1, p2 LatLon) float64 {
	return NormalizeRadians(InitialBearingRad(p2, p1) + math.Pi)
}

// SphericalPolygonArea returns the area in square metres enclosed by the
// spherical polygon whose vertices are given in order (the polygon is closed
// automatically), on a sphere of the given radius. The edges are treated as
// great-circle arcs and the area is computed exactly from the spherical excess
// via the Gauss-Bonnet theorem, so it is correct for arbitrarily large
// polygons (as long as they are simple and do not enclose a pole). Winding
// order does not affect the magnitude.
func SphericalPolygonArea(vertices []LatLon, radius float64) float64 {
	n := len(vertices)
	if n < 3 {
		return 0
	}
	sum := 0.0
	for i := 0; i < n; i++ {
		prev := (i - 1 + n) % n
		next := (i + 1) % n
		arrive := FinalBearingRad(vertices[prev], vertices[i])
		depart := InitialBearingRad(vertices[i], vertices[next])
		sum += WrapPi(depart - arrive)
	}
	excess := math.Abs(2*math.Pi - sum)
	if excess > 2*math.Pi {
		excess = 4*math.Pi - excess
	}
	return excess * radius * radius
}

// SphericalPolygonPerimeter returns the great-circle perimeter in metres of the
// polygon (closed automatically) on a sphere of the given radius.
func SphericalPolygonPerimeter(vertices []LatLon, radius float64) float64 {
	n := len(vertices)
	if n < 2 {
		return 0
	}
	total := 0.0
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		total += HaversineDistanceR(vertices[i], vertices[j], radius)
	}
	return total
}

// SphericalExcess returns the spherical excess (E = A + B + C - π) in radians of
// a spherical triangle with the three given interior angles in radians. The
// triangle area on a sphere of radius R is E·R².
func SphericalExcess(a, b, c float64) float64 {
	return a + b + c - math.Pi
}

// SphericalTriangleArea returns the area in square metres of the spherical
// triangle with the three given vertices, on a sphere of the given radius.
func SphericalTriangleArea(p1, p2, p3 LatLon, radius float64) float64 {
	return SphericalPolygonArea([]LatLon{p1, p2, p3}, radius)
}

// PointInSphericalPolygon reports whether p lies inside the spherical polygon
// defined by vertices (in order, closed automatically), using a winding /
// ray-crossing test in longitude at the point's latitude. It assumes the
// polygon does not enclose a pole and edges are treated as rhumb-like segments
// in the lat/lon plane, which is adequate for modest polygons.
func PointInSphericalPolygon(p LatLon, vertices []LatLon) bool {
	n := len(vertices)
	if n < 3 {
		return false
	}
	inside := false
	j := n - 1
	for i := 0; i < n; i++ {
		yi, xi := vertices[i].Lat, vertices[i].Lon
		yj, xj := vertices[j].Lat, vertices[j].Lon
		if (yi > p.Lat) != (yj > p.Lat) {
			xcross := (xj-xi)*(p.Lat-yi)/(yj-yi) + xi
			if p.Lon < xcross {
				inside = !inside
			}
		}
		j = i
	}
	return inside
}

// clamp constrains x to the closed interval [lo, hi].
func clamp(x, lo, hi float64) float64 {
	if x < lo {
		return lo
	}
	if x > hi {
		return hi
	}
	return x
}
