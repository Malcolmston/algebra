package geodesy

import "math"

// InverseResult holds the solution of the geodesic inverse problem: the
// geodesic distance between two points and the forward azimuths at each end.
type InverseResult struct {
	// Distance is the geodesic (shortest-path) distance in metres.
	Distance float64
	// InitialBearing is the forward azimuth at the first point (degrees).
	InitialBearing float64
	// FinalBearing is the forward azimuth at the second point (degrees).
	FinalBearing float64
}

// DirectResult holds the solution of the geodesic direct problem: the point
// reached from a start point along a geodesic, and the azimuth there.
type DirectResult struct {
	// Destination is the point reached.
	Destination LatLon
	// FinalBearing is the forward azimuth at the destination (degrees).
	FinalBearing float64
}

// VincentyInverse solves the inverse geodesic problem on the given ellipsoid
// using Vincenty's formula, returning the distance and initial/final bearings.
// It returns ErrNoConvergence for near-antipodal points where the algorithm
// does not converge.
func VincentyInverse(p1, p2 LatLon, e Ellipsoid) (InverseResult, error) {
	a := e.A
	f := e.F
	b := a * (1 - f)

	φ1, φ2 := rad(p1.Lat), rad(p2.Lat)
	L := rad(p2.Lon - p1.Lon)

	tanU1 := (1 - f) * math.Tan(φ1)
	cosU1 := 1 / math.Sqrt(1+tanU1*tanU1)
	sinU1 := tanU1 * cosU1
	tanU2 := (1 - f) * math.Tan(φ2)
	cosU2 := 1 / math.Sqrt(1+tanU2*tanU2)
	sinU2 := tanU2 * cosU2

	λ := L
	var sinλ, cosλ, sinσ, cosσ, σ, sinα, cos2α, cos2σm, C float64
	const maxIter = 200
	converged := false
	λʹ := 0.0
	for i := 0; i < maxIter; i++ {
		sinλ = math.Sin(λ)
		cosλ = math.Cos(λ)
		sinSqσ := (cosU2*sinλ)*(cosU2*sinλ) +
			(cosU1*sinU2-sinU1*cosU2*cosλ)*(cosU1*sinU2-sinU1*cosU2*cosλ)
		if sinSqσ == 0 {
			// Coincident points.
			return InverseResult{Distance: 0, InitialBearing: 0, FinalBearing: 0}, nil
		}
		sinσ = math.Sqrt(sinSqσ)
		cosσ = sinU1*sinU2 + cosU1*cosU2*cosλ
		σ = math.Atan2(sinσ, cosσ)
		sinα = cosU1 * cosU2 * sinλ / sinσ
		cos2α = 1 - sinα*sinα
		if cos2α != 0 {
			cos2σm = cosσ - 2*sinU1*sinU2/cos2α
		} else {
			cos2σm = 0 // equatorial line
		}
		C = f / 16 * cos2α * (4 + f*(4-3*cos2α))
		λʹ = λ
		λ = L + (1-C)*f*sinα*
			(σ+C*sinσ*(cos2σm+C*cosσ*(-1+2*cos2σm*cos2σm)))
		if math.Abs(λ-λʹ) < 1e-12 {
			converged = true
			break
		}
	}
	if !converged {
		return InverseResult{}, ErrNoConvergence
	}

	uSq := cos2α * (a*a - b*b) / (b * b)
	A := 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
	B := uSq / 1024 * (256 + uSq*(-128+uSq*(74-47*uSq)))
	Δσ := B * sinσ * (cos2σm + B/4*(cosσ*(-1+2*cos2σm*cos2σm)-
		B/6*cos2σm*(-3+4*sinσ*sinσ)*(-3+4*cos2σm*cos2σm)))
	s := b * A * (σ - Δσ)

	α1 := math.Atan2(cosU2*sinλ, cosU1*sinU2-sinU1*cosU2*cosλ)
	α2 := math.Atan2(cosU1*sinλ, -sinU1*cosU2+cosU1*sinU2*cosλ)

	return InverseResult{
		Distance:       s,
		InitialBearing: NormalizeDegrees(deg(α1)),
		FinalBearing:   NormalizeDegrees(deg(α2)),
	}, nil
}

// VincentyDirect solves the direct geodesic problem on the given ellipsoid:
// from a start point, travel distance metres along the initial azimuth
// (degrees) and return the destination and the azimuth there.
func VincentyDirect(start LatLon, bearingDeg, distance float64, e Ellipsoid) DirectResult {
	a := e.A
	f := e.F
	b := a * (1 - f)

	φ1 := rad(start.Lat)
	λ1 := rad(start.Lon)
	α1 := rad(bearingDeg)
	s := distance

	sinα1 := math.Sin(α1)
	cosα1 := math.Cos(α1)

	tanU1 := (1 - f) * math.Tan(φ1)
	cosU1 := 1 / math.Sqrt(1+tanU1*tanU1)
	sinU1 := tanU1 * cosU1

	σ1 := math.Atan2(tanU1, cosα1)
	sinα := cosU1 * sinα1
	cosSqα := 1 - sinα*sinα
	uSq := cosSqα * (a*a - b*b) / (b * b)
	A := 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
	B := uSq / 1024 * (256 + uSq*(-128+uSq*(74-47*uSq)))

	σ := s / (b * A)
	var sinσ, cosσ, cos2σm, Δσ float64
	const maxIter = 200
	σʹ := 0.0
	for i := 0; i < maxIter; i++ {
		cos2σm = math.Cos(2*σ1 + σ)
		sinσ = math.Sin(σ)
		cosσ = math.Cos(σ)
		Δσ = B * sinσ * (cos2σm + B/4*(cosσ*(-1+2*cos2σm*cos2σm)-
			B/6*cos2σm*(-3+4*sinσ*sinσ)*(-3+4*cos2σm*cos2σm)))
		σʹ = σ
		σ = s/(b*A) + Δσ
		if math.Abs(σ-σʹ) < 1e-12 {
			break
		}
	}

	x := sinU1*sinσ - cosU1*cosσ*cosα1
	φ2 := math.Atan2(sinU1*cosσ+cosU1*sinσ*cosα1,
		(1-f)*math.Hypot(sinα, x))
	λ := math.Atan2(sinσ*sinα1, cosU1*cosσ-sinU1*sinσ*cosα1)
	C := f / 16 * cosSqα * (4 + f*(4-3*cosSqα))
	L := λ - (1-C)*f*sinα*
		(σ+C*sinσ*(cos2σm+C*cosσ*(-1+2*cos2σm*cos2σm)))
	λ2 := λ1 + L
	α2 := math.Atan2(sinα, -x)

	return DirectResult{
		Destination:  LatLon{Lat: deg(φ2), Lon: NormalizeLongitude(deg(λ2))},
		FinalBearing: NormalizeDegrees(deg(α2)),
	}
}

// VincentyDistance returns the geodesic distance in metres between two points
// on the WGS-84 ellipsoid.
func VincentyDistance(p1, p2 LatLon) (float64, error) {
	r, err := VincentyInverse(p1, p2, WGS84)
	return r.Distance, err
}

// VincentyInitialBearing returns the initial azimuth (degrees) of the geodesic
// from p1 to p2 on the WGS-84 ellipsoid.
func VincentyInitialBearing(p1, p2 LatLon) (float64, error) {
	r, err := VincentyInverse(p1, p2, WGS84)
	return r.InitialBearing, err
}

// VincentyFinalBearing returns the final azimuth (degrees) of the geodesic from
// p1 to p2 on the WGS-84 ellipsoid.
func VincentyFinalBearing(p1, p2 LatLon) (float64, error) {
	r, err := VincentyInverse(p1, p2, WGS84)
	return r.FinalBearing, err
}

// VincentyDestination returns the point reached from start along a WGS-84
// geodesic on the given initial bearing (degrees) for distance metres.
func VincentyDestination(start LatLon, bearingDeg, distance float64) LatLon {
	return VincentyDirect(start, bearingDeg, distance, WGS84).Destination
}

// GeodesicInterpolate returns the point at the given fraction of the geodesic
// distance from p1 (fraction 0) to p2 (fraction 1) on the given ellipsoid.
func GeodesicInterpolate(p1, p2 LatLon, fraction float64, e Ellipsoid) (LatLon, error) {
	inv, err := VincentyInverse(p1, p2, e)
	if err != nil {
		return LatLon{}, err
	}
	if inv.Distance == 0 {
		return p1, nil
	}
	d := VincentyDirect(p1, inv.InitialBearing, inv.Distance*fraction, e)
	return d.Destination, nil
}

// GeodesicWaypoints returns n+1 points evenly spaced by geodesic distance from
// p1 to p2 inclusive on the given ellipsoid. n must be >= 1.
func GeodesicWaypoints(p1, p2 LatLon, n int, e Ellipsoid) ([]LatLon, error) {
	if n < 1 {
		n = 1
	}
	inv, err := VincentyInverse(p1, p2, e)
	if err != nil {
		return nil, err
	}
	pts := make([]LatLon, 0, n+1)
	for i := 0; i <= n; i++ {
		frac := float64(i) / float64(n)
		if inv.Distance == 0 {
			pts = append(pts, p1)
			continue
		}
		d := VincentyDirect(p1, inv.InitialBearing, inv.Distance*frac, e)
		pts = append(pts, d.Destination)
	}
	return pts, nil
}
