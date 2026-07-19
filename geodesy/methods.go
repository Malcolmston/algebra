package geodesy

// DistanceTo returns the great-circle distance in metres from p to q on the
// WGS-84 mean sphere (Haversine).
func (p LatLon) DistanceTo(q LatLon) float64 { return HaversineDistance(p, q) }

// InitialBearingTo returns the initial great-circle bearing (degrees) from p to
// q.
func (p LatLon) InitialBearingTo(q LatLon) float64 { return InitialBearing(p, q) }

// FinalBearingTo returns the final great-circle bearing (degrees) from p to q.
func (p LatLon) FinalBearingTo(q LatLon) float64 { return FinalBearing(p, q) }

// MidpointTo returns the great-circle midpoint between p and q.
func (p LatLon) MidpointTo(q LatLon) LatLon { return Midpoint(p, q) }

// DestinationPoint returns the point reached from p along a great circle on the
// given bearing (degrees) for distance metres on the WGS-84 mean sphere.
func (p LatLon) DestinationPoint(bearing, distance float64) LatLon {
	return Destination(p, bearing, distance)
}

// RhumbDistanceTo returns the rhumb-line distance in metres from p to q on the
// WGS-84 mean sphere.
func (p LatLon) RhumbDistanceTo(q LatLon) float64 { return RhumbDistance(p, q) }

// RhumbBearingTo returns the constant rhumb-line bearing (degrees) from p to q.
func (p LatLon) RhumbBearingTo(q LatLon) float64 { return RhumbBearing(p, q) }

// GeodesicDistanceTo returns the geodesic distance in metres from p to q on the
// WGS-84 ellipsoid (Vincenty inverse).
func (p LatLon) GeodesicDistanceTo(q LatLon) (float64, error) {
	return VincentyDistance(p, q)
}

// CrossTrackDistanceTo returns the signed distance in metres of p from the
// great circle defined by start and end (WGS-84 mean sphere).
func (p LatLon) CrossTrackDistanceTo(start, end LatLon) float64 {
	return CrossTrackDistance(p, start, end)
}

// ToECEF converts p at the given ellipsoidal height (metres) to ECEF
// coordinates on the WGS-84 ellipsoid.
func (p LatLon) ToECEF(height float64) ECEF {
	return GeodeticToECEF(p.Lat, p.Lon, height, WGS84)
}

// ToUTM converts p to a UTM grid coordinate on WGS-84.
func (p LatLon) ToUTM() (UTMCoord, error) { return LatLonToUTM(p) }

// ToMGRS encodes p as an MGRS reference string with the given precision.
func (p LatLon) ToMGRS(precision int) (string, error) { return LatLonToMGRS(p, precision) }
