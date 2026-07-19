package geodesy

import "math"

// BoundingBox is an axis-aligned latitude/longitude rectangle, defined by its
// south-west (minimum) and north-east (maximum) corners. It does not handle
// boxes that span the ±180° antimeridian.
type BoundingBox struct {
	// Min is the south-west corner (minimum latitude and longitude).
	Min LatLon
	// Max is the north-east corner (maximum latitude and longitude).
	Max LatLon
}

// NewBoundingBox returns the bounding box with the given south-west and
// north-east corners, normalising the corner ordering.
func NewBoundingBox(sw, ne LatLon) BoundingBox {
	return BoundingBox{
		Min: LatLon{Lat: math.Min(sw.Lat, ne.Lat), Lon: math.Min(sw.Lon, ne.Lon)},
		Max: LatLon{Lat: math.Max(sw.Lat, ne.Lat), Lon: math.Max(sw.Lon, ne.Lon)},
	}
}

// BoundingBoxOfPoints returns the smallest bounding box containing all the given
// points. It returns the zero box for an empty slice.
func BoundingBoxOfPoints(points []LatLon) BoundingBox {
	if len(points) == 0 {
		return BoundingBox{}
	}
	minLat, maxLat := points[0].Lat, points[0].Lat
	minLon, maxLon := points[0].Lon, points[0].Lon
	for _, p := range points[1:] {
		minLat = math.Min(minLat, p.Lat)
		maxLat = math.Max(maxLat, p.Lat)
		minLon = math.Min(minLon, p.Lon)
		maxLon = math.Max(maxLon, p.Lon)
	}
	return BoundingBox{Min: LatLon{minLat, minLon}, Max: LatLon{maxLat, maxLon}}
}

// BoundingBoxAround returns a bounding box centred on the given point that
// encloses a circle of the given radius (metres), computed on the WGS-84 mean
// sphere. The box is a latitude/longitude rectangle tangent to that circle.
func BoundingBoxAround(center LatLon, radius float64) BoundingBox {
	dLat := deg(radius / EarthRadiusMean)
	cosLat := math.Cos(rad(center.Lat))
	dLon := dLat
	if math.Abs(cosLat) > 1e-12 {
		dLon = deg(radius / (EarthRadiusMean * cosLat))
	}
	return BoundingBox{
		Min: LatLon{Lat: center.Lat - dLat, Lon: center.Lon - dLon},
		Max: LatLon{Lat: center.Lat + dLat, Lon: center.Lon + dLon},
	}
}

// Contains reports whether the point lies within the box (inclusive of edges).
func (b BoundingBox) Contains(p LatLon) bool {
	return p.Lat >= b.Min.Lat && p.Lat <= b.Max.Lat &&
		p.Lon >= b.Min.Lon && p.Lon <= b.Max.Lon
}

// Center returns the centre point of the box.
func (b BoundingBox) Center() LatLon {
	return LatLon{Lat: (b.Min.Lat + b.Max.Lat) / 2, Lon: (b.Min.Lon + b.Max.Lon) / 2}
}

// Corners returns the four corners of the box in the order south-west,
// south-east, north-east, north-west.
func (b BoundingBox) Corners() [4]LatLon {
	return [4]LatLon{
		{Lat: b.Min.Lat, Lon: b.Min.Lon},
		{Lat: b.Min.Lat, Lon: b.Max.Lon},
		{Lat: b.Max.Lat, Lon: b.Max.Lon},
		{Lat: b.Max.Lat, Lon: b.Min.Lon},
	}
}

// SouthWest returns the south-west corner of the box.
func (b BoundingBox) SouthWest() LatLon { return b.Min }

// NorthEast returns the north-east corner of the box.
func (b BoundingBox) NorthEast() LatLon { return b.Max }

// HeightMeters returns the north-south extent of the box in metres on the
// WGS-84 mean sphere.
func (b BoundingBox) HeightMeters() float64 {
	return rad(b.Max.Lat-b.Min.Lat) * EarthRadiusMean
}

// WidthMeters returns the east-west extent of the box in metres at its centre
// latitude on the WGS-84 mean sphere.
func (b BoundingBox) WidthMeters() float64 {
	c := b.Center()
	return rad(b.Max.Lon-b.Min.Lon) * EarthRadiusMean * math.Cos(rad(c.Lat))
}

// Expand returns a copy of the box grown outward by the given amount of metres
// on all sides (WGS-84 mean sphere).
func (b BoundingBox) Expand(meters float64) BoundingBox {
	dLat := deg(meters / EarthRadiusMean)
	c := b.Center()
	cosLat := math.Cos(rad(c.Lat))
	dLon := dLat
	if math.Abs(cosLat) > 1e-12 {
		dLon = deg(meters / (EarthRadiusMean * cosLat))
	}
	return BoundingBox{
		Min: LatLon{Lat: b.Min.Lat - dLat, Lon: b.Min.Lon - dLon},
		Max: LatLon{Lat: b.Max.Lat + dLat, Lon: b.Max.Lon + dLon},
	}
}

// Union returns the smallest box containing both b and o.
func (b BoundingBox) Union(o BoundingBox) BoundingBox {
	return BoundingBox{
		Min: LatLon{Lat: math.Min(b.Min.Lat, o.Min.Lat), Lon: math.Min(b.Min.Lon, o.Min.Lon)},
		Max: LatLon{Lat: math.Max(b.Max.Lat, o.Max.Lat), Lon: math.Max(b.Max.Lon, o.Max.Lon)},
	}
}

// Intersects reports whether two boxes overlap (sharing an edge counts).
func (b BoundingBox) Intersects(o BoundingBox) bool {
	return b.Min.Lat <= o.Max.Lat && b.Max.Lat >= o.Min.Lat &&
		b.Min.Lon <= o.Max.Lon && b.Max.Lon >= o.Min.Lon
}
