package geodesy

import (
	"errors"
	"fmt"
	"math"
	"strconv"
	"strings"
)

// Common sentinel errors returned by the package.
var (
	// ErrNoConvergence is returned when an iterative solver (for example the
	// Vincenty inverse formula on near-antipodal points) fails to converge.
	ErrNoConvergence = errors.New("geodesy: iteration failed to converge")
	// ErrInvalidCoordinate indicates a latitude or longitude out of range.
	ErrInvalidCoordinate = errors.New("geodesy: invalid coordinate")
	// ErrInvalidMGRS indicates a malformed MGRS reference string.
	ErrInvalidMGRS = errors.New("geodesy: invalid MGRS reference")
	// ErrInvalidUTM indicates an out-of-range UTM zone or coordinate.
	ErrInvalidUTM = errors.New("geodesy: invalid UTM coordinate")
	// ErrPolarRegion indicates a latitude outside the UTM/MGRS band
	// (-80..84 degrees) where those grids are undefined.
	ErrPolarRegion = errors.New("geodesy: latitude outside UTM band")
)

// Ellipsoid is a biaxial reference ellipsoid defined by its semi-major axis A
// (metres) and flattening F.
type Ellipsoid struct {
	// A is the semi-major (equatorial) axis in metres.
	A float64
	// F is the flattening (a-b)/a.
	F float64
}

// Reference ellipsoids in common use.
var (
	// WGS84 is the World Geodetic System 1984 ellipsoid (GPS default).
	WGS84 = Ellipsoid{A: 6378137.0, F: 1 / 298.257223563}
	// GRS80 is the Geodetic Reference System 1980 ellipsoid.
	GRS80 = Ellipsoid{A: 6378137.0, F: 1 / 298.257222101}
	// Airy1830 is the Airy 1830 ellipsoid (Ordnance Survey Great Britain).
	Airy1830 = Ellipsoid{A: 6377563.396, F: 1 / 299.3249646}
	// International1924 is the Hayford / International 1924 ellipsoid.
	International1924 = Ellipsoid{A: 6378388.0, F: 1 / 297.0}
	// Clarke1866 is the Clarke 1866 ellipsoid (historic North America).
	Clarke1866 = Ellipsoid{A: 6378206.4, F: 1 / 294.9786982}
	// Bessel1841 is the Bessel 1841 ellipsoid.
	Bessel1841 = Ellipsoid{A: 6377397.155, F: 1 / 299.1528128}
	// SphereWGS84 is a sphere of the WGS-84 mean radius (F = 0).
	SphereWGS84 = Ellipsoid{A: EarthRadiusMean, F: 0}
)

// NewEllipsoid builds an ellipsoid from a semi-major axis and inverse
// flattening 1/f. Passing invFlattening <= 0 yields a sphere.
func NewEllipsoid(a, invFlattening float64) Ellipsoid {
	if invFlattening <= 0 {
		return Ellipsoid{A: a, F: 0}
	}
	return Ellipsoid{A: a, F: 1 / invFlattening}
}

// SemiMajorAxis returns the equatorial radius a in metres.
func (e Ellipsoid) SemiMajorAxis() float64 { return e.A }

// SemiMinorAxis returns the polar radius b = a(1-f) in metres.
func (e Ellipsoid) SemiMinorAxis() float64 { return e.A * (1 - e.F) }

// Flattening returns the flattening f = (a-b)/a.
func (e Ellipsoid) Flattening() float64 { return e.F }

// InverseFlattening returns 1/f, or +Inf for a sphere.
func (e Ellipsoid) InverseFlattening() float64 {
	if e.F == 0 {
		return math.Inf(1)
	}
	return 1 / e.F
}

// ThirdFlattening returns n = f/(2-f), used by the Krüger series.
func (e Ellipsoid) ThirdFlattening() float64 { return e.F / (2 - e.F) }

// FirstEccentricitySquared returns e² = f(2-f).
func (e Ellipsoid) FirstEccentricitySquared() float64 { return e.F * (2 - e.F) }

// FirstEccentricity returns the first eccentricity e.
func (e Ellipsoid) FirstEccentricity() float64 {
	return math.Sqrt(e.FirstEccentricitySquared())
}

// SecondEccentricitySquared returns e'² = e²/(1-e²).
func (e Ellipsoid) SecondEccentricitySquared() float64 {
	es := e.FirstEccentricitySquared()
	return es / (1 - es)
}

// SecondEccentricity returns the second eccentricity e'.
func (e Ellipsoid) SecondEccentricity() float64 {
	return math.Sqrt(e.SecondEccentricitySquared())
}

// IsSphere reports whether the ellipsoid is a perfect sphere (f == 0).
func (e Ellipsoid) IsSphere() bool { return e.F == 0 }

// MeanRadius returns the arithmetic mean radius R1 = (2a+b)/3.
func (e Ellipsoid) MeanRadius() float64 { return (2*e.A + e.SemiMinorAxis()) / 3 }

// AuthalicRadius returns the radius of a sphere with the same surface area.
func (e Ellipsoid) AuthalicRadius() float64 {
	if e.F == 0 {
		return e.A
	}
	es := e.FirstEccentricity()
	b := e.SemiMinorAxis()
	// R_A = sqrt( (a² + b²·atanh(e)/e) / 2 )
	term := b * b * math.Atanh(es) / es
	return math.Sqrt((e.A*e.A + term) / 2)
}

// VolumetricRadius returns the radius of a sphere with the same volume,
// (a²b)^(1/3).
func (e Ellipsoid) VolumetricRadius() float64 {
	return math.Cbrt(e.A * e.A * e.SemiMinorAxis())
}

// MeridianRadius returns the meridional radius of curvature M at the given
// geodetic latitude (degrees).
func (e Ellipsoid) MeridianRadius(latDeg float64) float64 {
	es := e.FirstEccentricitySquared()
	s := math.Sin(rad(latDeg))
	w := 1 - es*s*s
	return e.A * (1 - es) / (w * math.Sqrt(w))
}

// PrimeVerticalRadius returns the radius of curvature N in the prime vertical
// at the given geodetic latitude (degrees).
func (e Ellipsoid) PrimeVerticalRadius(latDeg float64) float64 {
	es := e.FirstEccentricitySquared()
	s := math.Sin(rad(latDeg))
	return e.A / math.Sqrt(1-es*s*s)
}

// GaussianRadius returns the Gaussian (geometric mean) radius of curvature
// sqrt(M·N) at the given geodetic latitude (degrees).
func (e Ellipsoid) GaussianRadius(latDeg float64) float64 {
	return math.Sqrt(e.MeridianRadius(latDeg) * e.PrimeVerticalRadius(latDeg))
}

// LatLon is a geodetic position given by latitude and longitude in degrees.
type LatLon struct {
	// Lat is the geodetic latitude in degrees, positive north.
	Lat float64
	// Lon is the longitude in degrees, positive east.
	Lon float64
}

// NewLatLon constructs a LatLon from latitude and longitude in degrees.
func NewLatLon(lat, lon float64) LatLon { return LatLon{Lat: lat, Lon: lon} }

// Valid reports whether the latitude is in [-90,90] and longitude in
// [-180,180].
func (p LatLon) Valid() bool {
	return p.Lat >= -90 && p.Lat <= 90 && p.Lon >= -180 && p.Lon <= 180
}

// Normalized returns the position with latitude clamped to [-90,90] and
// longitude wrapped to (-180,180].
func (p LatLon) Normalized() LatLon {
	return LatLon{Lat: NormalizeLatitude(p.Lat), Lon: NormalizeLongitude(p.Lon)}
}

// LatRad returns the latitude in radians.
func (p LatLon) LatRad() float64 { return rad(p.Lat) }

// LonRad returns the longitude in radians.
func (p LatLon) LonRad() float64 { return rad(p.Lon) }

// Antipode returns the point diametrically opposite on the sphere.
func (p LatLon) Antipode() LatLon {
	return LatLon{Lat: -p.Lat, Lon: NormalizeLongitude(p.Lon + 180)}
}

// Equal reports whether two positions are within tol degrees on both axes.
func (p LatLon) Equal(q LatLon, tol float64) bool {
	return math.Abs(p.Lat-q.Lat) <= tol && math.Abs(p.Lon-q.Lon) <= tol
}

// String renders the position in signed decimal degrees.
func (p LatLon) String() string {
	return fmt.Sprintf("%.6f, %.6f", p.Lat, p.Lon)
}

// sprintfDMS renders a degrees/minutes/seconds triple with a hemisphere
// suffix. It is used by FormatDMS.
func sprintfDMS(d, m int, s float64, hemi byte, secPrec int) string {
	return fmt.Sprintf("%d°%02d'%.*f\"%c", d, m, secPrec, s, hemi)
}

// ParseDMS parses a coordinate written in degrees/minutes/seconds into decimal
// degrees. It accepts flexible input such as "51°28'40.12\"N", "51 28 40.12 N",
// "-51.5", or "2°20′56″E". Hemisphere letters N/E give a positive result and
// S/W a negative one; a leading sign is also honoured.
func ParseDMS(s string) (float64, error) {
	orig := s
	s = strings.TrimSpace(s)
	if s == "" {
		return 0, fmt.Errorf("geodesy: empty DMS string")
	}
	// Determine hemisphere sign from a trailing/leading N/S/E/W.
	sign := 1.0
	hemiFound := false
	upper := strings.ToUpper(s)
	for _, h := range []byte{'N', 'S', 'E', 'W'} {
		if strings.IndexByte(upper, h) >= 0 {
			if h == 'S' || h == 'W' {
				sign = -1
			}
			hemiFound = true
			break
		}
	}
	// Replace all non-numeric separators with spaces.
	var b strings.Builder
	for _, r := range s {
		switch {
		case r >= '0' && r <= '9':
			b.WriteRune(r)
		case r == '.' || r == '-' || r == '+':
			b.WriteRune(r)
		default:
			b.WriteByte(' ')
		}
	}
	fields := strings.Fields(b.String())
	if len(fields) == 0 {
		return 0, fmt.Errorf("geodesy: no numbers in DMS string %q", orig)
	}
	vals := make([]float64, 0, len(fields))
	for _, f := range fields {
		v, err := strconv.ParseFloat(f, 64)
		if err != nil {
			return 0, fmt.Errorf("geodesy: bad number %q in %q", f, orig)
		}
		vals = append(vals, v)
	}
	// A leading explicit sign combined with a hemisphere letter would be
	// contradictory; the hemisphere letter wins only when no sign was given.
	explicitNeg := vals[0] < 0
	mag := math.Abs(vals[0])
	if len(vals) > 1 {
		mag += vals[1] / 60
	}
	if len(vals) > 2 {
		mag += vals[2] / 3600
	}
	if !hemiFound && explicitNeg {
		sign = -1
	}
	return sign * mag, nil
}
