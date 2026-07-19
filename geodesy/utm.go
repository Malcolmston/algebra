package geodesy

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// UTM grid constants.
const (
	// UTMScaleFactor is the central-meridian scale factor k0 for UTM (0.9996).
	UTMScaleFactor = 0.9996
	// UTMFalseEasting is the false easting applied to every UTM zone (metres).
	UTMFalseEasting = 500000.0
	// UTMFalseNorthingSouth is the false northing for the southern hemisphere.
	UTMFalseNorthingSouth = 10000000.0
)

// UTMCoord is a Universal Transverse Mercator grid coordinate.
type UTMCoord struct {
	// Zone is the UTM zone number (1..60).
	Zone int
	// Hemisphere is 'N' for northern or 'S' for southern.
	Hemisphere byte
	// Easting is the easting in metres.
	Easting float64
	// Northing is the northing in metres.
	Northing float64
}

// UTMZone returns the UTM zone number (1..60) for the given longitude and
// latitude (degrees), applying the Norway and Svalbard zone exceptions.
func UTMZone(lat, lon float64) int {
	lon = NormalizeLongitude(lon)
	zone := int(math.Floor((lon+180)/6)) + 1
	if zone < 1 {
		zone = 1
	}
	if zone > 60 {
		zone = 60
	}
	// Norway exception: zone 32 is widened for southern Norway.
	if lat >= 56 && lat < 64 && lon >= 3 && lon < 12 {
		zone = 32
	}
	// Svalbard exceptions.
	if lat >= 72 && lat < 84 {
		switch {
		case lon >= 0 && lon < 9:
			zone = 31
		case lon >= 9 && lon < 21:
			zone = 33
		case lon >= 21 && lon < 33:
			zone = 35
		case lon >= 33 && lon < 42:
			zone = 37
		}
	}
	return zone
}

// UTMCentralMeridian returns the central meridian (degrees) of the given UTM
// zone.
func UTMCentralMeridian(zone int) float64 {
	return float64((zone-1)*6 - 180 + 3)
}

// MGRSLatBand returns the MGRS/UTM latitude band letter for the given latitude
// (degrees). Bands run C (-80°) through X (84°), omitting I and O; latitudes
// outside [-80, 84) return a NUL byte.
func MGRSLatBand(lat float64) byte {
	if lat < -80 || lat >= 84 {
		return 0
	}
	bands := "CDEFGHJKLMNPQRSTUVWX"
	idx := int(math.Floor((lat + 80) / 8))
	if idx > 19 {
		idx = 19 // band X spans 12 degrees (72..84)
	}
	return bands[idx]
}

// LatLonToUTM converts a WGS-84 geodetic position to a UTM grid coordinate. It
// returns ErrPolarRegion for latitudes outside the UTM band [-80, 84).
func LatLonToUTM(p LatLon) (UTMCoord, error) {
	return LatLonToUTMEllipsoid(p, WGS84)
}

// LatLonToUTMEllipsoid converts a geodetic position to a UTM coordinate on the
// given ellipsoid.
func LatLonToUTMEllipsoid(p LatLon, e Ellipsoid) (UTMCoord, error) {
	if p.Lat < -80 || p.Lat >= 84 {
		return UTMCoord{}, ErrPolarRegion
	}
	zone := UTMZone(p.Lat, p.Lon)
	lon0 := UTMCentralMeridian(zone)
	hemi := byte('N')
	falseNorthing := 0.0
	if p.Lat < 0 {
		hemi = 'S'
		falseNorthing = UTMFalseNorthingSouth
	}
	easting, northing := TransverseMercatorForward(
		p.Lat, p.Lon, lon0, UTMScaleFactor, UTMFalseEasting, falseNorthing, e)
	return UTMCoord{Zone: zone, Hemisphere: hemi, Easting: easting, Northing: northing}, nil
}

// UTMToLatLon converts a UTM grid coordinate to a WGS-84 geodetic position.
func UTMToLatLon(u UTMCoord) (LatLon, error) {
	return UTMToLatLonEllipsoid(u, WGS84)
}

// UTMToLatLonEllipsoid converts a UTM coordinate to a geodetic position on the
// given ellipsoid.
func UTMToLatLonEllipsoid(u UTMCoord, e Ellipsoid) (LatLon, error) {
	if u.Zone < 1 || u.Zone > 60 {
		return LatLon{}, ErrInvalidUTM
	}
	lon0 := UTMCentralMeridian(u.Zone)
	falseNorthing := 0.0
	if u.Hemisphere == 'S' || u.Hemisphere == 's' {
		falseNorthing = UTMFalseNorthingSouth
	}
	lat, lon := TransverseMercatorInverse(
		u.Easting, u.Northing, lon0, UTMScaleFactor, UTMFalseEasting, falseNorthing, e)
	return LatLon{Lat: lat, Lon: NormalizeLongitude(lon)}, nil
}

// String renders a UTM coordinate as "zone hemisphere easting northing", for
// example "31 N 448252 5411932".
func (u UTMCoord) String() string {
	return fmt.Sprintf("%d %c %.0f %.0f", u.Zone, u.Hemisphere, u.Easting, u.Northing)
}

// ParseUTM parses a UTM string in the form produced by UTMCoord.String, for
// example "31 N 448252 5411932" or "31N 448252 5411932".
func ParseUTM(s string) (UTMCoord, error) {
	f := strings.Fields(strings.TrimSpace(s))
	// Allow the zone and hemisphere to be joined, e.g. "31N".
	if len(f) == 3 {
		// Split leading "31N" into "31" and "N".
		head := f[0]
		i := 0
		for i < len(head) && head[i] >= '0' && head[i] <= '9' {
			i++
		}
		if i > 0 && i < len(head) {
			f = append([]string{head[:i], head[i:]}, f[1:]...)
		}
	}
	if len(f) != 4 {
		return UTMCoord{}, ErrInvalidUTM
	}
	zone, err := strconv.Atoi(f[0])
	if err != nil {
		return UTMCoord{}, ErrInvalidUTM
	}
	hemi := strings.ToUpper(f[1])
	if hemi != "N" && hemi != "S" {
		return UTMCoord{}, ErrInvalidUTM
	}
	easting, err := strconv.ParseFloat(f[2], 64)
	if err != nil {
		return UTMCoord{}, ErrInvalidUTM
	}
	northing, err := strconv.ParseFloat(f[3], 64)
	if err != nil {
		return UTMCoord{}, ErrInvalidUTM
	}
	return UTMCoord{Zone: zone, Hemisphere: hemi[0], Easting: easting, Northing: northing}, nil
}
