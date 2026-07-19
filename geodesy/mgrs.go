package geodesy

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// Column letters for the MGRS 100 km squares, keyed by (zone-1) mod 3.
var mgrsE100k = [3]string{"ABCDEFGH", "JKLMNPQR", "STUVWXYZ"}

// Row letters for the MGRS 100 km squares, keyed by (zone-1) mod 2 (odd/even
// zone parity). The second string is the first rotated left by five letters.
var mgrsN100k = [2]string{"ABCDEFGHJKLMNPQRSTUV", "FGHJKLMNPQRSTUVABCDE"}

// mgrsBands is the ordered list of latitude band letters from C (-80°) to X.
const mgrsBands = "CDEFGHJKLMNPQRSTUVWX"

// UTMToMGRS encodes a UTM coordinate as an MGRS reference string with the given
// precision (number of digits per easting/northing, 1..5; 5 = 1 m resolution).
// The latitude band is derived from the supplied latitude (degrees), which must
// match the coordinate. Passing precision <= 0 yields the grid-square-only
// reference (for example "31UDQ").
func UTMToMGRS(u UTMCoord, lat float64, precision int) (string, error) {
	if u.Zone < 1 || u.Zone > 60 {
		return "", ErrInvalidUTM
	}
	band := MGRSLatBand(lat)
	if band == 0 {
		return "", ErrPolarRegion
	}
	col := int(math.Floor(u.Easting / 100000))
	if col < 1 || col > 8 {
		return "", ErrInvalidUTM
	}
	e100k := mgrsE100k[(u.Zone-1)%3][col-1]

	row := int(math.Floor(u.Northing/100000)) % 20
	if row < 0 {
		row += 20
	}
	n100k := mgrsN100k[(u.Zone-1)%2][row]

	var sb strings.Builder
	fmt.Fprintf(&sb, "%d%c%c%c", u.Zone, band, e100k, n100k)
	if precision > 0 {
		if precision > 5 {
			precision = 5
		}
		div := math.Pow(10, float64(5-precision))
		eLoc := int(math.Floor(math.Mod(u.Easting, 100000) / div))
		nLoc := int(math.Floor(math.Mod(u.Northing, 100000) / div))
		fmt.Fprintf(&sb, "%0*d%0*d", precision, eLoc, precision, nLoc)
	}
	return sb.String(), nil
}

// LatLonToMGRS encodes a WGS-84 geodetic position as an MGRS reference string
// with the given precision (1..5). It returns ErrPolarRegion outside the MGRS
// latitude band.
func LatLonToMGRS(p LatLon, precision int) (string, error) {
	u, err := LatLonToUTM(p)
	if err != nil {
		return "", err
	}
	return UTMToMGRS(u, p.Lat, precision)
}

// bandBottomNorthing returns the UTM northing (metres) of the southern edge of
// the given latitude band at the zone's central meridian.
func bandBottomNorthing(zone int, band byte) float64 {
	idx := strings.IndexByte(mgrsBands, band)
	if idx < 0 {
		return 0
	}
	bottomLat := float64(-80 + idx*8)
	lon0 := UTMCentralMeridian(zone)
	falseNorthing := 0.0
	if bottomLat < 0 {
		falseNorthing = UTMFalseNorthingSouth
	}
	_, northing := TransverseMercatorForward(
		bottomLat, lon0, lon0, UTMScaleFactor, UTMFalseEasting, falseNorthing, WGS84)
	return northing
}

// ParseMGRS decodes an MGRS reference string into a UTM coordinate. Whitespace
// is ignored, and the reference may omit the numeric location (grid square
// only), in which case the coordinate refers to the south-west corner of the
// 100 km square.
func ParseMGRS(s string) (UTMCoord, error) {
	s = strings.ToUpper(strings.Join(strings.Fields(s), ""))
	if len(s) < 5 {
		return UTMCoord{}, ErrInvalidMGRS
	}
	// Parse zone digits (1 or 2).
	i := 0
	for i < len(s) && s[i] >= '0' && s[i] <= '9' {
		i++
	}
	if i == 0 || i > 2 {
		return UTMCoord{}, ErrInvalidMGRS
	}
	zone, err := strconv.Atoi(s[:i])
	if err != nil || zone < 1 || zone > 60 {
		return UTMCoord{}, ErrInvalidMGRS
	}
	if i+3 > len(s) {
		return UTMCoord{}, ErrInvalidMGRS
	}
	band := s[i]
	if strings.IndexByte(mgrsBands, band) < 0 {
		return UTMCoord{}, ErrInvalidMGRS
	}
	e100k := s[i+1]
	n100k := s[i+2]
	digits := s[i+3:]
	if len(digits)%2 != 0 {
		return UTMCoord{}, ErrInvalidMGRS
	}

	col := strings.IndexByte(mgrsE100k[(zone-1)%3], e100k)
	if col < 0 {
		return UTMCoord{}, ErrInvalidMGRS
	}
	e100kNum := float64((col + 1) * 100000)

	row := strings.IndexByte(mgrsN100k[(zone-1)%2], n100k)
	if row < 0 {
		return UTMCoord{}, ErrInvalidMGRS
	}
	n100kNum := float64(row * 100000)

	var eLoc, nLoc float64
	if len(digits) > 0 {
		half := len(digits) / 2
		ev, err1 := strconv.Atoi(digits[:half])
		nv, err2 := strconv.Atoi(digits[half:])
		if err1 != nil || err2 != nil {
			return UTMCoord{}, ErrInvalidMGRS
		}
		mul := math.Pow(10, float64(5-half))
		eLoc = float64(ev) * mul
		nLoc = float64(nv) * mul
	}

	easting := e100kNum + eLoc

	nBand := math.Floor(bandBottomNorthing(zone, band)/100000) * 100000
	n2M := 0.0
	for n2M+n100kNum+nLoc < nBand {
		n2M += 2000000
	}
	northing := n2M + n100kNum + nLoc

	hemi := byte('N')
	if band < 'N' {
		hemi = 'S'
	}
	return UTMCoord{Zone: zone, Hemisphere: hemi, Easting: easting, Northing: northing}, nil
}

// MGRSToLatLon decodes an MGRS reference string into a WGS-84 geodetic position
// (the south-west corner of the referenced square when the numeric location is
// omitted, or of the referenced cell otherwise).
func MGRSToLatLon(s string) (LatLon, error) {
	u, err := ParseMGRS(s)
	if err != nil {
		return LatLon{}, err
	}
	return UTMToLatLon(u)
}
