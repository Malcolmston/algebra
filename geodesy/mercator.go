package geodesy

import "math"

// WebMercatorMaxLatitude is the latitude limit (~85.0511°) of the square Web
// Mercator (EPSG:3857) projection.
const WebMercatorMaxLatitude = 85.05112877980659

// MercatorForward projects a spherical latitude/longitude (degrees) to planar
// (x, y) coordinates in metres on a sphere of the given radius, with the
// projection centred on the prime meridian. x is eastings, y is northings.
func MercatorForward(lat, lon, radius float64) (x, y float64) {
	φ := rad(lat)
	x = radius * rad(lon)
	y = radius * math.Log(math.Tan(math.Pi/4+φ/2))
	return x, y
}

// MercatorInverse recovers spherical latitude/longitude (degrees) from planar
// Mercator (x, y) coordinates in metres on a sphere of the given radius.
func MercatorInverse(x, y, radius float64) (lat, lon float64) {
	lon = deg(x / radius)
	lat = deg(2*math.Atan(math.Exp(y/radius)) - math.Pi/2)
	return lat, lon
}

// WebMercatorForward projects WGS-84 latitude/longitude (degrees) to Web
// Mercator (EPSG:3857 / "Pseudo-Mercator") coordinates in metres, as used by
// most web map tile services. Latitudes are clamped to ±WebMercatorMaxLatitude.
func WebMercatorForward(lat, lon float64) (x, y float64) {
	if lat > WebMercatorMaxLatitude {
		lat = WebMercatorMaxLatitude
	} else if lat < -WebMercatorMaxLatitude {
		lat = -WebMercatorMaxLatitude
	}
	return MercatorForward(lat, lon, EarthRadiusEquatorial)
}

// WebMercatorInverse recovers WGS-84 latitude/longitude (degrees) from Web
// Mercator (EPSG:3857) coordinates in metres.
func WebMercatorInverse(x, y float64) (lat, lon float64) {
	return MercatorInverse(x, y, EarthRadiusEquatorial)
}

// EllipsoidalMercatorForward projects geodetic latitude/longitude (degrees) to
// the conformal ellipsoidal Mercator projection (metres) on the given
// ellipsoid, centred on the prime meridian.
func EllipsoidalMercatorForward(lat, lon float64, e Ellipsoid) (x, y float64) {
	φ := rad(lat)
	es := e.FirstEccentricity()
	esinφ := es * math.Sin(φ)
	x = e.A * rad(lon)
	y = e.A * math.Log(math.Tan(math.Pi/4+φ/2)*
		math.Pow((1-esinφ)/(1+esinφ), es/2))
	return x, y
}

// EllipsoidalMercatorInverse recovers geodetic latitude/longitude (degrees)
// from ellipsoidal Mercator coordinates (metres) on the given ellipsoid.
func EllipsoidalMercatorInverse(x, y float64, e Ellipsoid) (lat, lon float64) {
	es := e.FirstEccentricity()
	lon = deg(x / e.A)
	t := math.Exp(-y / e.A)
	φ := math.Pi/2 - 2*math.Atan(t)
	for i := 0; i < 20; i++ {
		esinφ := es * math.Sin(φ)
		φʹ := math.Pi/2 - 2*math.Atan(t*math.Pow((1-esinφ)/(1+esinφ), es/2))
		if math.Abs(φʹ-φ) < 1e-14 {
			φ = φʹ
			break
		}
		φ = φʹ
	}
	return deg(φ), lon
}

// LonLatToTile returns the slippy-map tile column and row (x, y) containing the
// given WGS-84 latitude/longitude (degrees) at the given zoom level, following
// the OpenStreetMap / Google/XYZ tiling scheme.
func LonLatToTile(lat, lon float64, zoom int) (tx, ty int) {
	n := math.Exp2(float64(zoom))
	tx = int(math.Floor((lon + 180) / 360 * n))
	φ := rad(lat)
	ty = int(math.Floor((1 - math.Log(math.Tan(φ)+1/math.Cos(φ))/math.Pi) / 2 * n))
	if tx < 0 {
		tx = 0
	}
	if maxT := int(n) - 1; tx > maxT {
		tx = maxT
	}
	if ty < 0 {
		ty = 0
	}
	if maxT := int(n) - 1; ty > maxT {
		ty = maxT
	}
	return tx, ty
}

// TileToLonLat returns the WGS-84 latitude/longitude (degrees) of the
// north-west corner of the slippy-map tile (tx, ty) at the given zoom level.
func TileToLonLat(tx, ty, zoom int) (lat, lon float64) {
	n := math.Exp2(float64(zoom))
	lon = float64(tx)/n*360 - 180
	φ := math.Atan(math.Sinh(math.Pi * (1 - 2*float64(ty)/n)))
	return deg(φ), lon
}
