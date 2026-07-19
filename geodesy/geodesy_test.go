package geodesy

import (
	"fmt"
	"math"
	"testing"
)

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func dms(d, m, s float64) float64 {
	sign := 1.0
	if d < 0 {
		sign, d = -1, -d
	}
	return sign * (d + m/60 + s/3600)
}

// --- angle helpers ---------------------------------------------------------

func TestNormalizeDegrees(t *testing.T) {
	tests := []struct{ in, want float64 }{
		{0, 0}, {360, 0}, {370, 10}, {-10, 350}, {720, 0}, {-350, 10}, {180, 180},
	}
	for _, tc := range tests {
		if got := NormalizeDegrees(tc.in); !approx(got, tc.want, 1e-9) {
			t.Errorf("NormalizeDegrees(%v)=%v want %v", tc.in, got, tc.want)
		}
	}
}

func TestNormalizeLongitude(t *testing.T) {
	tests := []struct{ in, want float64 }{
		{0, 0}, {180, -180}, {-180, -180}, {190, -170}, {-190, 170}, {540, -180},
	}
	for _, tc := range tests {
		if got := NormalizeLongitude(tc.in); !approx(got, tc.want, 1e-9) {
			t.Errorf("NormalizeLongitude(%v)=%v want %v", tc.in, got, tc.want)
		}
	}
}

func TestWrapPi(t *testing.T) {
	if got := WrapPi(3 * math.Pi); !approx(got, math.Pi, 1e-9) {
		t.Errorf("WrapPi(3pi)=%v want pi", got)
	}
	if got := WrapPi(-3 * math.Pi); !approx(got, math.Pi, 1e-9) {
		t.Errorf("WrapPi(-3pi)=%v want pi", got)
	}
}

func TestDegRadRoundtrip(t *testing.T) {
	for _, d := range []float64{0, 30, 45, 90, 123.456, -270} {
		if got := RadToDeg(DegToRad(d)); !approx(got, d, 1e-12) {
			t.Errorf("deg/rad roundtrip %v -> %v", d, got)
		}
	}
}

func TestDMS(t *testing.T) {
	v := DMSToDegrees(51, 28, 40.12)
	if !approx(v, 51.477811111, 1e-8) {
		t.Fatalf("DMSToDegrees=%v", v)
	}
	d, m, s := DegreesToDMS(v)
	if d != 51 || m != 28 || !approx(s, 40.12, 1e-6) {
		t.Errorf("DegreesToDMS=%d %d %v", d, m, s)
	}
	p, err := ParseDMS("51°28'40.12\"N")
	if err != nil || !approx(p, 51.477811111, 1e-8) {
		t.Errorf("ParseDMS=%v err=%v", p, err)
	}
	w, err := ParseDMS("2°20'56\"W")
	if err != nil || !approx(w, -2.348888888, 1e-8) {
		t.Errorf("ParseDMS W=%v err=%v", w, err)
	}
	if neg, _ := ParseDMS("-3.5"); !approx(neg, -3.5, 1e-9) {
		t.Errorf("ParseDMS signed=%v", neg)
	}
}

// --- ellipsoid -------------------------------------------------------------

func TestEllipsoid(t *testing.T) {
	if !approx(WGS84.SemiMinorAxis(), 6356752.314245, 1e-5) {
		t.Errorf("b=%v", WGS84.SemiMinorAxis())
	}
	if !approx(WGS84.FirstEccentricitySquared(), 0.00669437999014, 1e-12) {
		t.Errorf("e2=%v", WGS84.FirstEccentricitySquared())
	}
	if !approx(WGS84.InverseFlattening(), 298.257223563, 1e-9) {
		t.Errorf("1/f=%v", WGS84.InverseFlattening())
	}
	if !approx(WGS84.MeanRadius(), EarthRadiusMean, 1e-3) {
		t.Errorf("R1=%v", WGS84.MeanRadius())
	}
	if !SphereWGS84.IsSphere() {
		t.Error("SphereWGS84 should be a sphere")
	}
	// Prime-vertical radius at the equator equals a, at the pole equals a/sqrt(1-e2).
	if !approx(WGS84.PrimeVerticalRadius(0), WGS84.A, 1e-6) {
		t.Errorf("N(0)=%v", WGS84.PrimeVerticalRadius(0))
	}
}

// --- great circle ----------------------------------------------------------

func TestHaversine(t *testing.T) {
	// One degree of latitude on the mean sphere.
	d := HaversineDistance(LatLon{0, 0}, LatLon{1, 0})
	if !approx(d, EarthRadiusMean*rad(1), 1e-6) {
		t.Errorf("1deg lat=%v", d)
	}
	// JFK to LHR, known ~5540 km on the mean sphere.
	jfk, lhr := LatLon{40.6413, -73.7781}, LatLon{51.4700, -0.4543}
	if got := HaversineDistance(jfk, lhr); !approx(got, 5540019, 2000) {
		t.Errorf("JFK-LHR=%v", got)
	}
	// Law of cosines should agree with Haversine for a moderate distance.
	if lc := SphericalLawOfCosinesDistance(jfk, lhr, EarthRadiusMean); !approx(lc, HaversineDistance(jfk, lhr), 1) {
		t.Errorf("law-of-cosines disagrees: %v", lc)
	}
}

func TestBearings(t *testing.T) {
	if b := InitialBearing(LatLon{0, 0}, LatLon{0, 10}); !approx(b, 90, 1e-6) {
		t.Errorf("east bearing=%v", b)
	}
	if b := InitialBearing(LatLon{0, 0}, LatLon{10, 0}); !approx(b, 0, 1e-6) {
		t.Errorf("north bearing=%v", b)
	}
	if b := InitialBearing(LatLon{0, 0}, LatLon{0, -10}); !approx(b, 270, 1e-6) {
		t.Errorf("west bearing=%v", b)
	}
	// Final bearing along a great circle differs from initial (except on meridians/equator).
	fb := FinalBearing(LatLon{40.6413, -73.7781}, LatLon{51.47, -0.4543})
	if fb <= 90 || fb >= 180 {
		t.Errorf("final bearing out of expected range: %v", fb)
	}
}

func TestMidpointAndInterpolate(t *testing.T) {
	m := Midpoint(LatLon{0, 0}, LatLon{0, 90})
	if !approx(m.Lat, 0, 1e-9) || !approx(m.Lon, 45, 1e-9) {
		t.Errorf("midpoint=%v", m)
	}
	// Intermediate at 0.5 equals midpoint.
	mi := IntermediatePoint(LatLon{0, 0}, LatLon{0, 90}, 0.5)
	if !approx(mi.Lat, m.Lat, 1e-9) || !approx(mi.Lon, m.Lon, 1e-9) {
		t.Errorf("intermediate 0.5=%v", mi)
	}
	pts := InterpolatePoints(LatLon{0, 0}, LatLon{0, 90}, 3)
	if len(pts) != 4 || !approx(pts[3].Lon, 90, 1e-9) {
		t.Errorf("interpolate=%v", pts)
	}
}

func TestDestination(t *testing.T) {
	// East along the equator a quarter of the way round.
	d := DestinationR(LatLon{0, 0}, 90, EarthRadiusMean*math.Pi/2, EarthRadiusMean)
	if !approx(d.Lat, 0, 1e-6) || !approx(d.Lon, 90, 1e-6) {
		t.Errorf("dest=%v", d)
	}
	// Veness worked example.
	d2 := DestinationR(LatLon{51.47788, -0.00147}, 300.7, 7794, 6371000)
	if !approx(d2.Lat, 51.51363, 1e-4) || !approx(d2.Lon, -0.09832, 1e-4) {
		t.Errorf("dest2=%v", d2)
	}
	// Destination then reverse bearing returns to start.
	back := DestinationR(d2, BackAzimuth(300.7), 7794, 6371000)
	if !approx(back.Lat, 51.47788, 1e-3) {
		t.Errorf("dest roundtrip lat=%v", back.Lat)
	}
}

func TestCrossAndAlongTrack(t *testing.T) {
	p := LatLon{53.2611, -0.7972}
	a := LatLon{53.3206, -1.7297}
	b := LatLon{53.1887, 0.1334}
	if xtd := CrossTrackDistanceR(p, a, b, 6371000); !approx(xtd, -307.55, 1) {
		t.Errorf("xtd=%v", xtd)
	}
	if atd := AlongTrackDistanceR(p, a, b, 6371000); !approx(atd, 62331.49, 5) {
		t.Errorf("atd=%v", atd)
	}
}

func TestIntersection(t *testing.T) {
	i, err := Intersection(LatLon{51.8853, 0.2545}, 108.547, LatLon{49.0034, 2.5735}, 32.435)
	if err != nil {
		t.Fatalf("intersection err=%v", err)
	}
	if !approx(i.Lat, 50.9078, 1e-3) || !approx(i.Lon, 4.5084, 1e-3) {
		t.Errorf("intersection=%v", i)
	}
}

func TestSphericalArea(t *testing.T) {
	// Spherical octant (three right-angle vertices): area = pi/2 * R^2.
	oct := SphericalPolygonArea([]LatLon{{0, 0}, {0, 90}, {90, 0}}, 1)
	if !approx(oct, math.Pi/2, 1e-9) {
		t.Errorf("octant area=%v want pi/2", oct)
	}
	// A small triangle checked against L'Huilier's theorem.
	tri := []LatLon{{40, -100}, {45, -90}, {35, -95}}
	want := lhuilier(tri[0], tri[1], tri[2])
	if got := SphericalPolygonArea(tri, 1); !approx(got, want, 1e-9) {
		t.Errorf("triangle area=%v want %v", got, want)
	}
	// Perimeter of a square of side ~10deg.
	sq := []LatLon{{0, 0}, {0, 10}, {10, 10}, {10, 0}}
	per := SphericalPolygonPerimeter(sq, EarthRadiusMean)
	if per <= 0 {
		t.Errorf("perimeter=%v", per)
	}
}

func lhuilier(p1, p2, p3 LatLon) float64 {
	a := AngularDistance(p2, p3)
	b := AngularDistance(p1, p3)
	c := AngularDistance(p1, p2)
	s := (a + b + c) / 2
	tn := math.Tan(s/2) * math.Tan((s-a)/2) * math.Tan((s-b)/2) * math.Tan((s-c)/2)
	if tn < 0 {
		tn = 0
	}
	return 4 * math.Atan(math.Sqrt(tn))
}

func TestAntipode(t *testing.T) {
	p := LatLon{45, 90}
	a := p.Antipode()
	if !approx(a.Lat, -45, 1e-9) || !approx(a.Lon, -90, 1e-9) {
		t.Errorf("antipode=%v", a)
	}
}

// --- rhumb -----------------------------------------------------------------

func TestRhumb(t *testing.T) {
	if b := RhumbBearing(LatLon{0, 0}, LatLon{0, 90}); !approx(b, 90, 1e-9) {
		t.Errorf("rhumb bearing eq=%v", b)
	}
	d := RhumbDistanceR(LatLon{51.127, 1.338}, LatLon{50.964, 1.853}, 6371000)
	if !approx(d, 40307.7, 5) {
		t.Errorf("rhumb dover=%v", d)
	}
	b := RhumbBearing(LatLon{51.127, 1.338}, LatLon{50.964, 1.853})
	if !approx(b, 116.7219, 1e-3) {
		t.Errorf("rhumb bearing dover=%v", b)
	}
	// Rhumb destination inverts distance/bearing.
	dst := RhumbDestinationR(LatLon{51.127, 1.338}, b, d, 6371000)
	if !approx(dst.Lat, 50.964, 1e-3) || !approx(dst.Lon, 1.853, 1e-3) {
		t.Errorf("rhumb dest=%v", dst)
	}
	mid := RhumbMidpoint(LatLon{51.127, 1.338}, LatLon{50.964, 1.853})
	if mid.Lat < 50.9 || mid.Lat > 51.13 {
		t.Errorf("rhumb midpoint=%v", mid)
	}
}

// --- Vincenty --------------------------------------------------------------

func TestVincenty(t *testing.T) {
	fp := LatLon{dms(-37, 57, 3.72030), dms(144, 25, 29.52440)}
	bn := LatLon{dms(-37, 39, 10.15610), dms(143, 55, 35.38390)}
	inv, err := VincentyInverse(fp, bn, WGS84)
	if err != nil {
		t.Fatalf("inverse err=%v", err)
	}
	if !approx(inv.Distance, 54972.271, 1e-3) {
		t.Errorf("distance=%v want 54972.271", inv.Distance)
	}
	if !approx(inv.InitialBearing, 306.868159, 1e-5) {
		t.Errorf("initial bearing=%v", inv.InitialBearing)
	}
	if !approx(inv.FinalBearing, 307.173631, 1e-5) {
		t.Errorf("final bearing=%v", inv.FinalBearing)
	}
	// Direct problem should reproduce the second point.
	dir := VincentyDirect(fp, inv.InitialBearing, inv.Distance, WGS84)
	if !approx(dir.Destination.Lat, bn.Lat, 1e-9) || !approx(dir.Destination.Lon, bn.Lon, 1e-9) {
		t.Errorf("direct=%v want %v", dir.Destination, bn)
	}
	// Coincident points.
	if got, _ := VincentyDistance(fp, fp); !approx(got, 0, 1e-6) {
		t.Errorf("coincident distance=%v", got)
	}
	// Vincenty on WGS-84 should exceed the spherical Haversine distance slightly
	// but stay within a fraction of a percent for this line.
	hv := HaversineDistance(fp, bn)
	if math.Abs(hv-inv.Distance)/inv.Distance > 0.01 {
		t.Errorf("haversine %v vs vincenty %v too different", hv, inv.Distance)
	}
}

func TestGeodesicInterpolate(t *testing.T) {
	p1, p2 := LatLon{50, 0}, LatLon{52, 5}
	mid, err := GeodesicInterpolate(p1, p2, 0.5, WGS84)
	if err != nil {
		t.Fatalf("interp err=%v", err)
	}
	// Distances from each end to the midpoint should be equal.
	d1, _ := VincentyDistance(p1, mid)
	d2, _ := VincentyDistance(mid, p2)
	if !approx(d1, d2, 1e-3) {
		t.Errorf("geodesic midpoint not centered: %v vs %v", d1, d2)
	}
	wps, err := GeodesicWaypoints(p1, p2, 4, WGS84)
	if err != nil || len(wps) != 5 {
		t.Fatalf("waypoints=%v err=%v", wps, err)
	}
}

// --- ECEF / ENU ------------------------------------------------------------

func TestECEFRoundtrip(t *testing.T) {
	cases := []struct {
		lat, lon, h float64
	}{
		{48.8583, 2.2945, 100},
		{-33.8688, 151.2093, 58},
		{0, 0, 0},
		{89.9, 45, 1000},
	}
	for _, c := range cases {
		e := GeodeticToECEF(c.lat, c.lon, c.h, WGS84)
		lat, lon, h := ECEFToGeodetic(e, WGS84)
		if !approx(lat, c.lat, 1e-7) || !approx(lon, c.lon, 1e-7) || !approx(h, c.h, 1e-4) {
			t.Errorf("ECEF roundtrip %v -> %v %v %v", c, lat, lon, h)
		}
	}
	// Known ECEF magnitude near the surface ~ Earth radius.
	e := GeodeticToECEF(0, 0, 0, WGS84)
	if !approx(e.X, WGS84.A, 1e-3) || !approx(e.Y, 0, 1e-6) || !approx(e.Z, 0, 1e-6) {
		t.Errorf("ECEF(0,0,0)=%v", e)
	}
}

func TestENUAndAER(t *testing.T) {
	ref := LatLon{45, 9}
	// A pure eastward longitude step maps to ENU with ~zero north component and
	// an east component matching the parallel arc length.
	tgt := LatLon{ref.Lat, ref.Lon + 0.01}
	enu := GeodeticToENU(tgt.Lat, tgt.Lon, 0, ref.Lat, ref.Lon, 0, WGS84)
	if math.Abs(enu.N) > 0.5 {
		t.Errorf("east step north component=%v", enu.N)
	}
	if want := ParallelArcLength(ref.Lat, 0.01, WGS84); !approx(enu.E, want, 0.5) {
		t.Errorf("east ENU=%v want ~%v", enu.E, want)
	}
	// ENU -> ECEF -> geodetic roundtrip.
	lat, lon, _ := ENUToGeodetic(enu, ref.Lat, ref.Lon, 0, WGS84)
	if !approx(lat, tgt.Lat, 1e-7) || !approx(lon, tgt.Lon, 1e-7) {
		t.Errorf("ENU roundtrip=%v %v", lat, lon)
	}
	// AER roundtrip.
	aer := ENUToAER(ENU{E: 300, N: 400, U: 500})
	got := AERToENU(aer)
	if !approx(got.E, 300, 1e-6) || !approx(got.N, 400, 1e-6) || !approx(got.U, 500, 1e-6) {
		t.Errorf("AER roundtrip=%v", got)
	}
	if !approx(aer.Range, math.Sqrt(300*300+400*400+500*500), 1e-6) {
		t.Errorf("AER range=%v", aer.Range)
	}
}

// --- Mercator / tiles ------------------------------------------------------

func TestWebMercator(t *testing.T) {
	x, y := WebMercatorForward(0, 0)
	if !approx(x, 0, 1e-6) || !approx(y, 0, 1e-6) {
		t.Errorf("webmerc 0,0=%v %v", x, y)
	}
	// Full extent at max latitude ~ +/- pi*a.
	_, ymax := WebMercatorForward(WebMercatorMaxLatitude, 0)
	if !approx(ymax, math.Pi*EarthRadiusEquatorial, 1) {
		t.Errorf("webmerc ymax=%v", ymax)
	}
	// Roundtrip.
	lat, lon := WebMercatorInverse(WebMercatorForward(48.8583, 2.2945))
	if !approx(lat, 48.8583, 1e-9) || !approx(lon, 2.2945, 1e-9) {
		t.Errorf("webmerc roundtrip=%v %v", lat, lon)
	}
	// Ellipsoidal Mercator roundtrip.
	ex, ey := EllipsoidalMercatorForward(48.8583, 2.2945, WGS84)
	elat, elon := EllipsoidalMercatorInverse(ex, ey, WGS84)
	if !approx(elat, 48.8583, 1e-7) || !approx(elon, 2.2945, 1e-9) {
		t.Errorf("ell merc roundtrip=%v %v", elat, elon)
	}
}

func TestTiles(t *testing.T) {
	// Zoom 0 tile 0 north-west corner.
	lat, lon := TileToLonLat(0, 0, 0)
	if !approx(lon, -180, 1e-9) || !approx(lat, WebMercatorMaxLatitude, 1e-6) {
		t.Errorf("tile 0/0/0=%v %v", lat, lon)
	}
	tx, ty := LonLatToTile(0, 0, 1)
	if tx != 1 || ty != 1 {
		t.Errorf("tile of 0,0 at z1=%d %d", tx, ty)
	}
}

// --- Transverse Mercator / UTM / MGRS --------------------------------------

func TestTransverseMercatorMeridian(t *testing.T) {
	// On the central meridian, northing = k0 * meridian arc.
	lat := 48.8583
	_, n := TransverseMercatorForward(lat, 3, 3, 0.9996, 500000, 0, WGS84)
	want := 0.9996 * MeridianArcLength(lat, WGS84)
	if !approx(n, want, 1e-3) {
		t.Errorf("TM meridian northing=%v want %v", n, want)
	}
}

func TestUTM(t *testing.T) {
	eiffel := LatLon{48.8583, 2.2945}
	u, err := LatLonToUTM(eiffel)
	if err != nil {
		t.Fatalf("UTM err=%v", err)
	}
	if u.Zone != 31 || u.Hemisphere != 'N' {
		t.Errorf("UTM zone/hemi=%d %c", u.Zone, u.Hemisphere)
	}
	if !approx(u.Easting, 448252, 1) {
		t.Errorf("UTM easting=%v", u.Easting)
	}
	// Roundtrip.
	back, err := UTMToLatLon(u)
	if err != nil {
		t.Fatalf("UTM inverse err=%v", err)
	}
	if !approx(back.Lat, eiffel.Lat, 1e-7) || !approx(back.Lon, eiffel.Lon, 1e-7) {
		t.Errorf("UTM roundtrip=%v", back)
	}
	// Southern hemisphere.
	syd := LatLon{-33.8688, 151.2093}
	us, _ := LatLonToUTM(syd)
	if us.Hemisphere != 'S' || us.Zone != 56 {
		t.Errorf("Sydney UTM=%d %c", us.Zone, us.Hemisphere)
	}
	bs, _ := UTMToLatLon(us)
	if !approx(bs.Lat, syd.Lat, 1e-7) || !approx(bs.Lon, syd.Lon, 1e-7) {
		t.Errorf("Sydney roundtrip=%v", bs)
	}
	// Zone exceptions.
	if UTMZone(60, 5) != 32 {
		t.Errorf("Norway zone=%d", UTMZone(60, 5))
	}
	if UTMZone(75, 20) != 33 {
		t.Errorf("Svalbard zone=%d", UTMZone(75, 20))
	}
	// Parse/format roundtrip.
	pu, err := ParseUTM(u.String())
	if err != nil || pu.Zone != u.Zone {
		t.Errorf("ParseUTM=%v err=%v", pu, err)
	}
}

func TestUTMPolarError(t *testing.T) {
	if _, err := LatLonToUTM(LatLon{85, 0}); err != ErrPolarRegion {
		t.Errorf("expected polar error, got %v", err)
	}
}

func TestMGRS(t *testing.T) {
	eiffel := LatLon{48.8583, 2.2945}
	m, err := LatLonToMGRS(eiffel, 5)
	if err != nil {
		t.Fatalf("MGRS err=%v", err)
	}
	// Grid-zone designator and 100 km square for the Eiffel Tower.
	if m[:5] != "31UDQ" {
		t.Errorf("MGRS prefix=%q (full %q)", m[:5], m)
	}
	// Decode should land within a metre of the encoded cell corner.
	ll, err := MGRSToLatLon(m)
	if err != nil {
		t.Fatalf("MGRS decode err=%v", err)
	}
	if !approx(ll.Lat, eiffel.Lat, 1e-4) || !approx(ll.Lon, eiffel.Lon, 1e-4) {
		t.Errorf("MGRS roundtrip=%v from %q", ll, m)
	}
	// Latitude bands.
	if MGRSLatBand(0) != 'N' {
		t.Errorf("band at 0=%c", MGRSLatBand(0))
	}
	if MGRSLatBand(-80) != 'C' || MGRSLatBand(72) != 'X' {
		t.Errorf("band edges: %c %c", MGRSLatBand(-80), MGRSLatBand(72))
	}
	// Southern hemisphere roundtrip.
	syd := LatLon{-33.8688, 151.2093}
	ms, _ := LatLonToMGRS(syd, 5)
	lls, err := MGRSToLatLon(ms)
	if err != nil {
		t.Fatalf("Sydney MGRS decode err=%v", err)
	}
	if !approx(lls.Lat, syd.Lat, 1e-4) || !approx(lls.Lon, syd.Lon, 1e-4) {
		t.Errorf("Sydney MGRS roundtrip=%v from %q", lls, ms)
	}
	// Whitespace tolerance.
	if _, err := ParseMGRS("31U DQ 48251 11932"); err != nil {
		t.Errorf("ParseMGRS whitespace err=%v", err)
	}
	if _, err := ParseMGRS("not-mgrs"); err == nil {
		t.Error("expected error for bad MGRS")
	}
}

// --- conversions & compass -------------------------------------------------

func TestConversions(t *testing.T) {
	if !approx(NauticalMilesToMeters(1), 1852, 1e-9) {
		t.Error("nautical mile")
	}
	if !approx(MetersToFeet(FeetToMeters(100)), 100, 1e-9) {
		t.Error("feet roundtrip")
	}
	if !approx(DegToGrad(90), 100, 1e-9) {
		t.Error("grad")
	}
	if !approx(BackAzimuth(90), 270, 1e-9) {
		t.Error("back azimuth")
	}
	if !approx(RelativeBearing(350, 10), 20, 1e-9) {
		t.Errorf("relative bearing=%v", RelativeBearing(350, 10))
	}
}

func TestCompass(t *testing.T) {
	tests := []struct {
		b    float64
		want string
	}{
		{0, "N"}, {45, "NE"}, {90, "E"}, {225, "SW"},
	}
	for _, tc := range tests {
		if got := BearingToCompass8(tc.b); got != tc.want {
			t.Errorf("compass8(%v)=%q want %q", tc.b, got, tc.want)
		}
	}
	if BearingToCompass16(22.5) != "NNE" {
		t.Errorf("compass16(22.5)=%q", BearingToCompass16(22.5))
	}
	b, err := CompassToBearing("SW")
	if err != nil || !approx(b, 225, 1e-9) {
		t.Errorf("CompassToBearing SW=%v err=%v", b, err)
	}
	if _, err := CompassToBearing("XY"); err == nil {
		t.Error("expected error for bad compass point")
	}
}

// --- auxiliary latitudes & arcs --------------------------------------------

func TestArcs(t *testing.T) {
	if !approx(QuarterMeridian(WGS84), 10001965.729, 1e-3) {
		t.Errorf("quarter meridian=%v", QuarterMeridian(WGS84))
	}
	if !approx(LengthOfDegreeOfLatitude(0, WGS84), 110574.3, 0.5) {
		t.Errorf("deg lat @0=%v", LengthOfDegreeOfLatitude(0, WGS84))
	}
	if !approx(LengthOfDegreeOfLongitude(60, WGS84), 55800, 5) {
		t.Errorf("deg lon @60=%v", LengthOfDegreeOfLongitude(60, WGS84))
	}
	if !approx(MeridianArcBetween(0, 90, WGS84), QuarterMeridian(WGS84), 1e-6) {
		t.Errorf("arc between 0..90 mismatch")
	}
}

func TestAuxiliaryLatitudes(t *testing.T) {
	// All auxiliary latitudes agree with geodetic at 0 and 90 degrees.
	for _, lat := range []float64{0, 90, -90} {
		for _, f := range []func(float64, Ellipsoid) float64{
			ReducedLatitude, GeocentricLatitude, ConformalLatitude,
			AuthalicLatitude, RectifyingLatitude,
		} {
			if got := f(lat, WGS84); !approx(got, lat, 1e-4) {
				t.Errorf("aux lat at %v = %v", lat, got)
			}
		}
	}
	// Ordering at 45N: reduced > authalic > rectifying > conformal ~ geocentric.
	if !(ReducedLatitude(45, WGS84) > GeocentricLatitude(45, WGS84)) {
		t.Error("reduced should exceed geocentric at 45N")
	}
	// Geocentric roundtrip.
	if got := GeocentricToGeodeticLatitude(GeocentricLatitude(45, WGS84), WGS84); !approx(got, 45, 1e-9) {
		t.Errorf("geocentric roundtrip=%v", got)
	}
	// On a sphere all auxiliary latitudes equal the geodetic latitude.
	if got := ConformalLatitude(37, SphereWGS84); !approx(got, 37, 1e-9) {
		t.Errorf("conformal on sphere=%v", got)
	}
}

// --- grid scale & convergence ----------------------------------------------

func TestGridScaleConvergence(t *testing.T) {
	// Scale on the central meridian equals k0.
	if s := UTMPointScale(LatLon{45, 3}); !approx(s, 0.9996, 1e-5) {
		t.Errorf("scale on CM=%v", s)
	}
	// Convergence on the central meridian is zero.
	if c := UTMGridConvergence(LatLon{45, 3}); !approx(c, 0, 1e-4) {
		t.Errorf("convergence on CM=%v", c)
	}
	// Scale grows away from the central meridian.
	if s := UTMPointScale(LatLon{45, 6}); s <= 0.9996 {
		t.Errorf("scale off CM should exceed k0, got %v", s)
	}
}

// --- bounding box ----------------------------------------------------------

func TestBoundingBox(t *testing.T) {
	pts := []LatLon{{10, 20}, {12, 25}, {8, 22}}
	bb := BoundingBoxOfPoints(pts)
	if bb.Min.Lat != 8 || bb.Max.Lat != 12 || bb.Min.Lon != 20 || bb.Max.Lon != 25 {
		t.Errorf("bbox=%v", bb)
	}
	for _, p := range pts {
		if !bb.Contains(p) {
			t.Errorf("bbox should contain %v", p)
		}
	}
	c := bb.Center()
	if !approx(c.Lat, 10, 1e-9) || !approx(c.Lon, 22.5, 1e-9) {
		t.Errorf("center=%v", c)
	}
	around := BoundingBoxAround(LatLon{45, 9}, 1000)
	if !around.Contains(LatLon{45, 9}) {
		t.Error("bbox around should contain center")
	}
	if len(around.Corners()) != 4 {
		t.Error("corners")
	}
	other := NewBoundingBox(LatLon{11, 24}, LatLon{20, 30})
	if !bb.Intersects(other) {
		t.Error("boxes should intersect")
	}
	u := bb.Union(other)
	if u.Max.Lat != 20 || u.Max.Lon != 30 {
		t.Errorf("union=%v", u)
	}
}

// --- convenience methods ---------------------------------------------------

func TestLatLonMethods(t *testing.T) {
	p := LatLon{51.5, -0.1}
	q := LatLon{48.85, 2.35}
	if !approx(p.DistanceTo(q), HaversineDistance(p, q), 1e-9) {
		t.Error("DistanceTo")
	}
	if !approx(p.InitialBearingTo(q), InitialBearing(p, q), 1e-9) {
		t.Error("InitialBearingTo")
	}
	gd, err := p.GeodesicDistanceTo(q)
	if err != nil || gd <= 0 {
		t.Errorf("GeodesicDistanceTo=%v err=%v", gd, err)
	}
	if _, err := p.ToUTM(); err != nil {
		t.Errorf("ToUTM err=%v", err)
	}
	if !p.Valid() {
		t.Error("p should be valid")
	}
	if (LatLon{100, 0}).Valid() {
		t.Error("invalid latitude reported valid")
	}
}

// --- example ---------------------------------------------------------------

func ExampleHaversineDistance() {
	// Distance from London to Paris on the WGS-84 mean sphere.
	london := LatLon{51.5074, -0.1278}
	paris := LatLon{48.8566, 2.3522}
	km := MetersToKilometers(HaversineDistance(london, paris))
	fmt.Printf("%.1f km\n", km)
	// Output: 343.6 km
}

func ExampleVincentyInverse() {
	london := LatLon{51.5074, -0.1278}
	paris := LatLon{48.8566, 2.3522}
	inv, _ := VincentyInverse(london, paris, WGS84)
	fmt.Printf("%.0f m at %.1f°\n", inv.Distance, inv.InitialBearing)
	// Output: 343923 m at 148.0°
}
