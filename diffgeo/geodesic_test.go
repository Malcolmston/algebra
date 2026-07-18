package diffgeo

import (
	"math"
	"testing"
)

func TestGeodesicPlaneIsStraight(t *testing.T) {
	s := SurfacePlane(Vec3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0})
	start := GeodesicState{U: 0, V: 0, DU: 1, DV: 2}
	path := GeodesicPath(s, start, 3.0, 100)
	end := path[len(path)-1]
	// On a plane geodesics are straight lines: position = t*velocity.
	if math.Abs(end.U-3.0) > 1e-6 || math.Abs(end.V-6.0) > 1e-6 {
		t.Errorf("plane geodesic end = (%v,%v) want (3,6)", end.U, end.V)
	}
	// Velocity is unchanged.
	if math.Abs(end.DU-1) > 1e-6 || math.Abs(end.DV-2) > 1e-6 {
		t.Errorf("plane geodesic velocity drifted: (%v,%v)", end.DU, end.DV)
	}
}

func TestGeodesicSphereEquator(t *testing.T) {
	s := Sphere(1.0)
	// Start on the equator (v=0) heading along it (DV=0): must stay on equator.
	start := GeodesicState{U: 0, V: 0, DU: 1, DV: 0}
	path := GeodesicPath(s, start, 2.0, 400)
	for _, st := range path {
		if math.Abs(st.V) > 1e-4 {
			t.Errorf("equator geodesic left v=0: v=%v at u=%v", st.V, st.U)
		}
	}
}

func TestGeodesicGreatCircleLength(t *testing.T) {
	R := 2.0
	s := Sphere(R)
	// Unit-speed motion along the equator: arc length ~= R * (parameter span).
	start := GeodesicState{U: 0, V: 0, DU: 1.0 / R, DV: 0}
	path := GeodesicPath(s, start, math.Pi*R, 500) // quarter+ turn worth
	got := GeodesicLength(s, path)
	want := math.Pi * R // unit speed over parameter length pi*R
	if math.Abs(got-want) > 1e-2 {
		t.Errorf("great-circle length = %v want %v", got, want)
	}
}

func TestGeodesicAccelerationPlaneZero(t *testing.T) {
	s := SurfacePlane(Vec3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0})
	au, av := GeodesicAcceleration(s, GeodesicState{U: 1, V: 1, DU: 0.5, DV: -0.7})
	if math.Abs(au) > 1e-4 || math.Abs(av) > 1e-4 {
		t.Errorf("plane geodesic acceleration = (%v,%v) want 0", au, av)
	}
}

func TestParallelTransportPlaneConstant(t *testing.T) {
	s := SurfacePlane(Vec3{}, Vec3{1, 0, 0}, Vec3{0, 1, 0})
	start := GeodesicState{U: 0, V: 0, DU: 1, DV: 1}
	_, vecs := ParallelTransport(s, start, 0.3, -0.8, 2.0, 100)
	final := vecs[len(vecs)-1]
	if math.Abs(final[0]-0.3) > 1e-6 || math.Abs(final[1]+0.8) > 1e-6 {
		t.Errorf("plane transport changed vector: %v want (0.3,-0.8)", final)
	}
}

func TestParallelTransportPreservesNorm(t *testing.T) {
	s := Sphere(1.5)
	start := GeodesicState{U: 0.2, V: 0.4, DU: 1.0, DV: 0.3}
	path, vecs := ParallelTransport(s, start, 1.0, -0.5, 1.5, 600)
	metricNorm := func(st GeodesicState, w [2]float64) float64 {
		I := FirstFundamental(s, st.U, st.V)
		return math.Sqrt(I.E*w[0]*w[0] + 2*I.F*w[0]*w[1] + I.G*w[1]*w[1])
	}
	n0 := metricNorm(path[0], vecs[0])
	nEnd := metricNorm(path[len(path)-1], vecs[len(vecs)-1])
	if math.Abs(n0-nEnd) > 1e-2 {
		t.Errorf("transport did not preserve metric norm: %v -> %v", n0, nEnd)
	}
}

func BenchmarkParallelTransportSphere(b *testing.B) {
	s := Sphere(2.0)
	start := GeodesicState{U: 0.3, V: 0.5, DU: 1.0, DV: 0.2}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _ = ParallelTransport(s, start, 1.0, 0.0, 2.0, 500)
	}
}
