package fractal

import (
	"math"
	"testing"
)

func TestKochCurveIteration1(t *testing.T) {
	pts := KochCurve(Point2D{0, 0}, Point2D{1, 0}, 1)
	want := []Point2D{
		{0, 0},
		{1.0 / 3, 0},
		{0.5, math.Sqrt(3) / 6},
		{2.0 / 3, 0},
		{1, 0},
	}
	if len(pts) != len(want) {
		t.Fatalf("point count: got %d want %d", len(pts), len(want))
	}
	for i := range want {
		if pts[i].Dist(want[i]) > 1e-12 {
			t.Errorf("point %d: got %+v want %+v", i, pts[i], want[i])
		}
	}
}

func TestKochCurvePointCount(t *testing.T) {
	// After n iterations there are 4^n + 1 points.
	for n := 0; n <= 4; n++ {
		pts := KochCurve(Point2D{0, 0}, Point2D{1, 0}, n)
		want := int(math.Pow(4, float64(n))) + 1
		if len(pts) != want {
			t.Errorf("n=%d: got %d points want %d", n, len(pts), want)
		}
	}
}

func TestKochCurveLength(t *testing.T) {
	// Length of the Koch curve after n iterations is (4/3)^n times the base.
	for n := 0; n <= 4; n++ {
		pts := KochCurve(Point2D{0, 0}, Point2D{1, 0}, n)
		total := 0.0
		for i := 0; i < len(pts)-1; i++ {
			total += pts[i].Dist(pts[i+1])
		}
		want := math.Pow(4.0/3.0, float64(n))
		approx(t, total, want, 1e-9, "koch length")
	}
}

func TestKochSnowflake(t *testing.T) {
	sf := KochSnowflake(Point2D{0, 0}, 1, 2)
	// Closed polyline.
	if sf[0] != sf[len(sf)-1] {
		t.Errorf("snowflake not closed")
	}
	// Perimeter = 3 * side * (4/3)^n, side = r*sqrt(3).
	total := 0.0
	for i := 0; i < len(sf)-1; i++ {
		total += sf[i].Dist(sf[i+1])
	}
	side := math.Sqrt(3.0)
	want := 3 * side * math.Pow(4.0/3.0, 2)
	approx(t, total, want, 1e-9, "snowflake perimeter")
}

func TestSierpinskiTriangle(t *testing.T) {
	a, b, c := Point2D{0, 0}, Point2D{1, 0}, Point2D{0.5, 1}
	for d := 0; d <= 5; d++ {
		tris := SierpinskiTriangle(a, b, c, d)
		want := int(math.Pow(3, float64(d)))
		if len(tris) != want {
			t.Errorf("depth %d: got %d triangles want %d", d, len(tris), want)
		}
	}
}

func TestCantorSet(t *testing.T) {
	ivs := CantorSet(0, 1, 2)
	want := []Interval{
		{0, 1.0 / 9}, {2.0 / 9, 3.0 / 9}, {6.0 / 9, 7.0 / 9}, {8.0 / 9, 1},
	}
	if len(ivs) != len(want) {
		t.Fatalf("interval count: got %d want %d", len(ivs), len(want))
	}
	var total float64
	for i := range want {
		approx(t, ivs[i].Start, want[i].Start, 1e-12, "cantor start")
		approx(t, ivs[i].End, want[i].End, 1e-12, "cantor end")
		total += ivs[i].Length()
	}
	// Total remaining length is (2/3)^2 = 4/9.
	approx(t, total, 4.0/9.0, 1e-12, "cantor total length")

	// General: n iterations -> 2^n intervals, total (2/3)^n.
	for n := 0; n <= 6; n++ {
		set := CantorSet(0, 1, n)
		if len(set) != int(math.Pow(2, float64(n))) {
			t.Errorf("cantor n=%d: got %d intervals", n, len(set))
		}
		var tl float64
		for _, iv := range set {
			tl += iv.Length()
		}
		approx(t, tl, math.Pow(2.0/3.0, float64(n)), 1e-12, "cantor length n")
	}
}
