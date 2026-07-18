package fractal

import (
	"math"
	"testing"
)

func TestFitLine(t *testing.T) {
	// y = 2x + 1 exactly.
	xs := []float64{0, 1, 2, 3, 4}
	ys := []float64{1, 3, 5, 7, 9}
	m, b := FitLine(xs, ys)
	approx(t, m, 2, 1e-12, "slope")
	approx(t, b, 1, 1e-12, "intercept")
	approx(t, FitSlope(xs, ys), 2, 1e-12, "FitSlope")
}

func TestHausdorffDimensionSelfSimilar(t *testing.T) {
	cases := []struct {
		n     int
		scale float64
		want  float64
	}{
		{3, 0.5, math.Log(3) / math.Log(2)},       // Sierpinski triangle ~1.585
		{4, 1.0 / 3, math.Log(4) / math.Log(3)},   // Koch curve ~1.262
		{2, 1.0 / 3, math.Log(2) / math.Log(3)},   // Cantor set ~0.631
		{8, 1.0 / 3, math.Log(8) / math.Log(3)},   // Sierpinski carpet ~1.893
		{20, 1.0 / 3, math.Log(20) / math.Log(3)}, // Menger-like
	}
	for _, c := range cases {
		got := HausdorffDimensionSelfSimilar(c.n, c.scale)
		approx(t, got, c.want, 1e-12, "hausdorff")
	}
}

func TestBoxCountingDimensionFromCounts(t *testing.T) {
	// Exact power law: count = size^(-D). The recovered slope must equal D.
	D := math.Log(3) / math.Log(2)
	sizes := []float64{0.5, 0.25, 0.125, 0.0625, 0.03125}
	counts := make([]float64, len(sizes))
	for i, s := range sizes {
		counts[i] = math.Pow(1/s, D)
	}
	got := BoxCountingDimensionFromCounts(sizes, counts)
	approx(t, got, D, 1e-9, "boxdim from counts")
}

func TestBoxCount(t *testing.T) {
	// Four points in distinct unit boxes.
	pts := []Point2D{{0.1, 0.1}, {1.1, 0.1}, {0.1, 1.1}, {1.1, 1.1}}
	if n := BoxCount(pts, 1.0); n != 4 {
		t.Errorf("box size 1: got %d want 4", n)
	}
	// With a box big enough all four fall in one box.
	if n := BoxCount(pts, 10.0); n != 1 {
		t.Errorf("box size 10: got %d want 1", n)
	}
}

func TestBoxCountingDimensionLine(t *testing.T) {
	// A dense straight segment has box-counting dimension close to 1.
	var pts []Point2D
	for i := 0; i <= 4000; i++ {
		x := float64(i) / 4000
		pts = append(pts, Point2D{x, 0})
	}
	sizes := []float64{0.1, 0.05, 0.025, 0.0125, 0.00625}
	d := BoxCountingDimension(pts, sizes)
	approx(t, d, 1.0, 0.1, "line box dimension")
}

func TestBoxCountingDimensionSierpinski(t *testing.T) {
	// The Sierpinski triangle attractor should have box-counting dimension
	// near the theoretical log2(3) ~ 1.585.
	pts := SierpinskiTriangleIFS().ChaosGame(60000, 12345)
	sizes := []float64{0.1, 0.05, 0.025, 0.0125}
	d := BoxCountingDimension(pts, sizes)
	want := math.Log(3) / math.Log(2)
	approx(t, d, want, 0.15, "sierpinski box dimension")
}
