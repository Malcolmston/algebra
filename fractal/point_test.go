package fractal

import (
	"math"
	"testing"
)

func TestPoint2DOps(t *testing.T) {
	p := Point2D{1, 2}
	q := Point2D{3, 5}
	if p.Add(q) != (Point2D{4, 7}) {
		t.Errorf("Add wrong")
	}
	if q.Sub(p) != (Point2D{2, 3}) {
		t.Errorf("Sub wrong")
	}
	if p.Scale(2) != (Point2D{2, 4}) {
		t.Errorf("Scale wrong")
	}
	if p.Midpoint(q) != (Point2D{2, 3.5}) {
		t.Errorf("Midpoint wrong")
	}
	approx(t, p.Dist(q), math.Hypot(2, 3), 1e-12, "Dist")
	approx(t, Point2D{3, 4}.Norm(), 5, 1e-12, "Norm")
	if p.Lerp(q, 0) != p || p.Lerp(q, 1) != q {
		t.Errorf("Lerp endpoints wrong")
	}
}

func TestPoint2DRotate(t *testing.T) {
	r := Point2D{1, 0}.Rotate(math.Pi / 2)
	approx(t, r.X, 0, 1e-12, "rotate X")
	approx(t, r.Y, 1, 1e-12, "rotate Y")
}

func TestBoundingBox(t *testing.T) {
	pts := []Point2D{{1, -1}, {-2, 3}, {0, 0}}
	min, max := BoundingBox(pts)
	if min != (Point2D{-2, -1}) || max != (Point2D{1, 3}) {
		t.Errorf("bbox got min %+v max %+v", min, max)
	}
	min, max = BoundingBox(nil)
	if min != (Point2D{}) || max != (Point2D{}) {
		t.Errorf("empty bbox should be zero")
	}
}

func TestGrid(t *testing.T) {
	g := NewGrid(2, 2)
	g.Set(0, 0, 1)
	g.Set(1, 0, 2)
	g.Set(0, 1, 3)
	g.Set(1, 1, 4)
	if g.At(1, 1) != 4 {
		t.Errorf("At wrong")
	}
	min, max := g.MinMax()
	if min != 1 || max != 4 {
		t.Errorf("MinMax got %v %v", min, max)
	}
	approx(t, g.Mean(), 2.5, 1e-12, "Mean")
	n := g.Normalize()
	if n.At(0, 0) != 0 || n.At(1, 1) != 1 {
		t.Errorf("Normalize endpoints wrong: %v %v", n.At(0, 0), n.At(1, 1))
	}
	// Original unchanged.
	if g.At(0, 0) != 1 {
		t.Errorf("Normalize mutated receiver")
	}
	if g.CountFinite() != 4 {
		t.Errorf("CountFinite wrong")
	}
}
