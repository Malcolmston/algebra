package fractal

import "testing"

func TestAffineMapApply(t *testing.T) {
	m := NewAffineMap(2, 0, 0, 3, 1, -1)
	got := m.Apply(Point2D{1, 1})
	if got.X != 3 || got.Y != 2 {
		t.Errorf("apply: got %+v want {3 2}", got)
	}
}

func TestChaosGameDeterministic(t *testing.T) {
	s := SierpinskiTriangleIFS()
	a := s.ChaosGame(1000, 42)
	b := s.ChaosGame(1000, 42)
	if len(a) != len(b) {
		t.Fatalf("length mismatch")
	}
	for i := range a {
		if a[i] != b[i] {
			t.Fatalf("same seed produced different output at %d", i)
		}
	}
	// A different seed generally differs.
	c := s.ChaosGame(1000, 43)
	diff := false
	for i := range a {
		if a[i] != c[i] {
			diff = true
			break
		}
	}
	if !diff {
		t.Errorf("different seeds produced identical output")
	}
}

func TestSierpinskiChaosGameInsideTriangle(t *testing.T) {
	// Attractor lies inside the triangle (0,0),(1,0),(0.5,0.5):
	// y>=0, y<=x, y<=1-x.
	pts := SierpinskiTriangleIFS().ChaosGame(20000, 7)
	const eps = 1e-9
	for _, p := range pts {
		if p.Y < -eps || p.Y > p.X+eps || p.Y > 1-p.X+eps {
			t.Fatalf("point %+v outside Sierpinski triangle", p)
		}
	}
}

func TestBarnsleyFernBounds(t *testing.T) {
	pts := BarnsleyFern().ChaosGame(50000, 99)
	for _, p := range pts {
		if p.X < -2.5 || p.X > 3.0 || p.Y < -0.1 || p.Y > 10.1 {
			t.Fatalf("fern point %+v out of expected bounds", p)
		}
	}
	// The fern is not degenerate: it spans a meaningful height.
	min, max := BoundingBox(pts)
	if max.Y-min.Y < 8 {
		t.Errorf("fern height too small: %v", max.Y-min.Y)
	}
}

func TestSierpinskiCarpetIFS(t *testing.T) {
	s := SierpinskiCarpetIFS()
	if len(s.Maps) != 8 {
		t.Fatalf("carpet should have 8 maps, got %d", len(s.Maps))
	}
	pts := s.ChaosGame(20000, 5)
	const eps = 1e-9
	for _, p := range pts {
		if p.X < -eps || p.X > 1+eps || p.Y < -eps || p.Y > 1+eps {
			t.Fatalf("carpet point %+v outside unit square", p)
		}
	}
}
