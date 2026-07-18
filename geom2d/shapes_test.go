package geom2d

import (
	"math"
	"math/rand"
	"testing"
)

const shapesTol = 1e-9

func approxShapes(a, b float64) bool { return math.Abs(a-b) <= shapesTol }

func pointApproxTol(p, q Point2, tol float64) bool {
	return math.Abs(p.X-q.X) <= tol && math.Abs(p.Y-q.Y) <= tol
}

func TestPolygonAreaShoelace(t *testing.T) {
	tests := []struct {
		name string
		pts  []Point2
		want float64
	}{
		{"unit square", []Point2{{0, 0}, {1, 0}, {1, 1}, {0, 1}}, 1},
		{"3-4-5 triangle", []Point2{{0, 0}, {4, 0}, {0, 3}}, 6},
		{"rectangle 2x5", []Point2{{0, 0}, {5, 0}, {5, 2}, {0, 2}}, 10},
		{"clockwise square", []Point2{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, 1},
		{"L-shape", []Point2{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}, 3},
		{"degenerate line", []Point2{{0, 0}, {1, 1}, {2, 2}}, 0},
	}
	for _, tc := range tests {
		if got := PolygonArea(tc.pts); !approxShapes(got, tc.want) {
			t.Errorf("%s: PolygonArea = %v, want %v", tc.name, got, tc.want)
		}
	}
}

func TestPolygonSignedAreaWinding(t *testing.T) {
	ccw := []Point2{{0, 0}, {1, 0}, {1, 1}, {0, 1}}
	if a := PolygonSignedArea(ccw); a <= 0 {
		t.Errorf("CCW signed area should be positive, got %v", a)
	}
	cw := Polygon(ccw).Reverse()
	if a := cw.SignedArea(); a >= 0 {
		t.Errorf("CW signed area should be negative, got %v", a)
	}
}

func TestPolygonPerimeter(t *testing.T) {
	square := []Point2{{0, 0}, {3, 0}, {3, 3}, {0, 3}}
	if got := PolygonPerimeter(square); !approxShapes(got, 12) {
		t.Errorf("perimeter = %v, want 12", got)
	}
	tri := []Point2{{0, 0}, {4, 0}, {0, 3}}
	if got := PolygonPerimeter(tri); !approxShapes(got, 12) { // 3+4+5
		t.Errorf("triangle perimeter = %v, want 12", got)
	}
}

func TestPolygonCentroid(t *testing.T) {
	square := []Point2{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	c, ok := PolygonCentroid(square)
	if !ok || !pointApproxTol(c, Point2{1, 1}, shapesTol) {
		t.Errorf("square centroid = %v (ok=%v), want (1,1)", c, ok)
	}
	// Triangle centroid is the vertex average.
	tri := []Point2{{0, 0}, {6, 0}, {0, 9}}
	c, ok = PolygonCentroid(tri)
	if !ok || !pointApproxTol(c, Point2{2, 3}, shapesTol) {
		t.Errorf("triangle centroid = %v (ok=%v), want (2,3)", c, ok)
	}
}

func TestPointInPolygonRay(t *testing.T) {
	square := []Point2{{0, 0}, {4, 0}, {4, 4}, {0, 4}}
	inside := []Point2{{2, 2}, {1, 1}, {3, 3}, {0.1, 0.1}}
	outside := []Point2{{5, 5}, {-1, 2}, {2, 5}, {4.1, 2}}
	for _, p := range inside {
		if !PointInPolygon(p, square) {
			t.Errorf("point %v should be inside", p)
		}
	}
	for _, p := range outside {
		if PointInPolygon(p, square) {
			t.Errorf("point %v should be outside", p)
		}
	}
}

func TestPointInPolygonNonConvex(t *testing.T) {
	// A concave "arrowhead" polygon with a notch dipping to (2,1).
	poly := []Point2{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}
	inside := []Point2{{1, 1}, {3, 1}, {2, 0.5}}
	for _, p := range inside {
		if !PointInPolygon(p, poly) {
			t.Errorf("%v should be inside the concave polygon", p)
		}
	}
	if PointInPolygon(Point2{2, 3}, poly) {
		t.Errorf("(2,3) lies in the notch and should be outside")
	}
}

func TestWindingNumber(t *testing.T) {
	square := []Point2{{0, 0}, {4, 0}, {4, 4}, {0, 4}}
	if wn := WindingNumber(Point2{2, 2}, square); wn != 1 {
		t.Errorf("winding number inside = %d, want 1", wn)
	}
	if wn := WindingNumber(Point2{5, 5}, square); wn != 0 {
		t.Errorf("winding number outside = %d, want 0", wn)
	}
	if !PointInPolygonWinding(Point2{2, 2}, square) {
		t.Errorf("winding rule should report (2,2) inside")
	}
}

func TestIsConvexPolygon(t *testing.T) {
	convex := []Point2{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	if !IsConvexPolygon(convex) {
		t.Errorf("square should be convex")
	}
	concave := []Point2{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}
	if IsConvexPolygon(concave) {
		t.Errorf("arrowhead should be non-convex")
	}
}

func TestConvexHull(t *testing.T) {
	pts := []Point2{
		{0, 0}, {1, 1}, {2, 2}, {4, 0}, {0, 4}, {4, 4},
		{2, 2}, {3, 1}, {1, 3}, {2, 0}, {2, 4},
	}
	hull := ConvexHull(pts)
	// Hull is the 4x4 square, 4 vertices, area 16.
	if len(hull) != 4 {
		t.Fatalf("hull vertex count = %d, want 4: %v", len(hull), hull)
	}
	if a := hull.Area(); !approxShapes(a, 16) {
		t.Errorf("hull area = %v, want 16", a)
	}
	if !hull.IsCounterClockwise() {
		t.Errorf("hull should be counter-clockwise")
	}
	// Every input point must be inside or on the hull.
	for _, p := range pts {
		if !PointInPolygon(p, hull) && !OnPolygonBoundary(p, hull, 1e-9) {
			t.Errorf("point %v not contained by hull", p)
		}
	}
}

func TestConvexHullTriangle(t *testing.T) {
	pts := []Point2{{0, 0}, {5, 0}, {0, 5}, {1, 1}, {2, 1}}
	hull := ConvexHull(pts)
	if len(hull) != 3 {
		t.Fatalf("hull should have 3 vertices, got %d: %v", len(hull), hull)
	}
	if a := hull.Area(); !approxShapes(a, 12.5) {
		t.Errorf("hull area = %v, want 12.5", a)
	}
}

func TestTriangleAreaAndCentroid(t *testing.T) {
	tri := NewTriangle(Point2{0, 0}, Point2{4, 0}, Point2{0, 3})
	if a := tri.Area(); !approxShapes(a, 6) {
		t.Errorf("area = %v, want 6", a)
	}
	if p := tri.Perimeter(); !approxShapes(p, 12) {
		t.Errorf("perimeter = %v, want 12", p)
	}
	if c := tri.Centroid(); !pointApproxTol(c, Point2{4.0 / 3, 1}, shapesTol) {
		t.Errorf("centroid = %v, want (1.333,1)", c)
	}
}

func TestTriangleCircumcircle(t *testing.T) {
	// Right triangle: circumcenter at the hypotenuse midpoint, radius = half hyp.
	tri := NewTriangle(Point2{0, 0}, Point2{4, 0}, Point2{0, 3})
	c, ok := tri.Circumcircle()
	if !ok {
		t.Fatalf("circumcircle should exist")
	}
	if !pointApproxTol(c.Center, Point2{2, 1.5}, 1e-9) {
		t.Errorf("circumcenter = %v, want (2,1.5)", c.Center)
	}
	if !approxShapes(c.Radius, 2.5) {
		t.Errorf("circumradius = %v, want 2.5", c.Radius)
	}
	// Unit circle inscribed points.
	c2, ok := CircleFrom3Points(Point2{1, 0}, Point2{0, 1}, Point2{-1, 0})
	if !ok || !pointApproxTol(c2.Center, Point2{0, 0}, 1e-9) || !approxShapes(c2.Radius, 1) {
		t.Errorf("CircleFrom3Points = %+v (ok=%v), want unit circle", c2, ok)
	}
	_, ok = TriangleCircumcircle(Point2{0, 0}, Point2{1, 1}, Point2{2, 2})
	if ok {
		t.Errorf("collinear points must not yield a circumcircle")
	}
}

func TestTriangleIncircle(t *testing.T) {
	// 3-4-5 right triangle: inradius = (a+b-c)/2 = (3+4-5)/2 = 1.
	tri := NewTriangle(Point2{0, 0}, Point2{4, 0}, Point2{0, 3})
	ic, ok := tri.Incircle()
	if !ok {
		t.Fatalf("incircle should exist")
	}
	if !approxShapes(ic.Radius, 1) {
		t.Errorf("inradius = %v, want 1", ic.Radius)
	}
	if !pointApproxTol(ic.Center, Point2{1, 1}, 1e-9) {
		t.Errorf("incenter = %v, want (1,1)", ic.Center)
	}
}

func TestTriangleBarycentric(t *testing.T) {
	tri := NewTriangle(Point2{0, 0}, Point2{4, 0}, Point2{0, 4})
	// Centroid has barycentric (1/3,1/3,1/3).
	cen := tri.Centroid()
	u, v, w, ok := tri.Barycentric(cen)
	if !ok || !approxShapes(u, 1.0/3) || !approxShapes(v, 1.0/3) || !approxShapes(w, 1.0/3) {
		t.Errorf("barycentric of centroid = (%v,%v,%v), want thirds", u, v, w)
	}
	// Round trip.
	p := tri.PointFromBarycentric(0.2, 0.3, 0.5)
	u, v, w, _ = tri.Barycentric(p)
	if !approxShapes(u, 0.2) || !approxShapes(v, 0.3) || !approxShapes(w, 0.5) {
		t.Errorf("round-trip barycentric = (%v,%v,%v), want (0.2,0.3,0.5)", u, v, w)
	}
	if !tri.ContainsPoint(cen) {
		t.Errorf("centroid should be contained")
	}
	if tri.ContainsPoint(Point2{5, 5}) {
		t.Errorf("(5,5) should be outside the triangle")
	}
}

func TestTriangleOrthocenter(t *testing.T) {
	// Right triangle: orthocenter is at the right-angle vertex.
	tri := NewTriangle(Point2{0, 0}, Point2{4, 0}, Point2{0, 3})
	h, ok := tri.Orthocenter()
	if !ok || !pointApproxTol(h, Point2{0, 0}, 1e-9) {
		t.Errorf("orthocenter = %v (ok=%v), want (0,0)", h, ok)
	}
}

func TestCircleBasics(t *testing.T) {
	c := NewCircle(Point2{0, 0}, 2)
	if !approxShapes(c.Area(), math.Pi*4) {
		t.Errorf("area = %v, want 4pi", c.Area())
	}
	if !approxShapes(c.Circumference(), 4*math.Pi) {
		t.Errorf("circumference = %v, want 4pi", c.Circumference())
	}
	if !c.ContainsPoint(Point2{1, 1}) {
		t.Errorf("(1,1) should be inside radius-2 circle")
	}
	if c.ContainsPoint(Point2{2, 2}) {
		t.Errorf("(2,2) should be outside radius-2 circle")
	}
	if !c.OnBoundary(Point2{2, 0}, 1e-9) {
		t.Errorf("(2,0) should be on the boundary")
	}
}

func TestCircleLineIntersection(t *testing.T) {
	c := NewCircle(Point2{0, 0}, 1)
	// Horizontal line y=0 through the centre -> two points (-1,0),(1,0).
	pts := CircleLineIntersection(c, Line{A: Point2{-2, 0}, B: Point2{2, 0}})
	if len(pts) != 2 {
		t.Fatalf("expected 2 intersections, got %d: %v", len(pts), pts)
	}
	if !pointApproxTol(pts[0], Point2{-1, 0}, 1e-9) || !pointApproxTol(pts[1], Point2{1, 0}, 1e-9) {
		t.Errorf("intersections = %v, want (-1,0),(1,0)", pts)
	}
	// Tangent line y=1.
	tan := CircleLineIntersection(c, Line{A: Point2{-2, 1}, B: Point2{2, 1}})
	if len(tan) != 1 || !pointApproxTol(tan[0], Point2{0, 1}, 1e-7) {
		t.Errorf("tangent = %v, want single (0,1)", tan)
	}
	// Miss.
	if m := CircleLineIntersection(c, Line{A: Point2{-2, 2}, B: Point2{2, 2}}); len(m) != 0 {
		t.Errorf("expected no intersection, got %v", m)
	}
}

func TestCircleCircleIntersection(t *testing.T) {
	c1 := NewCircle(Point2{0, 0}, 2)
	c2 := NewCircle(Point2{3, 0}, 2)
	// Symmetric about x=1.5.
	pts, ok := CircleCircleIntersection(c1, c2)
	if !ok || len(pts) != 2 {
		t.Fatalf("expected 2 crossing points, got %d (ok=%v): %v", len(pts), ok, pts)
	}
	for _, p := range pts {
		if !approxShapes(Distance(p, c1.Center), 2) || !approxShapes(Distance(p, c2.Center), 2) {
			t.Errorf("intersection %v not on both circles", p)
		}
		if !approxShapes(p.X, 1.5) {
			t.Errorf("intersection x = %v, want 1.5", p.X)
		}
	}
	// Tangent externally.
	tan, ok := CircleCircleIntersection(NewCircle(Point2{0, 0}, 1), NewCircle(Point2{2, 0}, 1))
	if !ok || len(tan) != 1 || !pointApproxTol(tan[0], Point2{1, 0}, 1e-7) {
		t.Errorf("external tangent = %v (ok=%v), want single (1,0)", tan, ok)
	}
	// Separate.
	sep, ok := CircleCircleIntersection(NewCircle(Point2{0, 0}, 1), NewCircle(Point2{5, 0}, 1))
	if !ok || len(sep) != 0 {
		t.Errorf("separate circles = %v (ok=%v), want none", sep, ok)
	}
	// Coincident.
	if _, ok := CircleCircleIntersection(c1, c1); ok {
		t.Errorf("coincident circles should report ok=false")
	}
}

func TestBoundingBox(t *testing.T) {
	b, ok := BoundingBoxOf(Point2{1, 2}, Point2{-1, 5}, Point2{3, -2})
	if !ok {
		t.Fatalf("bounding box should exist")
	}
	if !pointApproxTol(b.Min, Point2{-1, -2}, 0) || !pointApproxTol(b.Max, Point2{3, 5}, 0) {
		t.Errorf("box = %+v, want min(-1,-2) max(3,5)", b)
	}
	if !approxShapes(b.Width(), 4) || !approxShapes(b.Height(), 7) || !approxShapes(b.Area(), 28) {
		t.Errorf("box dims wrong: w=%v h=%v a=%v", b.Width(), b.Height(), b.Area())
	}
	if !b.ContainsPoint(Point2{0, 0}) {
		t.Errorf("origin should be in the box")
	}
	if !pointApproxTol(b.Center(), Point2{1, 1.5}, shapesTol) {
		t.Errorf("center = %v, want (1,1.5)", b.Center())
	}
}

func TestBoundingBoxSetOps(t *testing.T) {
	a := NewBoundingBox(Point2{0, 0}, Point2{4, 4})
	c := NewBoundingBox(Point2{2, 2}, Point2{6, 6})
	if !a.Intersects(c) {
		t.Errorf("boxes should intersect")
	}
	inter, ok := a.Intersection(c)
	if !ok || !pointApproxTol(inter.Min, Point2{2, 2}, 0) || !pointApproxTol(inter.Max, Point2{4, 4}, 0) {
		t.Errorf("intersection = %+v (ok=%v), want min(2,2) max(4,4)", inter, ok)
	}
	uni := a.Union(c)
	if !pointApproxTol(uni.Min, Point2{0, 0}, 0) || !pointApproxTol(uni.Max, Point2{6, 6}, 0) {
		t.Errorf("union = %+v, want min(0,0) max(6,6)", uni)
	}
	d := NewBoundingBox(Point2{10, 10}, Point2{12, 12})
	if _, ok := a.Intersection(d); ok {
		t.Errorf("disjoint boxes should not intersect")
	}
}

func TestClipPolygon(t *testing.T) {
	// Clip a big square by a smaller square; result is the smaller square.
	subject := []Point2{{0, 0}, {4, 0}, {4, 4}, {0, 4}}
	clip := []Point2{{1, 1}, {3, 1}, {3, 3}, {1, 3}}
	res := ClipPolygon(subject, clip)
	if a := res.Area(); !approxShapes(a, 4) {
		t.Errorf("clipped area = %v, want 4", a)
	}
	// Clip a triangle poking outside the clip window.
	tri := []Point2{{-1, 0}, {5, 0}, {2, 5}}
	window := []Point2{{0, 0}, {4, 0}, {4, 4}, {0, 4}}
	clipped := ClipPolygon(tri, window)
	if len(clipped) < 3 {
		t.Fatalf("expected a non-degenerate clipped polygon, got %v", clipped)
	}
	// All resulting vertices must lie within the window (allowing boundary).
	for _, p := range clipped {
		if !window2Contains(p) {
			t.Errorf("clipped vertex %v outside window", p)
		}
	}
}

func window2Contains(p Point2) bool {
	return p.X >= -1e-9 && p.X <= 4+1e-9 && p.Y >= -1e-9 && p.Y <= 4+1e-9
}

func TestClipPolygonHalfPlane(t *testing.T) {
	// Keep the part of the unit square with x <= 1 (left of directed line
	// (1,0)->(1,1) is x<=1 for CCW-inside convention: inside is x<1... test half).
	square := []Point2{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	// Directed line from (1,0) to (1,2): inside (left) is x<1.
	res := ClipPolygonHalfPlane(square, Point2{1, 0}, Point2{1, 2})
	if a := res.Area(); !approxShapes(a, 2) {
		t.Errorf("half-plane clip area = %v, want 2", a)
	}
}

func TestConvexDiameter(t *testing.T) {
	tests := []struct {
		name string
		pts  []Point2
		want float64
	}{
		{"square", []Point2{{0, 0}, {1, 0}, {1, 1}, {0, 1}}, math.Sqrt2},
		{"segment", []Point2{{0, 0}, {3, 4}}, 5},
		{"triangle", []Point2{{0, 0}, {6, 0}, {0, 8}}, 10},
		{"collinear", []Point2{{0, 0}, {2, 0}, {5, 0}, {1, 0}}, 5},
	}
	for _, tc := range tests {
		if got := ConvexDiameter(tc.pts); !approxShapes(got, tc.want) {
			t.Errorf("%s: ConvexDiameter = %v, want %v", tc.name, got, tc.want)
		}
	}
	// Diameter points must reproduce the distance.
	p, q, d := ConvexDiameterPoints([]Point2{{0, 0}, {10, 0}, {5, 1}, {5, -1}})
	if !approxShapes(Distance(p, q), d) || !approxShapes(d, 10) {
		t.Errorf("diameter points %v,%v distance %v, want 10", p, q, d)
	}
}

func TestConvexDiameterMatchesBruteForce(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	for trial := 0; trial < 30; trial++ {
		n := 3 + rng.Intn(30)
		pts := make([]Point2, n)
		for i := range pts {
			pts[i] = Point2{X: rng.Float64() * 100, Y: rng.Float64() * 100}
		}
		var brute float64
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				if d := Distance(pts[i], pts[j]); d > brute {
					brute = d
				}
			}
		}
		if got := ConvexDiameter(pts); !approxShapes(got, brute) {
			t.Fatalf("trial %d: rotating calipers = %v, brute force = %v", trial, got, brute)
		}
	}
}

func BenchmarkConvexDiameterPoints(b *testing.B) {
	rng := rand.New(rand.NewSource(7))
	pts := make([]Point2, 2000)
	for i := range pts {
		pts[i] = Point2{X: rng.Float64() * 1000, Y: rng.Float64() * 1000}
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, _, _ = ConvexDiameterPoints(pts)
	}
}
