package voronoi

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b float64) bool { return math.Abs(a-b) <= 1e-7 }

func ptApprox(p, q Point) bool { return p.ApproxEqual(q, 1e-7) }

// ---- Point arithmetic ----

func TestPointArithmetic(t *testing.T) {
	p := NewPoint(3, 4)
	q := NewPoint(1, 2)
	if got := p.Add(q); !ptApprox(got, Point{4, 6}) {
		t.Errorf("Add = %v", got)
	}
	if got := p.Sub(q); !ptApprox(got, Point{2, 2}) {
		t.Errorf("Sub = %v", got)
	}
	if got := p.Scale(2); !ptApprox(got, Point{6, 8}) {
		t.Errorf("Scale = %v", got)
	}
	if got := p.Dot(q); !approx(got, 11) {
		t.Errorf("Dot = %v", got)
	}
	if got := p.Cross(q); !approx(got, 2) {
		t.Errorf("Cross = %v", got)
	}
	if got := p.Norm(); !approx(got, 5) {
		t.Errorf("Norm = %v", got)
	}
	if got := p.NormSq(); !approx(got, 25) {
		t.Errorf("NormSq = %v", got)
	}
	if got := p.Distance(q); !approx(got, math.Sqrt(8)) {
		t.Errorf("Distance = %v", got)
	}
	if got := p.Midpoint(q); !ptApprox(got, Point{2, 3}) {
		t.Errorf("Midpoint = %v", got)
	}
	if got := p.Lerp(q, 0.5); !ptApprox(got, Point{2, 3}) {
		t.Errorf("Lerp = %v", got)
	}
	if got := p.Perp(); !ptApprox(got, Point{-4, 3}) {
		t.Errorf("Perp = %v", got)
	}
	if got := p.Neg(); !ptApprox(got, Point{-3, -4}) {
		t.Errorf("Neg = %v", got)
	}
}

func TestPointRotate(t *testing.T) {
	p := NewPoint(1, 0)
	if got := p.Rotate(math.Pi / 2); !ptApprox(got, Point{0, 1}) {
		t.Errorf("Rotate 90 = %v", got)
	}
	c := NewPoint(1, 1)
	if got := NewPoint(2, 1).RotateAbout(c, math.Pi/2); !ptApprox(got, Point{1, 2}) {
		t.Errorf("RotateAbout = %v", got)
	}
	if got := NewPoint(3, 0).Normalize(); !ptApprox(got, Point{1, 0}) {
		t.Errorf("Normalize = %v", got)
	}
}

func TestCentroidAndDedupe(t *testing.T) {
	pts := []Point{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	if got := Centroid(pts); !ptApprox(got, Point{1, 1}) {
		t.Errorf("Centroid = %v", got)
	}
	dup := []Point{{0, 0}, {1, 1}, {0, 0}, {1, 1}, {2, 2}}
	if got := DedupePoints(dup); len(got) != 3 {
		t.Errorf("DedupePoints len = %d, want 3", len(got))
	}
}

// ---- Orientation / predicates ----

func TestOrientation(t *testing.T) {
	tests := []struct {
		a, b, c Point
		want    int
	}{
		{Point{0, 0}, Point{1, 0}, Point{0, 1}, 1},
		{Point{0, 0}, Point{0, 1}, Point{1, 0}, -1},
		{Point{0, 0}, Point{1, 1}, Point{2, 2}, 0},
	}
	for _, tc := range tests {
		if got := Orientation(tc.a, tc.b, tc.c, tol); got != tc.want {
			t.Errorf("Orientation(%v,%v,%v) = %d, want %d", tc.a, tc.b, tc.c, got, tc.want)
		}
	}
	if !IsCCW(Point{0, 0}, Point{1, 0}, Point{0, 1}) {
		t.Error("IsCCW false")
	}
	if !IsCW(Point{0, 0}, Point{0, 1}, Point{1, 0}) {
		t.Error("IsCW false")
	}
	if !Collinear(Point{0, 0}, Point{1, 1}, Point{2, 2}, tol) {
		t.Error("Collinear false")
	}
}

func TestInCircle(t *testing.T) {
	a := Point{0, 0}
	b := Point{1, 0}
	c := Point{0, 1}
	// Point clearly inside the circumcircle (centre 0.5,0.5 r ~0.707).
	if !InCircleTest(a, b, c, Point{0.4, 0.4}) {
		t.Error("expected inside")
	}
	// Point clearly outside.
	if InCircleTest(a, b, c, Point{2, 2}) {
		t.Error("expected outside")
	}
	// Also works with clockwise input.
	if !InCircleTest(a, c, b, Point{0.4, 0.4}) {
		t.Error("expected inside (cw input)")
	}
}

func TestCircumcenter(t *testing.T) {
	center, err := Circumcenter(Point{0, 0}, Point{2, 0}, Point{0, 2})
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if !ptApprox(center, Point{1, 1}) {
		t.Errorf("Circumcenter = %v, want (1,1)", center)
	}
	r, _ := Circumradius(Point{0, 0}, Point{2, 0}, Point{0, 2})
	if !approx(r, math.Sqrt2) {
		t.Errorf("Circumradius = %v, want sqrt2", r)
	}
	if _, err := Circumcenter(Point{0, 0}, Point{1, 1}, Point{2, 2}); err != ErrCollinear {
		t.Errorf("expected ErrCollinear, got %v", err)
	}
}

func TestTriangleMetrics(t *testing.T) {
	a, b, c := Point{0, 0}, Point{4, 0}, Point{0, 3}
	if got := TriangleArea(a, b, c); !approx(got, 6) {
		t.Errorf("TriangleArea = %v", got)
	}
	if got := TrianglePerimeter(a, b, c); !approx(got, 12) {
		t.Errorf("TrianglePerimeter = %v", got)
	}
	if got := TriangleCentroid(a, b, c); !ptApprox(got, Point{4.0 / 3, 1}) {
		t.Errorf("TriangleCentroid = %v", got)
	}
	// 3-4-5 right triangle inradius = area/s = 6/6 = 1.
	if got := Inradius(a, b, c); !approx(got, 1) {
		t.Errorf("Inradius = %v", got)
	}
	if got := Incenter(a, b, c); !ptApprox(got, Point{1, 1}) {
		t.Errorf("Incenter = %v", got)
	}
}

func TestBarycentricAndPointInTriangle(t *testing.T) {
	a, b, c := Point{0, 0}, Point{4, 0}, Point{0, 4}
	p := Point{1, 1}
	u, v, w, ok := Barycentric(p, a, b, c)
	if !ok || !approx(u+v+w, 1) {
		t.Fatalf("barycentric sum: %v %v %v ok=%v", u, v, w, ok)
	}
	rec := a.Scale(u).Add(b.Scale(v)).Add(c.Scale(w))
	if !ptApprox(rec, p) {
		t.Errorf("reconstruct = %v", rec)
	}
	if !PointInTriangle(p, a, b, c, tol) {
		t.Error("expected p inside")
	}
	if PointInTriangle(Point{5, 5}, a, b, c, tol) {
		t.Error("expected outside")
	}
}

// ---- Circle ----

func TestCircleBasics(t *testing.T) {
	c := NewCircle(Point{0, 0}, 2)
	if !approx(c.Area(), math.Pi*4) {
		t.Errorf("Area = %v", c.Area())
	}
	if !approx(c.Circumference(), 4*math.Pi) {
		t.Errorf("Circumference = %v", c.Circumference())
	}
	if !c.Contains(Point{1, 1}, tol) {
		t.Error("Contains false")
	}
	if c.ContainsStrict(Point{2, 0}, tol) {
		t.Error("ContainsStrict should be false on boundary")
	}
	if !c.OnBoundary(Point{2, 0}, tol) {
		t.Error("OnBoundary false")
	}
	d := CircleFromDiameter(Point{0, 0}, Point{4, 0})
	if !ptApprox(d.Center, Point{2, 0}) || !approx(d.Radius, 2) {
		t.Errorf("CircleFromDiameter = %v", d)
	}
}

func TestMinimumEnclosingCircle(t *testing.T) {
	pts := []Point{{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}}
	c := MinimumEnclosingCircle(pts)
	// Smallest circle for a 2x2 square is centred at (1,1) with r = sqrt2.
	if !ptApprox(c.Center, Point{1, 1}) {
		t.Errorf("MEC center = %v", c.Center)
	}
	if !approx(c.Radius, math.Sqrt2) {
		t.Errorf("MEC radius = %v, want sqrt2", c.Radius)
	}
	for _, p := range pts {
		if !c.Contains(p, 1e-6) {
			t.Errorf("MEC does not contain %v", p)
		}
	}
}

// ---- Rect ----

func TestBoundingBox(t *testing.T) {
	pts := []Point{{1, 2}, {-3, 5}, {4, -1}}
	r, err := BoundingBox(pts)
	if err != nil {
		t.Fatal(err)
	}
	if !ptApprox(r.Min, Point{-3, -1}) || !ptApprox(r.Max, Point{4, 5}) {
		t.Errorf("BoundingBox = %v", r)
	}
	if !approx(r.Width(), 7) || !approx(r.Height(), 6) {
		t.Errorf("dims = %v x %v", r.Width(), r.Height())
	}
	if !approx(r.Area(), 42) {
		t.Errorf("Area = %v", r.Area())
	}
	if _, err := BoundingBox(nil); err != ErrEmpty {
		t.Errorf("expected ErrEmpty, got %v", err)
	}
}

// ---- Convex hull ----

func TestConvexHull(t *testing.T) {
	pts := []Point{
		{0, 0}, {1, 1}, {2, 2}, {2, 0}, {0, 2}, {1, 0.5}, {2, 2},
	}
	hull := ConvexHull(pts)
	// Expected square corners; interior/collinear points removed.
	want := map[Point]bool{{0, 0}: true, {2, 0}: true, {2, 2}: true, {0, 2}: true}
	if len(hull) != 4 {
		t.Fatalf("hull len = %d (%v), want 4", len(hull), hull)
	}
	for _, h := range hull {
		if !want[h] {
			t.Errorf("unexpected hull vertex %v", h)
		}
	}
	if !approx(ConvexHullArea(pts), 4) {
		t.Errorf("hull area = %v, want 4", ConvexHullArea(pts))
	}
	if !approx(ConvexHullPerimeter(pts), 8) {
		t.Errorf("hull perim = %v, want 8", ConvexHullPerimeter(pts))
	}
	// CCW winding.
	if PolygonSignedArea(hull) <= 0 {
		t.Error("hull not CCW")
	}
}

func TestGrahamMatchesMonotone(t *testing.T) {
	pts := []Point{{0, 0}, {3, 0}, {3, 3}, {0, 3}, {1, 1}, {2, 1}}
	m := ConvexHull(pts)
	g := GrahamHullByAngle(pts)
	if len(m) != len(g) {
		t.Fatalf("hull sizes differ: %d vs %d", len(m), len(g))
	}
	if !approx(PolygonArea(m), PolygonArea(g)) {
		t.Errorf("hull areas differ: %v vs %v", PolygonArea(m), PolygonArea(g))
	}
}

func TestDiameter(t *testing.T) {
	pts := []Point{{0, 0}, {3, 0}, {3, 4}, {0, 4}}
	d, _, _ := Diameter(pts)
	if !approx(d, 5) {
		t.Errorf("Diameter = %v, want 5", d)
	}
}

// ---- Polygon ----

func TestPolygon(t *testing.T) {
	sq := []Point{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	if !approx(PolygonArea(sq), 4) {
		t.Errorf("area = %v", PolygonArea(sq))
	}
	if !approx(PolygonPerimeter(sq), 8) {
		t.Errorf("perim = %v", PolygonPerimeter(sq))
	}
	if !ptApprox(PolygonCentroid(sq), Point{1, 1}) {
		t.Errorf("centroid = %v", PolygonCentroid(sq))
	}
	if !PolygonIsCCW(sq) {
		t.Error("expected CCW")
	}
	if !PointInPolygon(Point{1, 1}, sq) {
		t.Error("expected inside")
	}
	if PointInPolygon(Point{3, 3}, sq) {
		t.Error("expected outside")
	}
	if !PointInConvexPolygon(Point{1, 1}, sq, tol) {
		t.Error("expected inside convex")
	}
	if PolygonWindingNumber(Point{1, 1}, sq) == 0 {
		t.Error("winding number should be nonzero inside")
	}
}

// ---- Delaunay ----

func TestTriangulateSquare(t *testing.T) {
	pts := []Point{{0, 0}, {1, 0}, {1, 1}, {0, 1}}
	tri, err := Triangulate(pts)
	if err != nil {
		t.Fatal(err)
	}
	if tri.NumTriangles() != 2 {
		t.Errorf("NumTriangles = %d, want 2", tri.NumTriangles())
	}
	if !approx(tri.TotalArea(), 1) {
		t.Errorf("TotalArea = %v, want 1", tri.TotalArea())
	}
	if !tri.IsDelaunay() {
		t.Error("triangulation is not Delaunay")
	}
	// Euler characteristic V - E + F = 2.
	if tri.EulerCharacteristic() != 2 {
		t.Errorf("Euler = %d, want 2", tri.EulerCharacteristic())
	}
}

func TestTriangulateGrid(t *testing.T) {
	var pts []Point
	for x := 0; x < 4; x++ {
		for y := 0; y < 4; y++ {
			pts = append(pts, Point{float64(x), float64(y)})
		}
	}
	tri, err := Triangulate(pts)
	if err != nil {
		t.Fatal(err)
	}
	// Convex hull is a 3x3 square, so total area = 9.
	if !approx(tri.TotalArea(), 9) {
		t.Errorf("TotalArea = %v, want 9", tri.TotalArea())
	}
	if !tri.IsDelaunay() {
		t.Error("grid triangulation not Delaunay")
	}
	// For n points with h points on the hull boundary (including collinear
	// ones): triangles = 2n - 2 - h. The 4x4 grid has 16 points and 12 on the
	// boundary, giving 2*16 - 2 - 12 = 18 triangles.
	if tri.NumTriangles() != 18 {
		t.Errorf("NumTriangles = %d, want 18", tri.NumTriangles())
	}
}

func TestTriangulateErrors(t *testing.T) {
	if _, err := Triangulate([]Point{{0, 0}, {1, 1}}); err != ErrTooFewPoints {
		t.Errorf("expected ErrTooFewPoints, got %v", err)
	}
	if _, err := Triangulate([]Point{{0, 0}, {1, 1}, {2, 2}, {3, 3}}); err != ErrDegenerate {
		t.Errorf("expected ErrDegenerate, got %v", err)
	}
}

func TestLocate(t *testing.T) {
	pts := []Point{{0, 0}, {4, 0}, {4, 4}, {0, 4}, {2, 2}}
	tri, err := Triangulate(pts)
	if err != nil {
		t.Fatal(err)
	}
	if ti := tri.Locate(Point{1, 1}); ti < 0 {
		t.Error("expected point inside to be located")
	}
	if ti := tri.Locate(Point{10, 10}); ti >= 0 {
		t.Errorf("expected outside point to be unlocated, got %d", ti)
	}
	if !tri.Contains(Point{2, 1}) {
		t.Error("Contains false for interior point")
	}
}

func TestTriangulationTopology(t *testing.T) {
	pts := []Point{{0, 0}, {4, 0}, {4, 4}, {0, 4}, {2, 2}}
	tri, _ := Triangulate(pts)
	// Boundary edges should form the 4-cycle of the square.
	if len(tri.BoundaryEdges()) != 4 {
		t.Errorf("BoundaryEdges = %d, want 4", len(tri.BoundaryEdges()))
	}
	// Neighbour lists must be symmetric.
	adj := tri.Neighbors()
	for i, ns := range adj {
		for _, j := range ns {
			found := false
			for _, k := range adj[j] {
				if k == i {
					found = true
				}
			}
			if !found {
				t.Errorf("neighbour relation not symmetric: %d-%d", i, j)
			}
		}
	}
	// Centre vertex (2,2) connects to all four corners: degree 4.
	if d := tri.VertexDegree(4); d != 4 {
		t.Errorf("VertexDegree(center) = %d, want 4", d)
	}
}

// ---- Voronoi ----

func TestVoronoiFourPoints(t *testing.T) {
	// Symmetric square of sites; the single Voronoi vertex is the centre.
	pts := []Point{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	vd, err := Voronoi(pts)
	if err != nil {
		t.Fatal(err)
	}
	if vd.NumSites() != 4 {
		t.Errorf("NumSites = %d", vd.NumSites())
	}
	// Two triangles share the diagonal; their circumcentres both equal (1,1).
	for i, v := range vd.Vertices() {
		if !ptApprox(v, Point{1, 1}) {
			t.Errorf("Voronoi vertex %d = %v, want (1,1)", i, v)
		}
	}
	// There is one finite Voronoi edge (degenerate, length 0 here).
	if vd.NumEdges() != 1 {
		t.Errorf("NumEdges = %d, want 1", vd.NumEdges())
	}
	// All four cells are unbounded (all sites on the hull).
	for _, c := range vd.Cells() {
		if c.Bounded {
			t.Errorf("cell %d unexpectedly bounded", c.Site)
		}
	}
}

func TestVoronoiBoundedCell(t *testing.T) {
	// A central site surrounded by a ring gives one bounded cell.
	pts := []Point{
		{0, 0},
		{4, 0}, {8, 0}, {8, 4}, {8, 8}, {4, 8}, {0, 8}, {0, 4},
		{4, 4}, // centre
	}
	vd, err := Voronoi(pts)
	if err != nil {
		t.Fatal(err)
	}
	// Find the centre site's cell.
	var centre VoronoiCell
	for _, c := range vd.Cells() {
		if ptApprox(c.SiteXY, Point{4, 4}) {
			centre = c
		}
	}
	if !centre.Bounded {
		t.Fatal("centre cell should be bounded")
	}
	if centre.NumVertices() < 3 {
		t.Errorf("centre cell has %d vertices", centre.NumVertices())
	}
	// The cell is symmetric about (4,4), so its centroid is the site.
	if !centre.Centroid().ApproxEqual(Point{4, 4}, 1e-6) {
		t.Errorf("centre cell centroid = %v", centre.Centroid())
	}
	if centre.Area() <= 0 || math.IsInf(centre.Area(), 1) {
		t.Errorf("centre cell area = %v", centre.Area())
	}
}

func TestNearestSite(t *testing.T) {
	pts := []Point{{0, 0}, {10, 0}, {0, 10}, {10, 10}}
	vd, _ := Voronoi(pts)
	i, d := vd.NearestSite(Point{1, 1})
	if !ptApprox(vd.Sites()[i], Point{0, 0}) {
		t.Errorf("NearestSite = %v", vd.Sites()[i])
	}
	if !approx(d, math.Sqrt2) {
		t.Errorf("distance = %v", d)
	}
}

// ---- Queries ----

func TestNearestNeighbor(t *testing.T) {
	pts := []Point{{0, 0}, {1, 0}, {5, 5}, {0, 2}}
	i, d := NearestNeighbor(pts, Point{0, 0}, 0)
	if i != 1 || !approx(d, 1) {
		t.Errorf("NearestNeighbor = %d, %v", i, d)
	}
	k := KNearest(pts, Point{0, 0}, 2, 0)
	if len(k) != 2 || k[0] != 1 {
		t.Errorf("KNearest = %v", k)
	}
	all := AllNearestNeighbors(pts)
	if all[0] != 1 {
		t.Errorf("AllNearestNeighbors[0] = %d, want 1", all[0])
	}
}

func TestClosestPair(t *testing.T) {
	pts := []Point{{0, 0}, {10, 10}, {11, 11}, {5, 5}, {0, 1}}
	a, b, d := ClosestPair(pts)
	// The closest pair is (0,0)-(0,1) at distance 1.
	if !approx(d, 1) {
		t.Errorf("ClosestPair dist = %v, want 1", d)
	}
	got := map[Point]bool{pts[a]: true, pts[b]: true}
	if !got[Point{0, 0}] || !got[Point{0, 1}] {
		t.Errorf("ClosestPair = %v,%v", pts[a], pts[b])
	}
}

func TestClosestPairMatchesBruteForce(t *testing.T) {
	pts := []Point{
		{3, 1}, {7, 8}, {2, 2}, {9, 3}, {5, 5}, {1, 9}, {8, 1}, {4, 6}, {6, 2}, {0, 0},
	}
	_, _, fast := ClosestPair(pts)
	brute := math.Inf(1)
	for i := range pts {
		for j := i + 1; j < len(pts); j++ {
			if d := pts[i].Distance(pts[j]); d < brute {
				brute = d
			}
		}
	}
	if !approx(fast, brute) {
		t.Errorf("ClosestPair = %v, brute = %v", fast, brute)
	}
}

func TestLargestEmptyCircle(t *testing.T) {
	// Four corners of a 10x10 square; the largest empty circle sits at the
	// centre with radius = half the diagonal.
	pts := []Point{{0, 0}, {10, 0}, {10, 10}, {0, 10}}
	c, err := LargestEmptyCircle(pts)
	if err != nil {
		t.Fatal(err)
	}
	if !ptApprox(c.Center, Point{5, 5}) {
		t.Errorf("LEC center = %v, want (5,5)", c.Center)
	}
	if !approx(c.Radius, math.Sqrt(50)) {
		t.Errorf("LEC radius = %v, want sqrt50", c.Radius)
	}
}

func TestGabrielAndRNG(t *testing.T) {
	pts := []Point{{0, 0}, {2, 0}, {1, 0.1}}
	g := GabrielEdges(pts)
	// The edge (0,0)-(2,0) has (1,0.1) inside its diametral circle, so it is
	// excluded from the Gabriel graph.
	for _, e := range g {
		if e.Equal(Edge{0, 1}) {
			t.Error("edge 0-1 should not be a Gabriel edge")
		}
	}
	rng := RelativeNeighborhoodEdges(pts)
	if len(rng) == 0 {
		t.Error("expected some RNG edges")
	}
}

func TestEuclideanMST(t *testing.T) {
	pts := []Point{{0, 0}, {1, 0}, {2, 0}, {3, 0}}
	mst := EuclideanMST(pts)
	if len(mst) != 3 {
		t.Fatalf("MST edges = %d, want 3", len(mst))
	}
	if !approx(TotalEdgeLength(pts, mst), 3) {
		t.Errorf("MST length = %v, want 3", TotalEdgeLength(pts, mst))
	}
}

// ---- Alpha shapes ----

func TestAlphaShape(t *testing.T) {
	pts := []Point{{0, 0}, {1, 0}, {1, 1}, {0, 1}}
	// Large alpha => full convex hull; its boundary is the 4 square edges.
	as, err := AlphaShapeOf(pts, 100)
	if err != nil {
		t.Fatal(err)
	}
	if as.NumEdges() != 4 {
		t.Errorf("large-alpha boundary edges = %d, want 4", as.NumEdges())
	}
	if !approx(as.Area(), 1) {
		t.Errorf("alpha area = %v, want 1", as.Area())
	}
	if !approx(as.Perimeter(), 4) {
		t.Errorf("alpha perimeter = %v, want 4", as.Perimeter())
	}
	// Tiny alpha => no triangle qualifies (circumradius ~0.707 > alpha).
	small, _ := AlphaShapeOf(pts, 0.1)
	if small.NumTriangles() != 0 {
		t.Errorf("small-alpha triangles = %d, want 0", small.NumTriangles())
	}
}

func TestCircumradiusSpectrum(t *testing.T) {
	pts := []Point{{0, 0}, {1, 0}, {1, 1}, {0, 1}}
	tri, _ := Triangulate(pts)
	spec := CircumradiusSpectrum(tri)
	if len(spec) != 2 {
		t.Fatalf("spectrum len = %d, want 2", len(spec))
	}
	for _, r := range spec {
		if !approx(r, math.Sqrt2/2) {
			t.Errorf("circumradius = %v, want sqrt2/2", r)
		}
	}
}

// ---- Topology helpers ----

func TestEdgeAndTriangle(t *testing.T) {
	e := NewEdge(3, 1)
	if c := e.Canonical(); c.A != 1 || c.B != 3 {
		t.Errorf("Canonical = %v", c)
	}
	if !e.Equal(NewEdge(1, 3)) {
		t.Error("Equal false")
	}
	if e.Other(3) != 1 {
		t.Errorf("Other = %d", e.Other(3))
	}
	tr := NewTriangle(0, 1, 2)
	if !tr.HasEdge(NewEdge(1, 2)) {
		t.Error("HasEdge false")
	}
	if tr.OppositeVertex(NewEdge(0, 1)) != 2 {
		t.Errorf("OppositeVertex = %d", tr.OppositeVertex(NewEdge(0, 1)))
	}
	if _, ok := tr.SharesEdge(NewTriangle(1, 2, 3)); !ok {
		t.Error("SharesEdge false")
	}
}

// ---- Runnable example ----

func ExampleVoronoi() {
	sites := []Point{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	vd, err := Voronoi(sites)
	if err != nil {
		panic(err)
	}
	fmt.Println("sites:", vd.NumSites())
	fmt.Println("vertex:", vd.Vertices()[0])
	i, _ := vd.NearestSite(Point{0.3, 0.3})
	fmt.Println("nearest site to (0.3,0.3):", vd.Sites()[i])
	// Output:
	// sites: 4
	// vertex: (1, 1)
	// nearest site to (0.3,0.3): (0, 0)
}

func ExampleTriangulate() {
	pts := []Point{{0, 0}, {1, 0}, {1, 1}, {0, 1}}
	tri, err := Triangulate(pts)
	if err != nil {
		panic(err)
	}
	fmt.Println("triangles:", tri.NumTriangles())
	fmt.Printf("total area: %.1f\n", tri.TotalArea())
	fmt.Println("delaunay:", tri.IsDelaunay())
	// Output:
	// triangles: 2
	// total area: 1.0
	// delaunay: true
}

func ExampleConvexHull() {
	pts := []Point{{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 1}}
	hull := ConvexHull(pts)
	fmt.Println("hull size:", len(hull))
	fmt.Printf("hull area: %.1f\n", PolygonArea(hull))
	// Output:
	// hull size: 4
	// hull area: 4.0
}
