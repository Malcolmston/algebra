package meshgen

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b float64) bool { return math.Abs(a-b) <= 1e-7 }

// ---- Vec2 ----

func TestVec2Arithmetic(t *testing.T) {
	p := NewVec2(3, 4)
	q := NewVec2(1, 2)
	tests := []struct {
		name string
		got  float64
		want float64
	}{
		{"Dot", p.Dot(q), 11},
		{"Cross", p.Cross(q), 2},
		{"Norm", p.Norm(), 5},
		{"NormSq", p.NormSq(), 25},
		{"Distance", p.Distance(q), math.Sqrt(8)},
	}
	for _, tc := range tests {
		if !approx(tc.got, tc.want) {
			t.Errorf("%s = %v, want %v", tc.name, tc.got, tc.want)
		}
	}
	if got := p.Add(q); !got.ApproxEqual(Vec2{4, 6}, tol) {
		t.Errorf("Add = %v", got)
	}
	if got := p.Sub(q); !got.ApproxEqual(Vec2{2, 2}, tol) {
		t.Errorf("Sub = %v", got)
	}
	if got := p.Normalize().Norm(); !approx(got, 1) {
		t.Errorf("Normalize norm = %v", got)
	}
	if got := p.Perp(); !got.ApproxEqual(Vec2{-4, 3}, tol) {
		t.Errorf("Perp = %v", got)
	}
	if got := NewVec2(1, 0).Rotate(math.Pi / 2); !got.ApproxEqual(Vec2{0, 1}, 1e-9) {
		t.Errorf("Rotate = %v", got)
	}
}

func TestCentroidAndBounds(t *testing.T) {
	pts := []Vec2{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	if got := CentroidVec2(pts); !got.ApproxEqual(Vec2{1, 1}, tol) {
		t.Errorf("centroid = %v", got)
	}
	min, max := BoundingBox2(pts)
	if !min.ApproxEqual(Vec2{0, 0}, tol) || !max.ApproxEqual(Vec2{2, 2}, tol) {
		t.Errorf("bounds = %v %v", min, max)
	}
}

// ---- Vec3 ----

func TestVec3Arithmetic(t *testing.T) {
	p := NewVec3(1, 2, 2)
	if !approx(p.Norm(), 3) {
		t.Errorf("Norm = %v", p.Norm())
	}
	x := NewVec3(1, 0, 0)
	y := NewVec3(0, 1, 0)
	if got := x.Cross(y); !got.ApproxEqual(Vec3{0, 0, 1}, tol) {
		t.Errorf("Cross = %v", got)
	}
	if got := x.AngleBetween(y); !approx(got, math.Pi/2) {
		t.Errorf("AngleBetween = %v", got)
	}
	if got := x.TripleProduct(y, NewVec3(0, 0, 1)); !approx(got, 1) {
		t.Errorf("TripleProduct = %v", got)
	}
}

// ---- Predicates ----

func TestPredicates(t *testing.T) {
	a := NewVec2(0, 0)
	b := NewVec2(1, 0)
	c := NewVec2(0, 1)
	if Orient2D(a, b, c) <= 0 {
		t.Errorf("Orient2D should be positive (CCW)")
	}
	if !IsCCW(a, b, c) {
		t.Errorf("IsCCW false")
	}
	if got := TriangleArea(a, b, c); !approx(got, 0.5) {
		t.Errorf("TriangleArea = %v", got)
	}
	if InCircle(a, b, c, NewVec2(0.2, 0.2)) <= 0 {
		t.Errorf("point should be inside circle")
	}
	if InCircle(a, b, c, NewVec2(5, 5)) >= 0 {
		t.Errorf("point should be outside circle")
	}
	center, r, err := Circumcircle(a, b, c)
	if err != nil {
		t.Fatalf("Circumcircle err: %v", err)
	}
	if !center.ApproxEqual(Vec2{0.5, 0.5}, 1e-9) || !approx(r, math.Sqrt(0.5)) {
		t.Errorf("circumcircle = %v r=%v", center, r)
	}
	if _, _, err := Circumcircle(a, b, NewVec2(2, 0)); err == nil {
		t.Errorf("expected ErrCollinear")
	}
	u, v, w, ok := Barycentric(NewVec2(0.25, 0.25), a, b, c)
	if !ok || !approx(u, 0.5) || !approx(v, 0.25) || !approx(w, 0.25) {
		t.Errorf("barycentric = %v %v %v", u, v, w)
	}
	if !PointInTriangle(NewVec2(0.25, 0.25), a, b, c, 1e-12) {
		t.Errorf("point should be inside triangle")
	}
	if PointInTriangle(NewVec2(1, 1), a, b, c, 1e-12) {
		t.Errorf("point should be outside triangle")
	}
}

func TestPolygonPredicates(t *testing.T) {
	sq := []Vec2{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	if got := PolygonArea(sq); !approx(got, 4) {
		t.Errorf("PolygonArea = %v", got)
	}
	if !PolygonIsCCW(sq) {
		t.Errorf("square should be CCW")
	}
	if got := PolygonPerimeter(sq); !approx(got, 8) {
		t.Errorf("perimeter = %v", got)
	}
	if got := PolygonCentroid(sq); !got.ApproxEqual(Vec2{1, 1}, 1e-9) {
		t.Errorf("centroid = %v", got)
	}
	if !PointInPolygon(NewVec2(1, 1), sq) {
		t.Errorf("center should be inside polygon")
	}
	if PointInPolygon(NewVec2(3, 3), sq) {
		t.Errorf("point should be outside polygon")
	}
}

func TestSegmentPredicates(t *testing.T) {
	if !SegmentsProperlyIntersect(NewVec2(0, 0), NewVec2(2, 2), NewVec2(0, 2), NewVec2(2, 0)) {
		t.Errorf("segments should cross")
	}
	if SegmentsProperlyIntersect(NewVec2(0, 0), NewVec2(1, 0), NewVec2(0, 1), NewVec2(1, 1)) {
		t.Errorf("segments should not cross")
	}
	p, ok := SegmentIntersection(NewVec2(0, 0), NewVec2(2, 2), NewVec2(0, 2), NewVec2(2, 0))
	if !ok || !p.ApproxEqual(Vec2{1, 1}, 1e-9) {
		t.Errorf("intersection = %v", p)
	}
	if got := DistancePointSegment(NewVec2(0, 1), NewVec2(0, 0), NewVec2(2, 0)); !approx(got, 1) {
		t.Errorf("distance = %v", got)
	}
}

func TestOrient3D(t *testing.T) {
	a := NewVec3(0, 0, 0)
	b := NewVec3(1, 0, 0)
	c := NewVec3(0, 1, 0)
	d := NewVec3(0, 0, 1)
	if got := TetraVolume(a, b, c, d); !approx(got, 1.0/6) {
		t.Errorf("TetraVolume = %v", got)
	}
	if got := TriangleArea3(a, b, c); !approx(got, 0.5) {
		t.Errorf("TriangleArea3 = %v", got)
	}
}

// ---- Delaunay ----

func TestDelaunaySquareCenter(t *testing.T) {
	pts := []Vec2{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}}
	tri, err := Triangulate(pts)
	if err != nil {
		t.Fatalf("Triangulate: %v", err)
	}
	if got := tri.NumTriangles(); got != 4 {
		t.Errorf("NumTriangles = %d, want 4", got)
	}
	if !tri.IsDelaunay(1e-9) {
		t.Errorf("triangulation not Delaunay")
	}
	m := tri.Mesh()
	if !approx(m.TotalArea(), 1) {
		t.Errorf("total area = %v, want 1", m.TotalArea())
	}
	if err := m.Validate(); err != nil {
		t.Errorf("validate: %v", err)
	}
}

func TestDelaunayGrid(t *testing.T) {
	pts := GridPoints(4, 4, 0, 0, 1, 1)
	m, err := DelaunayMesh(pts)
	if err != nil {
		t.Fatalf("DelaunayMesh: %v", err)
	}
	// 3x3 = 9 cells, each 2 triangles => 18 triangles covering a 3x3 area.
	if got := m.NumTriangles(); got != 18 {
		t.Errorf("NumTriangles = %d, want 18", got)
	}
	if !approx(m.TotalArea(), 9) {
		t.Errorf("area = %v, want 9", m.TotalArea())
	}
	if m.EulerCharacteristic() != 1 {
		t.Errorf("Euler = %d, want 1", m.EulerCharacteristic())
	}
}

func TestConvexHull(t *testing.T) {
	pts := []Vec2{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}, {0.2, 0.3}}
	hull := ConvexHull(pts)
	if len(hull) != 4 {
		t.Fatalf("hull size = %d, want 4", len(hull))
	}
	if !approx(PolygonArea(hull), 1) {
		t.Errorf("hull area = %v, want 1", PolygonArea(hull))
	}
	if !PolygonIsCCW(hull) {
		t.Errorf("hull should be CCW")
	}
}

// ---- CDT ----

func TestConstrainedEdge(t *testing.T) {
	// Convex quad, non-cyclic.
	pts := []Vec2{{0, 0}, {2, 0}, {3, 3}, {0, 2}}
	tri, err := Triangulate(pts)
	if err != nil {
		t.Fatalf("Triangulate: %v", err)
	}
	if err := tri.InsertConstraint(1, 3); err != nil {
		t.Fatalf("InsertConstraint: %v", err)
	}
	if !tri.hasEdge(1, 3) {
		t.Errorf("constraint edge 1-3 not present")
	}
	m := tri.Mesh()
	if !approx(m.TotalArea(), PolygonArea(pts)) {
		t.Errorf("area = %v, want %v", m.TotalArea(), PolygonArea(pts))
	}
}

func TestTriangulateConstrained(t *testing.T) {
	pts := []Vec2{{0, 0}, {4, 0}, {4, 4}, {0, 4}, {3, 1}}
	segs := [][2]int{{0, 2}}
	tri, err := TriangulateConstrained(pts, segs)
	if err != nil {
		t.Fatalf("TriangulateConstrained: %v", err)
	}
	if !tri.hasEdge(0, 2) {
		t.Errorf("constraint 0-2 not recovered")
	}
}

// ---- Refinement ----

func TestRefineImprovesQuality(t *testing.T) {
	// A random cloud with slivers: guarded Ruppert refinement must improve the
	// worst angle without ever degrading it, and must not change the covered
	// area (Steiner points are interior or on collinear boundary midpoints).
	pts := RandomPoints(40, 11, Vec2{0, 0}, Vec2{5, 5})
	tri, err := Triangulate(pts)
	if err != nil {
		t.Fatalf("Triangulate: %v", err)
	}
	before := tri.Mesh()
	beforeMin := before.MinAngleDeg()
	beforeArea := before.TotalArea()

	added := tri.Refine(25, 3000)
	after := tri.Mesh()

	if added == 0 {
		t.Fatalf("no Steiner points added")
	}
	if after.MinAngleDeg() < beforeMin-1e-9 {
		t.Errorf("min angle degraded: before %v after %v", beforeMin, after.MinAngleDeg())
	}
	if after.MinAngleDeg() <= beforeMin {
		t.Errorf("min angle did not improve: before %v after %v", beforeMin, after.MinAngleDeg())
	}
	if math.Abs(after.TotalArea()-beforeArea) > 0.02 {
		t.Errorf("area changed: before %v after %v", beforeArea, after.TotalArea())
	}
	if err := after.Validate(); err != nil {
		t.Errorf("refined mesh invalid: %v", err)
	}
}

// ---- Quality ----

func TestTriangleQualityEquilateral(t *testing.T) {
	a := NewVec2(0, 0)
	b := NewVec2(1, 0)
	c := NewVec2(0.5, math.Sqrt(3)/2)
	q := TriangleQualityOf(a, b, c)
	if !approx(q.MinAngle, math.Pi/3) || !approx(q.MaxAngle, math.Pi/3) {
		t.Errorf("angles = %v", q.Angles)
	}
	if !approx(q.RadiusRatio, 1) {
		t.Errorf("radius ratio = %v, want 1", q.RadiusRatio)
	}
	if !approx(q.Shape, 1) {
		t.Errorf("shape = %v, want 1", q.Shape)
	}
	if !approx(q.Area, math.Sqrt(3)/4) {
		t.Errorf("area = %v", q.Area)
	}
	if !approx(q.EdgeRatio, 1) {
		t.Errorf("edge ratio = %v", q.EdgeRatio)
	}
}

func TestSliverDetection(t *testing.T) {
	a := NewVec2(0, 0)
	b := NewVec2(10, 0)
	c := NewVec2(5, 0.1)
	if !IsSliver(a, b, c, 20) {
		t.Errorf("should be a sliver")
	}
	if IsWellShaped(a, b, c, 20) {
		t.Errorf("should not be well shaped")
	}
}

// ---- Smoothing ----

func TestLaplacianSmooth(t *testing.T) {
	// Grid mesh; move the single interior vertex off-centre then smooth back.
	m := GridMesh(3, 3, 0, 0, 1, 1)
	interior := 4 // center node of 3x3
	m.Vertices[interior] = Vec2{1.4, 1.4}
	sm := LaplacianSmooth(m, 50, 1.0, true)
	if !sm.Vertices[interior].ApproxEqual(Vec2{1, 1}, 1e-6) {
		t.Errorf("smoothed interior = %v, want (1,1)", sm.Vertices[interior])
	}
	// Boundary vertices must be unchanged.
	if !sm.Vertices[0].ApproxEqual(Vec2{0, 0}, tol) {
		t.Errorf("boundary moved: %v", sm.Vertices[0])
	}
}

func TestLloydSmooth(t *testing.T) {
	m := GridMesh(3, 3, 0, 0, 1, 1)
	m.Vertices[4] = Vec2{1.3, 1.2}
	sm := LloydSmooth(m, 100, true)
	if !sm.Vertices[4].ApproxEqual(Vec2{1, 1}, 1e-3) {
		t.Errorf("lloyd interior = %v, want ~(1,1)", sm.Vertices[4])
	}
}

// ---- Connectivity ----

func TestConnectivity(t *testing.T) {
	m := GridMesh(2, 2, 0, 0, 1, 1)
	if m.NumTriangles() != 2 {
		t.Fatalf("triangles = %d", m.NumTriangles())
	}
	if m.NumEdges() != 5 {
		t.Errorf("edges = %d, want 5", m.NumEdges())
	}
	if got := len(m.BoundaryEdges()); got != 4 {
		t.Errorf("boundary edges = %d, want 4", got)
	}
	if got := len(m.InteriorEdges()); got != 1 {
		t.Errorf("interior edges = %d, want 1", got)
	}
	if !m.IsManifold() {
		t.Errorf("should be manifold")
	}
	if m.ConnectedComponents() != 1 {
		t.Errorf("components = %d, want 1", m.ConnectedComponents())
	}
	if m.EulerCharacteristic() != 1 {
		t.Errorf("Euler = %d, want 1", m.EulerCharacteristic())
	}
	if !approx(m.BoundaryLength(), 4) {
		t.Errorf("boundary length = %v, want 4", m.BoundaryLength())
	}
}

// ---- Marching squares ----

func TestMarchingSquares(t *testing.T) {
	g := SampleGrid2(3, 3, 0, 0, 0.5, 0.5, func(x, y float64) float64 { return x })
	segs := MarchingSquares(g, 0.25)
	if len(segs) != 2 {
		t.Fatalf("segments = %d, want 2", len(segs))
	}
	for _, s := range segs {
		if !approx(s.A.X, 0.25) || !approx(s.B.X, 0.25) {
			t.Errorf("segment not at x=0.25: %v-%v", s.A, s.B)
		}
	}
	if got := TotalSegmentLength(segs); !approx(got, 1) {
		t.Errorf("total length = %v, want 1", got)
	}
}

func TestMarchingSquaresCircle(t *testing.T) {
	// Circle radius 1: f = x^2+y^2, iso=1 over [-1.5,1.5]^2.
	g := SampleGrid2(31, 31, -1.5, -1.5, 0.1, 0.1, func(x, y float64) float64 { return x*x + y*y })
	segs := MarchingSquares(g, 1.0)
	if len(segs) == 0 {
		t.Fatalf("no segments")
	}
	// Every endpoint should lie near the unit circle.
	for _, s := range segs {
		if math.Abs(s.A.Norm()-1) > 0.05 {
			t.Errorf("endpoint off circle: %v (r=%v)", s.A, s.A.Norm())
		}
	}
	// The contour length should approximate 2*pi.
	if l := TotalSegmentLength(segs); math.Abs(l-2*math.Pi) > 0.2 {
		t.Errorf("contour length = %v, want ~%v", l, 2*math.Pi)
	}
}

// ---- Marching cubes / tetrahedra ----

func TestMarchingTetrahedron(t *testing.T) {
	p := [4]Vec3{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
	// One corner inside.
	tris := MarchingTetrahedron(p, [4]float64{1, -1, -1, -1}, 0)
	if len(tris) != 1 {
		t.Errorf("one-in case: %d tris, want 1", len(tris))
	}
	// Two-two split gives a quad (two triangles).
	tris = MarchingTetrahedron(p, [4]float64{1, 1, -1, -1}, 0)
	if len(tris) != 2 {
		t.Errorf("two-in case: %d tris, want 2", len(tris))
	}
	// All inside or all outside gives nothing.
	if len(MarchingTetrahedron(p, [4]float64{1, 1, 1, 1}, 0)) != 0 {
		t.Errorf("all-in should give 0 tris")
	}
	if len(MarchingTetrahedron(p, [4]float64{-1, -1, -1, -1}, 0)) != 0 {
		t.Errorf("all-out should give 0 tris")
	}
}

func TestMarchingCubesSphere(t *testing.T) {
	// f = 1 - (x^2+y^2+z^2); inside sphere f>0. iso=0 => unit sphere.
	f := func(x, y, z float64) float64 { return 1 - (x*x + y*y + z*z) }
	g := SampleGrid3(21, 21, 21, -1.5, -1.5, -1.5, 0.15, 0.15, 0.15, f)
	tris := MarchingCubes(g, 0)
	if len(tris) == 0 {
		t.Fatalf("no triangles from sphere field")
	}
	// Every vertex should lie near the unit sphere.
	for _, tr := range tris {
		for _, v := range [3]Vec3{tr.A, tr.B, tr.C} {
			if math.Abs(v.Norm()-1) > 0.15 {
				t.Errorf("vertex off sphere: %v (r=%v)", v, v.Norm())
				break
			}
		}
	}
	// Surface area should approximate 4*pi.
	area := SurfaceArea(tris)
	if math.Abs(area-4*math.Pi) > 1.5 {
		t.Errorf("surface area = %v, want ~%v", area, 4*math.Pi)
	}
}

func TestMarchingCubesUniform(t *testing.T) {
	g := SampleGrid3(5, 5, 5, 0, 0, 0, 1, 1, 1, func(x, y, z float64) float64 { return 5 })
	if got := len(MarchingCubes(g, 0)); got != 0 {
		t.Errorf("uniform field gave %d tris, want 0", got)
	}
}

// ---- Grid / polygon meshing ----

func TestGridMeshAndEarClip(t *testing.T) {
	m := GridMesh(3, 3, 0, 0, 1, 1)
	if m.NumTriangles() != 8 {
		t.Errorf("grid triangles = %d, want 8", m.NumTriangles())
	}
	if !approx(m.TotalArea(), 4) {
		t.Errorf("grid area = %v, want 4", m.TotalArea())
	}

	// L-shaped non-convex polygon.
	poly := []Vec2{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}
	em, err := TriangulatePolygonEarClip(poly)
	if err != nil {
		t.Fatalf("ear clip: %v", err)
	}
	if em.NumTriangles() != len(poly)-2 {
		t.Errorf("ear clip triangles = %d, want %d", em.NumTriangles(), len(poly)-2)
	}
	if !approx(em.TotalArea(), PolygonArea(poly)) {
		t.Errorf("ear clip area = %v, want %v", em.TotalArea(), PolygonArea(poly))
	}
}

// ---- Quadtree ----

func TestQuadtree(t *testing.T) {
	root := NewQuad(Vec2{0, 0}, Vec2{1, 1})
	qt := BuildQuadtree(root, 2, func(cell Quad, depth int) bool { return true })
	if qt.Depth() != 2 {
		t.Errorf("depth = %d, want 2", qt.Depth())
	}
	if qt.LeafCount() != 16 {
		t.Errorf("leaf count = %d, want 16", qt.LeafCount())
	}
	// Total leaf area equals the root area.
	var area float64
	for _, q := range qt.Leaves() {
		area += q.Area()
	}
	if !approx(area, 1) {
		t.Errorf("total leaf area = %v, want 1", area)
	}
	if len(qt.LeafCorners()) != 25 {
		t.Errorf("leaf corners = %d, want 25", len(qt.LeafCorners()))
	}
}

// ---- Stats ----

func TestMeshStats(t *testing.T) {
	m := GridMesh(3, 3, 0, 0, 1, 1)
	s := m.Stats()
	if s.NumTriangles != 8 || s.NumVertices != 9 {
		t.Errorf("stats counts = %+v", s)
	}
	if !approx(s.TotalArea, 4) {
		t.Errorf("total area = %v", s.TotalArea)
	}
	if !approx(s.MinAngleDeg, 45) {
		t.Errorf("min angle = %v, want 45", s.MinAngleDeg)
	}
	if !approx(s.MaxAngleDeg, 90) {
		t.Errorf("max angle = %v, want 90", s.MaxAngleDeg)
	}
	areas := m.TriangleAreas()
	if !approx(MeanFloat(areas), 0.5) {
		t.Errorf("mean area = %v", MeanFloat(areas))
	}
	if !approx(Median(areas), 0.5) {
		t.Errorf("median area = %v", Median(areas))
	}
}

// ---- Seeded sampling ----

func TestRandomPointsDeterministic(t *testing.T) {
	a := RandomPoints(20, 42, Vec2{0, 0}, Vec2{1, 1})
	b := RandomPoints(20, 42, Vec2{0, 0}, Vec2{1, 1})
	if len(a) != 20 {
		t.Fatalf("len = %d", len(a))
	}
	for i := range a {
		if !a[i].Equal(b[i]) {
			t.Errorf("seed not deterministic at %d", i)
		}
		if a[i].X < 0 || a[i].X > 1 || a[i].Y < 0 || a[i].Y > 1 {
			t.Errorf("point out of bounds: %v", a[i])
		}
	}
	// Triangulate the random cloud; its total area should closely match the
	// area of the convex hull that bounds it.
	m, err := DelaunayMesh(a)
	if err != nil {
		t.Fatalf("DelaunayMesh: %v", err)
	}
	hullArea := PolygonArea(ConvexHull(a))
	if rel := math.Abs(m.TotalArea()-hullArea) / hullArea; rel > 0.01 {
		t.Errorf("delaunay area %v vs hull area %v (rel %v)", m.TotalArea(), hullArea, rel)
	}
}

// ---- Example ----

func ExampleDelaunayMesh() {
	pts := []Vec2{{0, 0}, {1, 0}, {1, 1}, {0, 1}}
	m, _ := DelaunayMesh(pts)
	fmt.Printf("triangles=%d area=%.2f\n", m.NumTriangles(), m.TotalArea())
	// Output: triangles=2 area=1.00
}

func ExampleMarchingSquares() {
	g := SampleGrid2(3, 3, 0, 0, 0.5, 0.5, func(x, y float64) float64 { return x })
	segs := MarchingSquares(g, 0.5)
	fmt.Printf("segments=%d length=%.2f\n", len(segs), TotalSegmentLength(segs))
	// Output: segments=2 length=1.00
}
