package hull3d

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool { return math.Abs(a-b) <= eps }

func TestVec3Algebra(t *testing.T) {
	a := NewVec3(1, 2, 3)
	b := NewVec3(4, -5, 6)
	tests := []struct {
		name string
		got  float64
		want float64
	}{
		{"dot", a.Dot(b), 4 - 10 + 18},
		{"lengthSq", a.LengthSq(), 14},
		{"length", a.Length(), math.Sqrt(14)},
		{"l1", b.L1(), 15},
		{"linf", b.LInf(), 6},
		{"distanceSq", a.DistanceSq(b), 9 + 49 + 9},
		{"triple", Triple(UnitX(), UnitY(), UnitZ()), 1},
	}
	for _, tc := range tests {
		if !approx(tc.got, tc.want, tol) {
			t.Errorf("%s = %v, want %v", tc.name, tc.got, tc.want)
		}
	}
	if c := a.Cross(b); !c.ApproxEqual(Vec3{2*6 - 3*-5, 3*4 - 1*6, 1*-5 - 2*4}, tol) {
		t.Errorf("cross = %v", c)
	}
	if u, err := a.Normalize(); err != nil || !approx(u.Length(), 1, tol) {
		t.Errorf("normalize failed: %v %v", u, err)
	}
	if _, err := Zero().Normalize(); err == nil {
		t.Error("expected error normalizing zero vector")
	}
	if p, err := a.ProjectOnto(UnitX()); err != nil || !p.ApproxEqual(Vec3{1, 0, 0}, tol) {
		t.Errorf("project = %v %v", p, err)
	}
	ang, _ := UnitX().AngleBetween(UnitY())
	if !approx(ang, math.Pi/2, tol) {
		t.Errorf("angle = %v", ang)
	}
	u, v, err := UnitZ().OrthonormalBasis()
	if err != nil || !approx(u.Dot(v), 0, tol) || !approx(u.Length(), 1, tol) {
		t.Errorf("orthobasis bad: %v %v %v", u, v, err)
	}
}

func TestOrientAndInSphere(t *testing.T) {
	a := Vec3{0, 0, 0}
	b := Vec3{1, 0, 0}
	c := Vec3{0, 1, 0}
	tests := []struct {
		d    Vec3
		sign int
	}{
		{Vec3{0, 0, 1}, -1}, // Orient3D negative for this configuration
		{Vec3{0, 0, -1}, 1},
		{Vec3{1, 1, 0}, 0}, // coplanar
	}
	for _, tc := range tests {
		if s := OrientSign(a, b, c, tc.d, tol); s != tc.sign {
			t.Errorf("OrientSign(%v) = %d, want %d", tc.d, s, tc.sign)
		}
	}
	// Circumsphere of the standard tetra: point at centroid is inside.
	d := Vec3{0, 0, 1}
	if !inCircumsphere(a, b, c, d, Vec3{0.2, 0.2, 0.2}) {
		t.Error("expected point inside circumsphere")
	}
	if inCircumsphere(a, b, c, d, Vec3{9, 9, 9}) {
		t.Error("expected far point outside circumsphere")
	}
}

func TestConvexHullCube(t *testing.T) {
	pts := []Vec3{
		{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
		{1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
		{0.5, 0.5, 0.5}, {0.3, 0.7, 0.2}, // interior points
	}
	h, err := ConvexHull(pts)
	if err != nil {
		t.Fatal(err)
	}
	cases := []struct {
		name string
		got  float64
		want float64
	}{
		{"volume", h.Volume(), 1},
		{"surface", h.SurfaceArea(), 6},
		{"nverts", float64(h.NumVertices()), 8},
		{"euler", float64(h.EulerCharacteristic()), 2},
		{"edges", float64(h.NumEdges()), 18}, // 12 cube edges + 6 face diagonals (triangulated)
	}
	for _, c := range cases {
		if !approx(c.got, c.want, 1e-9) {
			t.Errorf("cube %s = %v, want %v", c.name, c.got, c.want)
		}
	}
	if !h.IsClosed() {
		t.Error("hull not closed")
	}
	if cen := h.Centroid(); !cen.ApproxEqual(Vec3{0.5, 0.5, 0.5}, 1e-9) {
		t.Errorf("centroid = %v", cen)
	}
	if !h.Contains(Vec3{0.5, 0.5, 0.5}, 1e-9) {
		t.Error("should contain center")
	}
	if h.Contains(Vec3{2, 2, 2}, 1e-9) {
		t.Error("should not contain outside point")
	}
	// Cross-check with gift wrap.
	g, err := GiftWrapHull(pts)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(g.Volume(), 1, 1e-9) {
		t.Errorf("giftwrap volume = %v", g.Volume())
	}
	if g.NumVertices() != 8 {
		t.Errorf("giftwrap verts = %d", g.NumVertices())
	}
}

func TestHullDegenerate(t *testing.T) {
	// Coplanar points -> error.
	_, err := ConvexHull([]Vec3{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}})
	if err == nil {
		t.Error("expected error for coplanar points")
	}
	_, err = ConvexHull([]Vec3{{0, 0, 0}, {1, 0, 0}})
	if err == nil {
		t.Error("expected error for too few points")
	}
}

func TestHullRandomVsGiftWrap(t *testing.T) {
	pts := SamplePointsInBall(40, Zero(), 1, 12345)
	h, err := ConvexHull(pts)
	if err != nil {
		t.Fatal(err)
	}
	g, err := GiftWrapHull(pts)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(h.Volume(), g.Volume(), 1e-6) {
		t.Errorf("hull volume %v != giftwrap %v", h.Volume(), g.Volume())
	}
	if h.EulerCharacteristic() != 2 {
		t.Errorf("euler = %d", h.EulerCharacteristic())
	}
	// Every input point should be contained.
	for _, p := range pts {
		if !h.Contains(p, 1e-6) {
			t.Errorf("point %v not contained", p)
		}
	}
}

func TestPrimitives(t *testing.T) {
	oct := Octahedron(1)
	if !approx(oct.Volume(), 4.0/3.0, 1e-9) {
		t.Errorf("octahedron volume = %v, want %v", oct.Volume(), 4.0/3.0)
	}
	if !approx(oct.SurfaceArea(), 4*math.Sqrt(3), 1e-9) {
		t.Errorf("octahedron area = %v", oct.SurfaceArea())
	}
	if oct.NumVertices() != 6 || oct.NumFaces() != 8 {
		t.Errorf("octahedron V=%d F=%d", oct.NumVertices(), oct.NumFaces())
	}

	box := BoxPolytope(Zero(), Vec3{1, 2, 3})
	if !approx(box.Volume(), 2*4*6, 1e-9) {
		t.Errorf("box volume = %v", box.Volume())
	}

	ico := Icosahedron(1)
	if ico.NumVertices() != 12 || ico.NumFaces() != 20 {
		t.Errorf("icosahedron V=%d F=%d", ico.NumVertices(), ico.NumFaces())
	}
	for _, v := range ico.Vertices {
		if !approx(v.Length(), 1, 1e-9) {
			t.Errorf("icosahedron vertex off sphere: %v", v.Length())
		}
	}

	dod := Dodecahedron(1)
	if dod.NumVertices() != 20 {
		t.Errorf("dodecahedron V=%d", dod.NumVertices())
	}
	if dod.EulerCharacteristic() != 2 {
		t.Errorf("dodecahedron euler = %d", dod.EulerCharacteristic())
	}

	tet := RegularTetrahedron(1)
	for _, v := range tet.Vertices {
		if !approx(v.Length(), 1, 1e-9) {
			t.Errorf("tetra vertex off circumsphere: %v", v.Length())
		}
	}
}

func TestGJKDistance(t *testing.T) {
	tests := []struct {
		name string
		a, b ConvexShape
		want float64
	}{
		{"sphere-sphere", SphereShape{Vec3{0, 0, 0}, 1}, SphereShape{Vec3{5, 0, 0}, 1}, 3},
		{"box-box", BoxShape{Vec3{0, 0, 0}, Splat(1)}, BoxShape{Vec3{4, 0, 0}, Splat(1)}, 2},
		{"point-box", PointShape{Vec3{3, 0, 0}}, BoxShape{Zero(), Splat(1)}, 2},
		{"seg-seg", SegmentShape{Vec3{0, 0, 0}, Vec3{0, 0, 1}}, SegmentShape{Vec3{2, 0, 0}, Vec3{2, 1, 0}}, 2},
	}
	for _, tc := range tests {
		d, _, _, err := GJKDistance(tc.a, tc.b)
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		if !approx(d, tc.want, 1e-6) {
			t.Errorf("%s distance = %v, want %v", tc.name, d, tc.want)
		}
		if GJKIntersect(tc.a, tc.b) {
			t.Errorf("%s should not intersect", tc.name)
		}
	}
}

func TestGJKWitnessPoints(t *testing.T) {
	a := BoxShape{Vec3{0, 0, 0}, Splat(1)}
	b := BoxShape{Vec3{4, 0, 0}, Splat(1)}
	d, pa, pb, err := GJKDistance(a, b)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(pa.X, 1, 1e-6) || !approx(pb.X, 3, 1e-6) {
		t.Errorf("witness points pa=%v pb=%v", pa, pb)
	}
	if !approx(pa.Distance(pb), d, 1e-6) {
		t.Errorf("witness distance mismatch")
	}
}

func TestEPAPenetration(t *testing.T) {
	tests := []struct {
		name string
		a, b ConvexShape
		want float64
	}{
		{"spheres", SphereShape{Vec3{0, 0, 0}, 1}, SphereShape{Vec3{1.5, 0, 0}, 1}, 0.5},
		{"boxes", BoxShape{Vec3{0, 0, 0}, Splat(1)}, BoxShape{Vec3{1.5, 0, 0}, Splat(1)}, 0.5},
		{"boxes-y", BoxShape{Vec3{0, 0, 0}, Splat(1)}, BoxShape{Vec3{0, 1.75, 0}, Splat(1)}, 0.25},
	}
	for _, tc := range tests {
		if !GJKIntersect(tc.a, tc.b) {
			t.Fatalf("%s: expected intersection", tc.name)
		}
		pen, err := EPAPenetration(tc.a, tc.b)
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		if !approx(pen.Depth, tc.want, 1e-2) {
			t.Errorf("%s depth = %v, want %v", tc.name, pen.Depth, tc.want)
		}
		if !approx(pen.Normal.Length(), 1, 1e-6) {
			t.Errorf("%s normal not unit: %v", tc.name, pen.Normal.Length())
		}
	}
	// Non-overlapping shapes return an error.
	if _, err := EPAPenetration(SphereShape{Zero(), 1}, SphereShape{Vec3{5, 0, 0}, 1}); err == nil {
		t.Error("expected error for non-overlapping EPA")
	}
}

func TestMinkowskiSum(t *testing.T) {
	a := BoxPolytope(Zero(), Splat(0.5)) // unit cube
	b := BoxPolytope(Zero(), Splat(0.5))
	sum, err := MinkowskiSum(a, b)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(sum.Volume(), 8, 1e-9) { // cube of side 2
		t.Errorf("minkowski sum volume = %v, want 8", sum.Volume())
	}
	// Overlap test via Minkowski difference containing origin.
	c1 := BoxPolytope(Zero(), Splat(1))
	c2 := BoxPolytope(Vec3{1.5, 0, 0}, Splat(1))
	over, err := PolytopesOverlap(c1, c2, 1e-9)
	if err != nil {
		t.Fatal(err)
	}
	if !over {
		t.Error("boxes should overlap")
	}
	c3 := BoxPolytope(Vec3{5, 0, 0}, Splat(1))
	over, _ = PolytopesOverlap(c1, c3, 1e-9)
	if over {
		t.Error("distant boxes should not overlap")
	}
}

func TestHalfSpaceIntersection(t *testing.T) {
	cube := BoxPolytope(Zero(), Splat(1))
	hs := cube.HalfSpaces()
	rebuilt, err := HalfSpaceIntersection(hs, Zero())
	if err != nil {
		t.Fatal(err)
	}
	if !approx(rebuilt.Volume(), 8, 1e-6) {
		t.Errorf("reconstructed cube volume = %v, want 8", rebuilt.Volume())
	}
	// Chebyshev centre of a unit cube's half-spaces is the origin with radius 1.
	center, radius, err := ChebyshevCenterEstimate(hs)
	if err != nil {
		t.Fatal(err)
	}
	if !center.ApproxEqual(Zero(), 1e-6) || !approx(radius, 1, 1e-6) {
		t.Errorf("chebyshev center=%v radius=%v", center, radius)
	}
}

func TestDelaunay(t *testing.T) {
	pts := SamplePointsInBall(30, Zero(), 1, 999)
	del, err := Delaunay(pts)
	if err != nil {
		t.Fatal(err)
	}
	hullVol, err := HullVolume(pts)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(del.TotalVolume(), hullVol, 1e-6) {
		t.Errorf("delaunay volume %v != hull volume %v", del.TotalVolume(), hullVol)
	}
	if !del.IsDelaunay(1e-9 * hullVol) {
		t.Error("triangulation is not Delaunay")
	}
	if del.NumTets() == 0 {
		t.Error("no tetrahedra produced")
	}
}

func TestCircumcenter(t *testing.T) {
	// Circumcenter of a tetra with a right-angle corner at origin.
	c, r, err := Circumcenter(Vec3{0, 0, 0}, Vec3{2, 0, 0}, Vec3{0, 2, 0}, Vec3{0, 0, 2})
	if err != nil {
		t.Fatal(err)
	}
	if !c.ApproxEqual(Vec3{1, 1, 1}, 1e-9) {
		t.Errorf("circumcenter = %v, want (1,1,1)", c)
	}
	if !approx(r, math.Sqrt(3), 1e-9) {
		t.Errorf("circumradius = %v", r)
	}
}

func TestVoronoi(t *testing.T) {
	pts := []Vec3{
		{0, 0, 0},
		{2, 0, 0}, {-2, 0, 0}, {0, 2, 0}, {0, -2, 0}, {0, 0, 2}, {0, 0, -2},
		{2, 2, 2}, {-2, -2, -2}, {2, -2, 2}, {-2, 2, -2},
	}
	vor, err := VoronoiCells(pts)
	if err != nil {
		t.Fatal(err)
	}
	// The central site's nearest query near origin returns index 0.
	idx, err := vor.NearestSite(Vec3{0.1, 0.1, 0.1})
	if err != nil {
		t.Fatal(err)
	}
	if idx != 0 {
		t.Errorf("nearest site = %d, want 0", idx)
	}
	// The central cell should be bounded.
	cell, _ := vor.Cell(0)
	if !cell.Bounded {
		t.Error("central Voronoi cell should be bounded")
	}
	if !cell.Contains(Vec3{0.2, 0, 0}, 1e-9) {
		t.Error("origin-adjacent point should be in central cell")
	}
}

func TestQueries(t *testing.T) {
	cube := BoxPolytope(Zero(), Splat(1))
	// Ray from outside through center.
	tE, tX, ok := RayPolytopeIntersection(cube, Vec3{-3, 0, 0}, Vec3{1, 0, 0})
	if !ok || !approx(tE, 2, 1e-9) || !approx(tX, 4, 1e-9) {
		t.Errorf("ray hit tE=%v tX=%v ok=%v", tE, tX, ok)
	}
	// Ray missing the cube.
	if _, _, ok := RayPolytopeIntersection(cube, Vec3{-3, 5, 0}, Vec3{1, 0, 0}); ok {
		t.Error("ray should miss")
	}
	// Closest point to an exterior point.
	pt, d, err := ClosestPointInPolytope(cube, Vec3{3, 0, 0})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(d, 2, 1e-6) || !approx(pt.X, 1, 1e-6) {
		t.Errorf("closest point %v dist %v", pt, d)
	}
	// Separating plane between two disjoint cubes.
	a := BoxPolytope(Zero(), Splat(1))
	b := BoxPolytope(Vec3{4, 0, 0}, Splat(1))
	pl, err := SeparatingPlane(a, b, 1e-6)
	if err != nil {
		t.Fatal(err)
	}
	if pl.Eval(Zero()) > 0 || pl.Eval(Vec3{4, 0, 0}) < 0 {
		t.Errorf("separating plane misclassifies: %v", pl)
	}
	// Barycentric coordinates.
	w, err := BarycentricTetrahedron(Vec3{0.25, 0.25, 0.25},
		Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(w[0]+w[1]+w[2]+w[3], 1, 1e-9) {
		t.Errorf("barycentric sum = %v", w[0]+w[1]+w[2]+w[3])
	}
	if !approx(w[0], 0.25, 1e-9) {
		t.Errorf("barycentric w0 = %v", w[0])
	}
}

func TestSegmentClosest(t *testing.T) {
	cp, cq, d := LineSegmentClosest(
		Vec3{0, 0, 0}, Vec3{0, 0, 2},
		Vec3{1, 0, 0}, Vec3{1, 2, 0})
	if !approx(d, 1, 1e-9) {
		t.Errorf("segment distance = %v, want 1", d)
	}
	if !cp.ApproxEqual(Vec3{0, 0, 0}, 1e-9) || !cq.ApproxEqual(Vec3{1, 0, 0}, 1e-9) {
		t.Errorf("closest points cp=%v cq=%v", cp, cq)
	}
}

func TestPlaneOps(t *testing.T) {
	pl, err := PlaneFromPoints(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(math.Abs(pl.SignedDistance(Vec3{0, 0, 3})), 3, 1e-9) {
		t.Errorf("signed distance = %v", pl.SignedDistance(Vec3{0, 0, 3}))
	}
	proj, err := pl.Project(Vec3{2, 2, 5})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(proj.Z, 0, 1e-9) {
		t.Errorf("projection z = %v", proj.Z)
	}
}

func ExampleConvexHull() {
	// Convex hull of the eight corners of a unit cube (plus an interior point).
	pts := []Vec3{
		{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
		{1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
		{0.5, 0.5, 0.5},
	}
	h, err := ConvexHull(pts)
	if err != nil {
		panic(err)
	}
	fmt.Printf("vertices=%d faces=%d\n", h.NumVertices(), h.NumFaces())
	fmt.Printf("volume=%.3f surface=%.3f\n", h.Volume(), h.SurfaceArea())
	// Output:
	// vertices=8 faces=12
	// volume=1.000 surface=6.000
}

func ExampleGJKDistance() {
	a := SphereShape{Center: Vec3{0, 0, 0}, Radius: 1}
	b := SphereShape{Center: Vec3{5, 0, 0}, Radius: 1}
	d, _, _, _ := GJKDistance(a, b)
	fmt.Printf("gap=%.3f\n", d)
	// Output:
	// gap=3.000
}
