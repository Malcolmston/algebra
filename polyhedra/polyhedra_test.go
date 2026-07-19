package polyhedra

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b, eps float64) bool { return math.Abs(a-b) <= eps }

func TestVec3Arithmetic(t *testing.T) {
	a := NewVec3(1, 2, 3)
	b := NewVec3(4, -5, 6)
	if !a.Add(b).Equal(Vec3{5, -3, 9}, tol) {
		t.Errorf("Add = %v", a.Add(b))
	}
	if !a.Sub(b).Equal(Vec3{-3, 7, -3}, tol) {
		t.Errorf("Sub = %v", a.Sub(b))
	}
	if !a.Scale(2).Equal(Vec3{2, 4, 6}, tol) {
		t.Errorf("Scale = %v", a.Scale(2))
	}
	if got := a.Dot(b); !approx(got, 4-10+18, tol) {
		t.Errorf("Dot = %v want 12", got)
	}
	if !a.Cross(b).Equal(Vec3{2*6 - 3*(-5), 3*4 - 1*6, 1*(-5) - 2*4}, tol) {
		t.Errorf("Cross = %v", a.Cross(b))
	}
	if got := NewVec3(3, 4, 0).Len(); !approx(got, 5, tol) {
		t.Errorf("Len = %v want 5", got)
	}
	if u, ok := NewVec3(0, 0, 2).Normalize(); !ok || !u.Equal(ZAxis(), tol) {
		t.Errorf("Normalize = %v", u)
	}
	if _, ok := ZeroVec3().Normalize(); ok {
		t.Errorf("Normalize of zero should fail")
	}
	if got := XAxis().Angle(YAxis()); !approx(got, math.Pi/2, tol) {
		t.Errorf("Angle = %v want pi/2", got)
	}
	if got := NewVec3(1, 0, 0).Triple(YAxis(), ZAxis()); !approx(got, 1, tol) {
		t.Errorf("Triple = %v want 1", got)
	}
}

func TestVec3ProjectReject(t *testing.T) {
	v := NewVec3(2, 3, 0)
	onto := XAxis()
	if !v.Project(onto).Equal(Vec3{2, 0, 0}, tol) {
		t.Errorf("Project = %v", v.Project(onto))
	}
	if !v.Reject(onto).Equal(Vec3{0, 3, 0}, tol) {
		t.Errorf("Reject = %v", v.Reject(onto))
	}
	if !NewVec3(1, -1, 0).Reflect(YAxis()).Equal(Vec3{1, 1, 0}, tol) {
		t.Errorf("Reflect wrong")
	}
}

func TestTriangleAndTetraHelpers(t *testing.T) {
	a, b, c := Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}
	if got := TriangleArea(a, b, c); !approx(got, 0.5, tol) {
		t.Errorf("TriangleArea = %v want 0.5", got)
	}
	n, ok := TriangleUnitNormal(a, b, c)
	if !ok || !n.Equal(ZAxis(), tol) {
		t.Errorf("normal = %v", n)
	}
	d := Vec3{0, 0, 3}
	if got := TetraVolume(a, b, c, d); !approx(got, 0.5, tol) {
		t.Errorf("TetraVolume = %v want 0.5", got)
	}
	if got := SignedTetraVolume(a, b, c, d); !approx(got, 0.5, tol) {
		t.Errorf("SignedTetraVolume = %v want 0.5", got)
	}
}

func TestRegularPolygonHelpers(t *testing.T) {
	// Unit square.
	if got := RegularPolygonArea(4, 1); !approx(got, 1, tol) {
		t.Errorf("square area = %v want 1", got)
	}
	if got := RegularPolygonCircumradius(4, math.Sqrt2); !approx(got, 1, tol) {
		t.Errorf("square circumradius = %v want 1", got)
	}
	if got := RegularPolygonInradius(4, 2); !approx(got, 1, tol) {
		t.Errorf("square inradius = %v want 1", got)
	}
	if got := RegularPolygonInteriorAngle(6); !approx(got, 2*math.Pi/3, tol) {
		t.Errorf("hexagon interior angle = %v", got)
	}
	// Equilateral triangle area.
	if got := RegularPolygonArea(3, 2); !approx(got, math.Sqrt(3), tol) {
		t.Errorf("triangle area = %v want sqrt3", got)
	}
}

func TestPlatonicAnalytic(t *testing.T) {
	tests := []struct {
		name           string
		solid          PlatonicSolid
		v, e, f        int
		vol, area, dih float64
		circ, inr, mid float64
	}{
		{"tetra", NewTetrahedron(1), 4, 6, 4,
			0.1178511302, 1.7320508076, math.Acos(1.0 / 3),
			0.6123724357, 0.2041241452, 0.3535533906},
		{"cube", NewCube(1), 8, 12, 6,
			1, 6, math.Pi / 2,
			0.8660254038, 0.5, 0.7071067812},
		{"octa", NewOctahedron(1), 6, 12, 8,
			0.4714045208, 3.4641016151, math.Acos(-1.0 / 3),
			0.7071067812, 0.4082482905, 0.5},
		{"dodeca", NewDodecahedron(1), 20, 30, 12,
			7.6631189606, 20.6457288071, math.Acos(-1 / math.Sqrt(5)),
			1.4012585384, 1.1135163644, 1.3090169944},
		{"icosa", NewIcosahedron(1), 12, 30, 20,
			2.1816949907, 8.6602540378, math.Acos(-math.Sqrt(5) / 3),
			0.9510565163, 0.7557613141, 0.8090169944},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			s := tc.solid
			if s.NumV != tc.v || s.NumE != tc.e || s.NumF != tc.f {
				t.Errorf("VEF = %d,%d,%d want %d,%d,%d", s.NumV, s.NumE, s.NumF, tc.v, tc.e, tc.f)
			}
			if s.EulerCharacteristic() != 2 {
				t.Errorf("euler = %d want 2", s.EulerCharacteristic())
			}
			if !approx(s.Volume(), tc.vol, 1e-8) {
				t.Errorf("Volume = %v want %v", s.Volume(), tc.vol)
			}
			if !approx(s.SurfaceArea(), tc.area, 1e-8) {
				t.Errorf("SurfaceArea = %v want %v", s.SurfaceArea(), tc.area)
			}
			if !approx(s.DihedralAngle(), tc.dih, 1e-9) {
				t.Errorf("Dihedral = %v want %v", s.DihedralAngle(), tc.dih)
			}
			if !approx(s.Circumradius(), tc.circ, 1e-8) {
				t.Errorf("Circumradius = %v want %v", s.Circumradius(), tc.circ)
			}
			if !approx(s.Inradius(), tc.inr, 1e-8) {
				t.Errorf("Inradius = %v want %v", s.Inradius(), tc.inr)
			}
			if !approx(s.Midradius(), tc.mid, 1e-8) {
				t.Errorf("Midradius = %v want %v", s.Midradius(), tc.mid)
			}
		})
	}
}

// TestPlatonicMeshMatchesAnalytic verifies that the materialised meshes
// reproduce the closed-form volume, area, radii and dihedral angle.
func TestPlatonicMeshMatchesAnalytic(t *testing.T) {
	for _, s := range PlatonicSolids(1.3) {
		t.Run(s.Name, func(t *testing.T) {
			m := s.Mesh()
			if err := m.Validate(); err != nil {
				t.Fatalf("validate: %v", err)
			}
			if m.NumVertices() != s.NumV || m.NumEdges() != s.NumE || m.NumFaces() != s.NumF {
				t.Errorf("mesh VEF = %d,%d,%d want %d,%d,%d",
					m.NumVertices(), m.NumEdges(), m.NumFaces(), s.NumV, s.NumE, s.NumF)
			}
			if m.EulerCharacteristic() != 2 {
				t.Errorf("euler = %d", m.EulerCharacteristic())
			}
			if !m.IsClosed() {
				t.Errorf("mesh not closed")
			}
			if !IsConvex(m, 1e-6) {
				t.Errorf("mesh not convex")
			}
			if !approx(m.Volume(), s.Volume(), 1e-6) {
				t.Errorf("mesh volume = %v want %v", m.Volume(), s.Volume())
			}
			if !approx(m.SurfaceArea(), s.SurfaceArea(), 1e-6) {
				t.Errorf("mesh area = %v want %v", m.SurfaceArea(), s.SurfaceArea())
			}
			if !approx(m.Inradius(), s.Inradius(), 1e-6) {
				t.Errorf("mesh inradius = %v want %v", m.Inradius(), s.Inradius())
			}
			if !approx(m.Circumradius(), s.Circumradius(), 1e-6) {
				t.Errorf("mesh circumradius = %v want %v", m.Circumradius(), s.Circumradius())
			}
			if !approx(m.Midradius(), s.Midradius(), 1e-6) {
				t.Errorf("mesh midradius = %v want %v", m.Midradius(), s.Midradius())
			}
			d, ok := DihedralAngleAtEdge(m, m.EdgeList()[0])
			if !ok || !approx(d, s.DihedralAngle(), 1e-6) {
				t.Errorf("mesh dihedral = %v want %v", d, s.DihedralAngle())
			}
			// Centroid at the origin (meshes are centered).
			c, _ := m.Centroid()
			if !c.Equal(ZeroVec3(), 1e-6) {
				t.Errorf("centroid = %v want origin", c)
			}
		})
	}
}

func TestEdgeScaling(t *testing.T) {
	// Volume scales as edge^3, area as edge^2.
	a := NewCube(1)
	b := NewCube(2)
	if !approx(b.Volume(), 8*a.Volume(), tol) {
		t.Errorf("volume scaling wrong")
	}
	if !approx(b.SurfaceArea(), 4*a.SurfaceArea(), tol) {
		t.Errorf("area scaling wrong")
	}
}

func TestArchimedeanConsistency(t *testing.T) {
	for _, s := range ArchimedeanSolids() {
		t.Run(s.Name, func(t *testing.T) {
			if s.EulerCharacteristic() != 2 {
				t.Errorf("euler = %d want 2", s.EulerCharacteristic())
			}
			if s.TotalFaceCount() != s.NumF {
				t.Errorf("face count %d != NumF %d", s.TotalFaceCount(), s.NumF)
			}
			// Edge-count check: sum over faces of sides = 2E.
			sumSides := 0
			for _, ft := range s.Faces {
				sumSides += ft.Sides * ft.Count
			}
			if sumSides != 2*s.NumE {
				t.Errorf("sum of sides %d != 2E %d", sumSides, 2*s.NumE)
			}
			// Midradius consistency: R^2 = mid^2 + (a/2)^2.
			r := s.Circumradius(1)
			mid := s.Midradius(1)
			if !approx(r*r, mid*mid+0.25, 1e-6) {
				t.Errorf("mid/circum inconsistent: R=%v mid=%v", r, mid)
			}
		})
	}
}

func TestArchimedeanKnownValues(t *testing.T) {
	tests := []struct {
		ctor      func() ArchimedeanSolid
		area, vol float64
		circ      float64
	}{
		{NewTruncatedTetrahedron, 7 * math.Sqrt(3), 23 * math.Sqrt2 / 12, math.Sqrt(22) / 4},
		{NewCuboctahedron, 6 + 2*math.Sqrt(3), 5 * math.Sqrt2 / 3, 1},
		{NewTruncatedOctahedron, 6 + 12*math.Sqrt(3), 8 * math.Sqrt2, math.Sqrt(10) / 2},
		{NewIcosidodecahedron, 0, (45 + 17*math.Sqrt(5)) / 6, Phi},
		{NewTruncatedIcosahedron, 0, (125 + 43*math.Sqrt(5)) / 4, 0},
	}
	for _, tc := range tests {
		s := tc.ctor()
		if !approx(s.Volume(1), tc.vol, 1e-6) {
			t.Errorf("%s volume = %v want %v", s.Name, s.Volume(1), tc.vol)
		}
		if tc.area != 0 && !approx(s.SurfaceArea(1), tc.area, 1e-6) {
			t.Errorf("%s area = %v want %v", s.Name, s.SurfaceArea(1), tc.area)
		}
		if tc.circ != 0 && !approx(s.Circumradius(1), tc.circ, 1e-6) {
			t.Errorf("%s circ = %v want %v", s.Name, s.Circumradius(1), tc.circ)
		}
	}
	// Snub solids (numeric references).
	if !approx(NewSnubCube().SurfaceArea(1), 6+8*math.Sqrt(3), 1e-6) {
		t.Errorf("snub cube area wrong")
	}
	if !approx(NewSnubCube().Volume(1), 7.8894774, 1e-5) {
		t.Errorf("snub cube volume wrong")
	}
	if !approx(NewSnubDodecahedron().Volume(1), 37.6166499, 1e-4) {
		t.Errorf("snub dodeca volume wrong")
	}
}

func TestConvexHullCube(t *testing.T) {
	pts := []Vec3{
		{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0},
		{0, 0, 2}, {2, 0, 2}, {2, 2, 2}, {0, 2, 2},
		{1, 1, 1}, // interior point must be dropped
	}
	h, err := ConvexHull(pts)
	if err != nil {
		t.Fatalf("hull: %v", err)
	}
	if h.NumVertices() != 8 {
		t.Errorf("hull vertices = %d want 8", h.NumVertices())
	}
	if !approx(h.Volume(), 8, 1e-6) {
		t.Errorf("hull volume = %v want 8", h.Volume())
	}
	if !approx(h.SurfaceArea(), 24, 1e-6) {
		t.Errorf("hull area = %v want 24", h.SurfaceArea())
	}
	if !h.IsClosed() || !IsConvex(h, 1e-6) {
		t.Errorf("hull not closed/convex")
	}
	// Merged polygons: a cube has six quadrilateral faces.
	hp, err := ConvexHullPolygons(pts)
	if err != nil {
		t.Fatalf("hull polygons: %v", err)
	}
	if hp.NumFaces() != 6 {
		t.Errorf("merged faces = %d want 6", hp.NumFaces())
	}
	if hp.EulerCharacteristic() != 2 {
		t.Errorf("merged euler = %d want 2", hp.EulerCharacteristic())
	}
	if !approx(hp.Volume(), 8, 1e-6) {
		t.Errorf("merged volume = %v want 8", hp.Volume())
	}
}

func TestConvexHullMatchesIcosahedron(t *testing.T) {
	ico := NewIcosahedron(1).Mesh()
	h, err := ConvexHull(ico.Verts)
	if err != nil {
		t.Fatalf("hull: %v", err)
	}
	if h.NumVertices() != 12 {
		t.Errorf("hull vertices = %d want 12", h.NumVertices())
	}
	if !approx(h.Volume(), ico.Volume(), 1e-6) {
		t.Errorf("hull volume = %v want %v", h.Volume(), ico.Volume())
	}
	// GiftWrapHull agrees.
	g, err := GiftWrapHull(ico.Verts)
	if err != nil {
		t.Fatalf("giftwrap: %v", err)
	}
	if !approx(g.Volume(), h.Volume(), 1e-6) {
		t.Errorf("giftwrap volume disagrees: %v vs %v", g.Volume(), h.Volume())
	}
}

func TestConvexHullDegenerate(t *testing.T) {
	// Coplanar points -> no 3-D hull.
	pts := []Vec3{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 3, 0}}
	if _, err := ConvexHull(pts); err == nil {
		t.Errorf("expected degenerate error for coplanar points")
	}
	if _, err := ConvexHull([]Vec3{{0, 0, 0}, {1, 0, 0}}); err == nil {
		t.Errorf("expected error for too few points")
	}
}

func TestDualPolyhedra(t *testing.T) {
	tests := []struct {
		name       string
		solid      PlatonicSolid
		dv, de, df int
	}{
		{"cube->octa", NewCube(2), 6, 12, 8},
		{"octa->cube", NewOctahedron(2), 8, 12, 6},
		{"dodeca->icosa", NewDodecahedron(1), 12, 30, 20},
		{"icosa->dodeca", NewIcosahedron(1), 20, 30, 12},
		{"tetra->tetra", NewTetrahedron(1), 4, 6, 4},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			m := tc.solid.Mesh()
			d, err := m.Dual()
			if err != nil {
				t.Fatalf("dual: %v", err)
			}
			if d.NumVertices() != tc.dv || d.NumEdges() != tc.de || d.NumFaces() != tc.df {
				t.Errorf("dual VEF = %d,%d,%d want %d,%d,%d",
					d.NumVertices(), d.NumEdges(), d.NumFaces(), tc.dv, tc.de, tc.df)
			}
			if d.EulerCharacteristic() != 2 {
				t.Errorf("dual euler = %d", d.EulerCharacteristic())
			}
			if !d.IsClosed() {
				t.Errorf("dual not closed")
			}
			if !IsConvex(d, 1e-9) {
				t.Errorf("dual not convex")
			}
		})
	}
}

func TestMeshDivergenceTheorem(t *testing.T) {
	// A translated, unevenly scaled cube: volume via divergence theorem should
	// match width*depth*height. Build an axis-aligned box mesh by hand.
	box := boxMesh(2, 3, 5)
	if !approx(box.Volume(), 30, 1e-9) {
		t.Errorf("box volume = %v want 30", box.Volume())
	}
	if !approx(box.SurfaceArea(), 2*(2*3+2*5+3*5), 1e-9) {
		t.Errorf("box area = %v", box.SurfaceArea())
	}
	c, _ := box.Centroid()
	if !c.Equal(ZeroVec3(), 1e-9) {
		t.Errorf("box centroid = %v want origin", c)
	}
	// Translation invariance of volume.
	moved := box.Translate(Vec3{10, -4, 7})
	if !approx(moved.Volume(), 30, 1e-9) {
		t.Errorf("translated volume = %v want 30", moved.Volume())
	}
	// Face normals point outward: dot with (faceCentroid - center) > 0.
	for i := range moved.Faces {
		n, _ := moved.FaceNormal(i)
		fc, _ := moved.FaceCentroid(i)
		mc, _ := moved.Centroid()
		if fc.Sub(mc).Dot(n) <= 0 {
			t.Errorf("face %d normal not outward", i)
		}
	}
}

func boxMesh(w, d, h float64) *Polyhedron {
	x, y, z := w/2, d/2, h/2
	verts := []Vec3{
		{-x, -y, -z}, {x, -y, -z}, {x, y, -z}, {-x, y, -z},
		{-x, -y, z}, {x, -y, z}, {x, y, z}, {-x, y, z},
	}
	faces := [][]int{
		{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
		{1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7},
	}
	return NewPolyhedron(verts, faces).OrientOutward()
}

func TestPrismAntiprismPyramid(t *testing.T) {
	// Triangular prism volume = base area * height.
	p := Prism(3, 2, 5)
	if p.EulerCharacteristic() != 2 || !p.IsClosed() {
		t.Errorf("prism topology wrong")
	}
	if !approx(p.Volume(), PrismVolume(3, 2, 5), 1e-9) {
		t.Errorf("prism volume = %v want %v", p.Volume(), PrismVolume(3, 2, 5))
	}
	if !approx(p.SurfaceArea(), PrismSurfaceArea(3, 2, 5), 1e-9) {
		t.Errorf("prism area = %v want %v", p.SurfaceArea(), PrismSurfaceArea(3, 2, 5))
	}
	// Square antiprism: all edges equal length a, and it is closed & convex.
	ap := UniformAntiprism(4, 1)
	if !ap.IsClosed() || ap.EulerCharacteristic() != 2 {
		t.Errorf("antiprism topology wrong")
	}
	if !IsConvex(ap, 1e-9) {
		t.Errorf("antiprism not convex")
	}
	for _, l := range ap.EdgeLengths() {
		if !approx(l, 1, 1e-9) {
			t.Errorf("antiprism edge length = %v want 1", l)
		}
	}
	// Pyramid volume = base area * h / 3.
	py := Pyramid(6, 1, 4)
	if !approx(py.Volume(), PyramidVolume(6, 1, 4), 1e-9) {
		t.Errorf("pyramid volume = %v want %v", py.Volume(), PyramidVolume(6, 1, 4))
	}
	if py.EulerCharacteristic() != 2 {
		t.Errorf("pyramid euler = %d", py.EulerCharacteristic())
	}
	// Bipyramid is closed.
	bp := Bipyramid(5, 1, 2)
	if !bp.IsClosed() || bp.EulerCharacteristic() != 2 {
		t.Errorf("bipyramid topology wrong")
	}
}

func TestFaceAdjacencyAndDegrees(t *testing.T) {
	cube := NewCube(1).Mesh()
	adj := FaceAdjacency(cube)
	// Each face of a cube borders exactly four others.
	for i, a := range adj {
		if len(a) != 4 {
			t.Errorf("face %d has %d neighbors want 4", i, len(a))
		}
	}
	// Every cube vertex has degree 3.
	for i, deg := range cube.VertexDegrees() {
		if deg != 3 {
			t.Errorf("vertex %d degree = %d want 3", i, deg)
		}
	}
}

func TestSphereAndSphericity(t *testing.T) {
	if !approx(SphereVolume(2), 4.0/3*math.Pi*8, tol) {
		t.Errorf("sphere volume wrong")
	}
	if !approx(SphereSurfaceArea(3), 4*math.Pi*9, tol) {
		t.Errorf("sphere area wrong")
	}
	// A cube has sphericity pi^(1/3)/6^(1/3)*... < 1; specifically ~0.806.
	s := Sphericity(NewCube(1).Mesh())
	if s <= 0.8 || s >= 0.81 {
		t.Errorf("cube sphericity = %v want ~0.806", s)
	}
	// The icosahedron is rounder than the tetrahedron.
	if Sphericity(NewIcosahedron(1).Mesh()) <= Sphericity(NewTetrahedron(1).Mesh()) {
		t.Errorf("icosahedron should be rounder than tetrahedron")
	}
}

func TestDihedralFromSchlafli(t *testing.T) {
	cases := []struct {
		p, q int
		want float64
	}{
		{3, 3, math.Acos(1.0 / 3)},
		{4, 3, math.Pi / 2},
		{3, 4, math.Acos(-1.0 / 3)},
		{5, 3, math.Acos(-1 / math.Sqrt(5))},
		{3, 5, math.Acos(-math.Sqrt(5) / 3)},
	}
	for _, c := range cases {
		if got := DihedralFromSchlafli(c.p, c.q); !approx(got, c.want, 1e-9) {
			t.Errorf("{%d,%d} dihedral = %v want %v", c.p, c.q, got, c.want)
		}
	}
}

// ExamplePlatonicSolid demonstrates reading properties of a regular icosahedron
// and confirming its Euler characteristic from the materialised mesh.
func ExamplePlatonicSolid() {
	ico := NewIcosahedron(1)
	mesh := ico.Mesh()
	fmt.Printf("faces=%d vertices=%d edges=%d\n", ico.NumF, ico.NumV, ico.NumE)
	fmt.Printf("V-E+F=%d\n", mesh.EulerCharacteristic())
	fmt.Printf("volume=%.4f area=%.4f\n", mesh.Volume(), mesh.SurfaceArea())
	fmt.Printf("dihedral=%.2f deg\n", ico.DihedralAngle()*180/math.Pi)
	// Output:
	// faces=20 vertices=12 edges=30
	// V-E+F=2
	// volume=2.1817 area=8.6603
	// dihedral=138.19 deg
}

// ExampleConvexHull computes the convex hull of a point cloud whose interior
// point is discarded, recovering a unit cube of volume 1.
func ExampleConvexHull() {
	pts := []Vec3{
		{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
		{0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
		{0.5, 0.5, 0.5},
	}
	h, _ := ConvexHull(pts)
	fmt.Printf("vertices=%d volume=%.1f area=%.1f\n", h.NumVertices(), h.Volume(), h.SurfaceArea())
	// Output:
	// vertices=8 volume=1.0 area=6.0
}
