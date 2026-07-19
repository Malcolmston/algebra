package tilings

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
)

const tol = 1e-9

func approx(a, b float64) bool { return math.Abs(a-b) <= 1e-6 }

// ------------------------------------------------------------------
// Geometry.
// ------------------------------------------------------------------

func TestPointOps(t *testing.T) {
	p := NewPoint(3, 4)
	if got := p.Norm(); !approx(got, 5) {
		t.Errorf("Norm = %v, want 5", got)
	}
	if got := p.Norm2(); !approx(got, 25) {
		t.Errorf("Norm2 = %v, want 25", got)
	}
	q := NewPoint(1, 2)
	if got := p.Dot(q); !approx(got, 11) {
		t.Errorf("Dot = %v, want 11", got)
	}
	if got := p.Cross(q); !approx(got, 2) {
		t.Errorf("Cross = %v, want 2", got)
	}
	if got := p.Add(q); !got.ApproxEqual(NewPoint(4, 6), tol) {
		t.Errorf("Add = %v", got)
	}
	if got := p.Distance(q); !approx(got, math.Hypot(2, 2)) {
		t.Errorf("Distance = %v", got)
	}
	if got := NewPoint(1, 0).Rotate(math.Pi / 2); !got.ApproxEqual(NewPoint(0, 1), 1e-9) {
		t.Errorf("Rotate = %v, want (0,1)", got)
	}
}

func TestPolygonAreaAndCentroid(t *testing.T) {
	sq := []Point{{0, 0}, {2, 0}, {2, 2}, {0, 2}}
	if got := PolygonArea(sq); !approx(got, 4) {
		t.Errorf("PolygonArea = %v, want 4", got)
	}
	if got := PolygonCentroid(sq); !got.ApproxEqual(NewPoint(1, 1), 1e-9) {
		t.Errorf("PolygonCentroid = %v, want (1,1)", got)
	}
	if !PolygonContains(sq, NewPoint(1, 1)) {
		t.Error("PolygonContains should contain center")
	}
	if PolygonContains(sq, NewPoint(3, 3)) {
		t.Error("PolygonContains should exclude outside point")
	}
}

func TestAngleHelpers(t *testing.T) {
	if got := AngleBetween(NewPoint(1, 0), NewPoint(0, 1)); !approx(got, math.Pi/2) {
		t.Errorf("AngleBetween = %v", got)
	}
	if got := SignedAngle(NewPoint(1, 0), NewPoint(0, -1)); !approx(got, -math.Pi/2) {
		t.Errorf("SignedAngle = %v", got)
	}
	if got := NormalizeAngle(-math.Pi / 2); !approx(got, 3*math.Pi/2) {
		t.Errorf("NormalizeAngle = %v", got)
	}
}

// ------------------------------------------------------------------
// Isometries.
// ------------------------------------------------------------------

func TestIsometryClassify(t *testing.T) {
	tests := []struct {
		name string
		iso  Isometry
		kind IsometryKind
	}{
		{"identity", Identity(), KindIdentity},
		{"translation", Translation(2, 3), KindTranslation},
		{"rotation", RotationAbout(NewPoint(1, 1), math.Pi/2), KindRotation},
		{"reflection", ReflectionLine(NewPoint(0, 1), 0), KindReflection},
		{"glide", GlideReflection(Origin(), 0, 0.5), KindGlideReflection},
	}
	for _, tc := range tests {
		if got := tc.iso.Classify(1e-9).Kind; got != tc.kind {
			t.Errorf("%s: kind = %v, want %v", tc.name, got, tc.kind)
		}
	}
}

func TestRotationClassifyData(t *testing.T) {
	r := RotationAbout(NewPoint(1, 1), math.Pi/2)
	cl := r.Classify(1e-9)
	if !approx(cl.Angle, math.Pi/2) {
		t.Errorf("angle = %v, want pi/2", cl.Angle)
	}
	if !cl.Center.ApproxEqual(NewPoint(1, 1), 1e-9) {
		t.Errorf("center = %v, want (1,1)", cl.Center)
	}
	fp, ok := r.FixedPoint()
	if !ok || !fp.ApproxEqual(NewPoint(1, 1), 1e-9) {
		t.Errorf("FixedPoint = %v ok=%v", fp, ok)
	}
}

func TestGlideData(t *testing.T) {
	g := GlideReflection(Origin(), 0, 0.5)
	if got := g.Apply(NewPoint(0, 1)); !got.ApproxEqual(NewPoint(0.5, -1), 1e-9) {
		t.Errorf("glide apply = %v, want (0.5,-1)", got)
	}
	cl := g.Classify(1e-9)
	if !approx(cl.Glide, 0.5) || !approx(cl.AxisAngle, 0) {
		t.Errorf("glide classify = %+v", cl)
	}
}

func TestComposeInverse(t *testing.T) {
	r := RotationAbout(NewPoint(2, -1), 0.7)
	inv := r.Inverse()
	if got := r.Compose(inv); !got.IsIdentity(1e-9) {
		t.Errorf("r∘r^-1 = %v, want identity", got)
	}
	// Composition applies inner first.
	a := Translation(1, 0)
	b := Rotation(math.Pi / 2)
	got := b.Compose(a).Apply(Origin()) // translate then rotate
	if !got.ApproxEqual(NewPoint(0, 1), 1e-9) {
		t.Errorf("compose = %v, want (0,1)", got)
	}
	if r.Determinant() <= 0 != false {
		t.Error("rotation should be direct")
	}
}

func TestIsometryOrder(t *testing.T) {
	if got := Rotation(2*math.Pi/6).Order(12, 1e-9); got != 6 {
		t.Errorf("order = %v, want 6", got)
	}
	if got := Reflection(0.3).Order(12, 1e-9); got != 2 {
		t.Errorf("reflection order = %v, want 2", got)
	}
	if got := Translation(1, 0).Order(12, 1e-9); got != 0 {
		t.Errorf("translation order = %v, want 0", got)
	}
}

// ------------------------------------------------------------------
// Point groups.
// ------------------------------------------------------------------

func TestPointGroups(t *testing.T) {
	tests := []struct {
		n           int
		cyclicOrder int
		dihedOrder  int
	}{
		{1, 1, 2}, {2, 2, 4}, {3, 3, 6}, {4, 4, 8}, {6, 6, 12},
	}
	for _, tc := range tests {
		c := CyclicGroup(tc.n)
		if c.Order() != tc.cyclicOrder {
			t.Errorf("C%d order = %d, want %d", tc.n, c.Order(), tc.cyclicOrder)
		}
		if !c.IsClosed(1e-9) {
			t.Errorf("C%d not closed", tc.n)
		}
		if c.MaxRotationOrder() != tc.n {
			t.Errorf("C%d maxrot = %d", tc.n, c.MaxRotationOrder())
		}
		d := DihedralGroup(tc.n)
		if d.Order() != tc.dihedOrder {
			t.Errorf("D%d order = %d, want %d", tc.n, d.Order(), tc.dihedOrder)
		}
		if !d.IsClosed(1e-9) {
			t.Errorf("D%d not closed", tc.n)
		}
		if d.NumReflections() != tc.n {
			t.Errorf("D%d reflections = %d, want %d", tc.n, d.NumReflections(), tc.n)
		}
	}
}

func TestPointGroupOrbit(t *testing.T) {
	d := DihedralGroup(4)
	orb := d.Orbit(NewPoint(1, 0), 1e-9)
	if len(orb) != 4 {
		t.Errorf("orbit of (1,0) under D4 has %d points, want 4", len(orb))
	}
	if len(CrystallographicPointGroups()) != 10 {
		t.Error("want 10 crystallographic point groups")
	}
	for _, n := range []int{1, 2, 3, 4, 6} {
		if !IsCrystallographicRestriction(n) {
			t.Errorf("%d should satisfy crystallographic restriction", n)
		}
	}
	for _, n := range []int{5, 7, 8} {
		if IsCrystallographicRestriction(n) {
			t.Errorf("%d should violate crystallographic restriction", n)
		}
	}
}

// ------------------------------------------------------------------
// Lattices.
// ------------------------------------------------------------------

func TestLattices(t *testing.T) {
	if got := SquareLattice(1).FundamentalArea(); !approx(got, 1) {
		t.Errorf("square area = %v", got)
	}
	if got := HexagonalLattice(1).FundamentalArea(); !approx(got, math.Sqrt(3)/2) {
		t.Errorf("hex area = %v", got)
	}
	tests := []struct {
		name string
		lat  Lattice2D
		want LatticeType
	}{
		{"square", SquareLattice(1), Square},
		{"rect", RectangularLattice(1, 1.4), Rectangular},
		{"hex", HexagonalLattice(1), Hexagonal},
		{"rhombic", RhombicLattice(1, 0.8), CenteredRectangular},
		{"oblique", ObliqueLattice(1, 1.35, Deg2Rad(72)), Oblique},
	}
	for _, tc := range tests {
		if got := tc.lat.Classify(1e-6); got != tc.want {
			t.Errorf("%s classify = %v, want %v", tc.name, got, tc.want)
		}
	}
	// Dual of the square lattice is itself.
	d := SquareLattice(2).Dual()
	if !approx(d.V1.Norm(), 0.5) {
		t.Errorf("dual square |V1| = %v, want 0.5", d.V1.Norm())
	}
	// Shortest vector of the hexagonal lattice has length 1.
	if got := HexagonalLattice(1).ShortestLength(); !approx(got, 1) {
		t.Errorf("hex shortest = %v, want 1", got)
	}
	// Reduce lands inside the fundamental cell.
	l := SquareLattice(1)
	r := l.Reduce(NewPoint(2.3, -1.7))
	if !r.ApproxEqual(NewPoint(0.3, 0.3), 1e-9) {
		t.Errorf("reduce = %v, want (0.3,0.3)", r)
	}
}

func TestLatticeHolohedry(t *testing.T) {
	tests := []struct {
		typ   LatticeType
		order int
	}{
		{Oblique, 2}, {Rectangular, 4}, {CenteredRectangular, 4}, {Square, 8}, {Hexagonal, 12},
	}
	for _, tc := range tests {
		if got := tc.typ.PointGroupOrder(); got != tc.order {
			t.Errorf("%v holohedry = %d, want %d", tc.typ, got, tc.order)
		}
	}
}

// ------------------------------------------------------------------
// Wallpaper groups.
// ------------------------------------------------------------------

func TestWallpaperGroups(t *testing.T) {
	if len(WallpaperGroups()) != 17 {
		t.Fatalf("want 17 wallpaper groups, got %d", len(WallpaperGroups()))
	}
	tests := []struct {
		iuc      string
		pgOrder  int
		maxRot   int
		refl     bool
		glide    bool
		orbifold string
	}{
		{"p1", 1, 1, false, false, "o"},
		{"p2", 2, 2, false, false, "2222"},
		{"pm", 2, 1, true, false, "**"},
		{"pg", 2, 1, false, true, "xx"},
		{"cm", 2, 1, true, false, "*x"},
		{"pmm", 4, 2, true, false, "*2222"},
		{"pmg", 4, 2, true, true, "22*"},
		{"pgg", 4, 2, false, true, "22x"},
		{"cmm", 4, 2, true, false, "2*22"},
		{"p4", 4, 4, false, false, "442"},
		{"p4m", 8, 4, true, false, "*442"},
		{"p4g", 8, 4, true, true, "4*2"},
		{"p3", 3, 3, false, false, "333"},
		{"p3m1", 6, 3, true, false, "*333"},
		{"p31m", 6, 3, true, false, "3*3"},
		{"p6", 6, 6, false, false, "632"},
		{"p6m", 12, 6, true, false, "*632"},
	}
	for _, tc := range tests {
		g, ok := WallpaperGroupByIUC(tc.iuc)
		if !ok {
			t.Errorf("%s not found", tc.iuc)
			continue
		}
		if g.PointGroupOrder() != tc.pgOrder {
			t.Errorf("%s point group order = %d, want %d", tc.iuc, g.PointGroupOrder(), tc.pgOrder)
		}
		if g.MaxRotationOrder() != tc.maxRot {
			t.Errorf("%s max rotation = %d, want %d", tc.iuc, g.MaxRotationOrder(), tc.maxRot)
		}
		if g.HasReflection() != tc.refl {
			t.Errorf("%s reflection = %v, want %v", tc.iuc, g.HasReflection(), tc.refl)
		}
		if g.HasGlideReflection() != tc.glide {
			t.Errorf("%s glide = %v, want %v", tc.iuc, g.HasGlideReflection(), tc.glide)
		}
		if g.Orbifold != tc.orbifold {
			t.Errorf("%s orbifold = %s, want %s", tc.iuc, g.Orbifold, tc.orbifold)
		}
	}
}

func TestWallpaperOrbitLatticeInvariance(t *testing.T) {
	g, _ := WallpaperGroupByIUC("p4")
	orb := g.Orbit(NewPoint(0.2, 0.1), 1, 1e-9)
	if len(orb) == 0 {
		t.Fatal("empty orbit")
	}
	// Every group element must map the lattice onto itself.
	lat := g.TranslationLattice()
	for _, e := range g.Cosets() {
		img := e.ApplyVec(lat.V1)
		// image of a basis vector must be an integer combination.
		a, b := lat.ToFractional(img)
		if math.Abs(a-math.Round(a)) > 1e-9 || math.Abs(b-math.Round(b)) > 1e-9 {
			t.Errorf("coset does not preserve lattice: image %v -> frac (%v,%v)", img, a, b)
		}
	}
}

// ------------------------------------------------------------------
// Frieze groups.
// ------------------------------------------------------------------

func TestFriezeGroups(t *testing.T) {
	if len(FriezeGroups()) != 7 {
		t.Fatalf("want 7 frieze groups, got %d", len(FriezeGroups()))
	}
	tests := []struct {
		iuc     string
		pgOrder int
		refl    bool
		glide   bool
		rot     bool
	}{
		{"p1", 1, false, false, false},
		{"p11g", 2, false, true, false},
		{"p1m1", 2, true, false, false},
		{"p2", 2, false, false, true},
		{"p2mg", 4, true, true, true},
		{"p11m", 2, true, false, false},
		{"p2mm", 4, true, false, true},
	}
	for _, tc := range tests {
		g, ok := FriezeGroupByIUC(tc.iuc)
		if !ok {
			t.Errorf("%s not found", tc.iuc)
			continue
		}
		if g.PointGroupOrder() != tc.pgOrder {
			t.Errorf("%s order = %d, want %d", tc.iuc, g.PointGroupOrder(), tc.pgOrder)
		}
		if g.HasReflection() != tc.refl {
			t.Errorf("%s refl = %v, want %v", tc.iuc, g.HasReflection(), tc.refl)
		}
		if g.HasGlideReflection() != tc.glide {
			t.Errorf("%s glide = %v, want %v", tc.iuc, g.HasGlideReflection(), tc.glide)
		}
		if g.HasRotation() != tc.rot {
			t.Errorf("%s rot = %v, want %v", tc.iuc, g.HasRotation(), tc.rot)
		}
	}
}

// ------------------------------------------------------------------
// Orbifold utilities.
// ------------------------------------------------------------------

func TestOrbifoldCost(t *testing.T) {
	all := append(WallpaperIUCNames()[:0:0], WallpaperIUCNames()...)
	_ = all
	for _, g := range WallpaperGroups() {
		c, err := OrbifoldCost(g.Orbifold)
		if err != nil {
			t.Errorf("%s: %v", g.Orbifold, err)
			continue
		}
		if !approx(c, 2) {
			t.Errorf("%s cost = %v, want 2", g.Orbifold, c)
		}
		chi, _ := OrbifoldEulerCharacteristic(g.Orbifold)
		if !approx(chi, 0) {
			t.Errorf("%s euler = %v, want 0", g.Orbifold, chi)
		}
	}
	for _, g := range FriezeGroups() {
		if !IsWallpaperSignature(g.Orbifold, 1e-9) {
			t.Errorf("frieze %s should have cost 2", g.Orbifold)
		}
	}
	if _, err := OrbifoldCost("bad!"); err == nil {
		t.Error("expected error for invalid symbol")
	}
	if orders := GyrationOrders("442"); len(orders) != 3 || orders[0] != 4 || orders[1] != 4 || orders[2] != 2 {
		t.Errorf("GyrationOrders(442) = %v", orders)
	}
	if corners := KaleidoscopeCorners("*632"); len(corners) != 3 {
		t.Errorf("KaleidoscopeCorners(*632) = %v", corners)
	}
	if !HasMirror("*442") || HasMirror("442") {
		t.Error("HasMirror wrong")
	}
}

// ------------------------------------------------------------------
// Regular / Archimedean tilings.
// ------------------------------------------------------------------

func TestInteriorAngles(t *testing.T) {
	tests := []struct {
		n   int
		deg float64
	}{{3, 60}, {4, 90}, {6, 120}, {8, 135}, {12, 150}}
	for _, tc := range tests {
		if got := InteriorAngleDeg(tc.n); !approx(got, tc.deg) {
			t.Errorf("interior angle of %d-gon = %v, want %v", tc.n, got, tc.deg)
		}
	}
}

func TestUniformTilings(t *testing.T) {
	if len(RegularTilings()) != 3 {
		t.Error("want 3 regular tilings")
	}
	if len(SemiregularTilings()) != 8 {
		t.Error("want 8 semiregular tilings")
	}
	if len(ArchimedeanTilings()) != 11 {
		t.Error("want 11 uniform tilings")
	}
	for _, tl := range ArchimedeanTilings() {
		if !tl.Config.IsPlanar(1e-9) {
			t.Errorf("%s vertex angle sum = %v, want 360", tl.Name, tl.Config.AngleSumDeg())
		}
		if !IsArchimedeanConfig(tl.Config) {
			t.Errorf("%s config not recognised as Archimedean", tl.Name)
		}
	}
	for _, tl := range RegularTilings() {
		if !tl.Config.IsRegular(1e-9) {
			t.Errorf("%s should be regular", tl.Name)
		}
	}
}

func TestEnumerateVertexTypes(t *testing.T) {
	types := EnumerateVertexTypes()
	if len(types) != 17 {
		t.Errorf("want 17 vertex types, got %d", len(types))
	}
	for _, v := range types {
		if !v.IsPlanar(1e-9) {
			t.Errorf("vertex type %v does not sum to 360 (got %v)", v, v.AngleSumDeg())
		}
		if v.Degree() < 3 {
			t.Errorf("vertex type %v has degree < 3", v)
		}
	}
	// The regular ones must appear.
	found := map[string]bool{}
	for _, v := range types {
		found[v.String()] = true
	}
	for _, want := range []string{"3.3.3.3.3.3", "4.4.4.4", "6.6.6"} {
		if !found[want] {
			t.Errorf("missing regular vertex type %s", want)
		}
	}
}

func TestDualTilingName(t *testing.T) {
	if d, _ := DualTilingName("triangular"); d != "hexagonal" {
		t.Errorf("dual of triangular = %s", d)
	}
	if d, _ := DualTilingName("square"); d != "square" {
		t.Errorf("square should be self-dual, got %s", d)
	}
}

// ------------------------------------------------------------------
// Substitution tilings.
// ------------------------------------------------------------------

func TestPenroseDeflation(t *testing.T) {
	seed := PenroseSun(1)
	if len(seed) != 10 {
		t.Fatalf("sun seed size = %d, want 10", len(seed))
	}
	area0 := 0.0
	for _, tr := range seed {
		area0 += tr.Area()
	}
	d3 := PenroseP3(3, 1)
	if len(d3) != 130 {
		t.Errorf("P3 after 3 deflations = %d triangles, want 130", len(d3))
	}
	area3 := 0.0
	for _, tr := range d3 {
		area3 += tr.Area()
	}
	if math.Abs(area0-area3) > 1e-9 {
		t.Errorf("area not conserved: %v vs %v", area0, area3)
	}
	// Rhombus assembly yields some tiles.
	if len(PenroseRhombi(d3)) == 0 {
		t.Error("expected some Penrose rhombi")
	}
	// Edge lengths appear in golden ratio: longest/shortest edge of any triangle.
	tr := d3[0]
	e := []float64{tr.A.Distance(tr.B), tr.B.Distance(tr.C), tr.C.Distance(tr.A)}
	max, min := e[0], e[0]
	for _, x := range e {
		max = math.Max(max, x)
		min = math.Min(min, x)
	}
	if !approx(max/min, Phi) {
		t.Errorf("edge ratio = %v, want phi=%v", max/min, Phi)
	}
}

func TestChairSubstitution(t *testing.T) {
	base := NewChair()
	if !approx(base.Area(), 3) {
		t.Errorf("chair area = %v, want 3", base.Area())
	}
	for n := 0; n <= 4; n++ {
		tiles := ChairTiling(n)
		want := 1
		for i := 0; i < n; i++ {
			want *= 4
		}
		if len(tiles) != want {
			t.Errorf("chair n=%d count = %d, want %d", n, len(tiles), want)
		}
		total := 0.0
		for _, c := range tiles {
			total += c.Area()
		}
		if !approx(total, 3) {
			t.Errorf("chair n=%d total area = %v, want 3", n, total)
		}
	}
	// Coverage: the four children partition the parent (no overlap, full cover).
	children := base.Substitute()
	rng := rand.New(rand.NewSource(1))
	parent := base.Polygon()
	samples, covered := 0, 0
	for i := 0; i < 4000; i++ {
		p := NewPoint(rng.Float64()*2, rng.Float64()*2)
		if !PolygonContains(parent, p) {
			continue
		}
		samples++
		count := 0
		for _, ch := range children {
			if PolygonContains(ch.Polygon(), p) {
				count++
			}
		}
		if count == 1 {
			covered++
		}
	}
	if samples == 0 || covered < samples-samples/50 {
		t.Errorf("chair coverage %d/%d not (nearly) exact", covered, samples)
	}
}

func TestPinwheelSubstitution(t *testing.T) {
	base := CanonicalPinwheel()
	if !approx(base.Area(), 1) {
		t.Errorf("pinwheel area = %v, want 1", base.Area())
	}
	for n := 0; n <= 3; n++ {
		tiles := PinwheelTiling(n)
		want := 1
		for i := 0; i < n; i++ {
			want *= 5
		}
		if len(tiles) != want {
			t.Errorf("pinwheel n=%d count = %d, want %d", n, len(tiles), want)
		}
		total := 0.0
		for _, tr := range tiles {
			total += tr.Area()
		}
		if !approx(total, 1) {
			t.Errorf("pinwheel n=%d total area = %v, want 1", n, total)
		}
	}
	// Each child is a right triangle with legs in ratio 2:1.
	for _, ch := range base.Substitute() {
		leg1 := ch.B.Sub(ch.A)
		leg2 := ch.C.Sub(ch.A)
		if math.Abs(leg1.Dot(leg2)) > 1e-9 {
			t.Errorf("child right angle not preserved: dot = %v", leg1.Dot(leg2))
		}
		if !approx(leg1.Norm()/leg2.Norm(), 2) {
			t.Errorf("child leg ratio = %v, want 2", leg1.Norm()/leg2.Norm())
		}
		if !approx(ch.Area(), 0.2) {
			t.Errorf("child area = %v, want 0.2", ch.Area())
		}
	}
}

func TestAffine2(t *testing.T) {
	a := ScaleAffine(2).Compose(TranslateAffine(1, 0))
	if got := a.Apply(NewPoint(1, 1)); !got.ApproxEqual(NewPoint(4, 2), 1e-9) {
		t.Errorf("affine apply = %v, want (4,2)", got)
	}
	if !approx(a.Determinant(), 4) {
		t.Errorf("determinant = %v, want 4", a.Determinant())
	}
}

// ------------------------------------------------------------------
// Examples.
// ------------------------------------------------------------------

func ExampleVertexConfiguration() {
	v := VertexConfiguration{4, 8, 8}
	fmt.Printf("%s sum=%.0f planar=%v\n", v, v.AngleSumDeg(), v.IsPlanar(1e-9))
	// Output: 4.8.8 sum=360 planar=true
}

func ExampleWallpaperGroupByIUC() {
	g, _ := WallpaperGroupByIUC("p4m")
	fmt.Println(g.Orbifold, g.PointGroupOrder(), g.Lattice)
	// Output: *442 8 square
}

func ExampleIsometry_Classify() {
	r := RotationAbout(NewPoint(1, 1), math.Pi/2)
	cl := r.Classify(1e-9)
	fmt.Printf("%v by %.0f deg about (%.0f,%.0f)\n", cl.Kind, Rad2Deg(cl.Angle), cl.Center.X, cl.Center.Y)
	// Output: rotation by 90 deg about (1,1)
}
