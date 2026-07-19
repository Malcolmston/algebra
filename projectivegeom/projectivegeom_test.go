package projectivegeom

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
)

const tol = 1e-9

func approx(a, b float64) bool { return math.Abs(a-b) <= 1e-7 }

func mustAffine(t *testing.T, p Point) (float64, float64) {
	t.Helper()
	x, y, err := p.Affine()
	if err != nil {
		t.Fatalf("Affine(%v): %v", p, err)
	}
	return x, y
}

func TestVec3Ops(t *testing.T) {
	a := NewVec3(1, 2, 3)
	b := NewVec3(4, 5, 6)
	if got := a.Dot(b); got != 32 {
		t.Errorf("Dot = %v want 32", got)
	}
	if got := a.Cross(b); !got.ApproxEqual(NewVec3(-3, 6, -3), tol) {
		t.Errorf("Cross = %v want [-3 6 -3]", got)
	}
	if got := a.Add(b); !got.ApproxEqual(NewVec3(5, 7, 9), tol) {
		t.Errorf("Add = %v", got)
	}
	if got := b.Sub(a); !got.ApproxEqual(NewVec3(3, 3, 3), tol) {
		t.Errorf("Sub = %v", got)
	}
	if !approx(NewVec3(3, 4, 0).Norm(), 5) {
		t.Errorf("Norm wrong")
	}
	if !NewVec3(1, 0, 0).Parallel(NewVec3(-2, 0, 0), tol) {
		t.Errorf("Parallel failed")
	}
	if !approx(NewVec3(1, 0, 0).AngleBetween(NewVec3(0, 1, 0)), math.Pi/2) {
		t.Errorf("angle wrong")
	}
}

func TestJoinMeet(t *testing.T) {
	// x-axis through (0,0) and (1,0).
	l, err := Join(PointFromAffine(0, 0), PointFromAffine(1, 0))
	if err != nil {
		t.Fatal(err)
	}
	if !l.Equal(NewLine(0, 1, 0), tol) {
		t.Errorf("x-axis line = %v want [0 1 0]", l.V)
	}
	// Meet of y=0 and x=2.
	p, err := Meet(NewLine(0, 1, 0), NewLine(1, 0, -2))
	if err != nil {
		t.Fatal(err)
	}
	x, y := mustAffine(t, p)
	if !approx(x, 2) || !approx(y, 0) {
		t.Errorf("meet = (%v,%v) want (2,0)", x, y)
	}
	// Degenerate join.
	if _, err := Join(PointFromAffine(1, 1), PointFromAffine(1, 1)); err != ErrDegenerate {
		t.Errorf("expected ErrDegenerate, got %v", err)
	}
	// Incidence: (2,0) on x-axis.
	if !OnLine(PointFromAffine(2, 0), NewLine(0, 1, 0), tol) {
		t.Errorf("OnLine failed")
	}
}

func TestParallelMeetAtInfinity(t *testing.T) {
	l1 := NewLine(0, 1, 0)  // y = 0
	l2 := NewLine(0, 1, -3) // y = 3
	p, err := Meet(l1, l2)
	if err != nil {
		t.Fatal(err)
	}
	if !p.IsAtInfinity(tol) {
		t.Errorf("parallel lines should meet at infinity, got %v", p.V)
	}
	if !ParallelLines(l1, l2, tol) {
		t.Errorf("ParallelLines false")
	}
}

func TestCollinearConcurrent(t *testing.T) {
	if !Collinear(PointFromAffine(0, 0), PointFromAffine(1, 1), PointFromAffine(2, 2), tol) {
		t.Errorf("collinear points reported non-collinear")
	}
	if Collinear(PointFromAffine(0, 0), PointFromAffine(1, 0), PointFromAffine(0, 1), tol) {
		t.Errorf("triangle reported collinear")
	}
	// Three concurrent lines through origin.
	if !Concurrent(NewLine(1, 0, 0), NewLine(0, 1, 0), NewLine(1, 1, 0), tol) {
		t.Errorf("concurrent lines reported non-concurrent")
	}
}

func TestCrossRatio(t *testing.T) {
	tests := []struct {
		a, b, c, d float64
		want       float64
	}{
		{0, 1, 2, 3, 4.0 / 3.0},
		{0, 1, 2, 4, (2 * 3) / (1 * 4.0)},
		{-1, 1, 0, 2, -1.0 / 3.0},
	}
	for _, tc := range tests {
		a := PointFromAffine(tc.a, 0)
		b := PointFromAffine(tc.b, 0)
		c := PointFromAffine(tc.c, 0)
		d := PointFromAffine(tc.d, 0)
		got, err := CrossRatio(a, b, c, d, tol)
		if err != nil {
			t.Fatal(err)
		}
		if !approx(got, tc.want) {
			t.Errorf("CrossRatio(%v)= %v want %v", tc, got, tc.want)
		}
		// Invariance under a homography.
		h := Compose(Translation(3, -2), Compose(Rotation(0.7), Scaling(2.5)))
		g, err := CrossRatio(h.Apply(a), h.Apply(b), h.Apply(c), h.Apply(d), 1e-6)
		if err != nil {
			t.Fatal(err)
		}
		if !approx(g, tc.want) {
			t.Errorf("cross ratio not invariant: %v vs %v", g, tc.want)
		}
	}
	// Not collinear.
	if _, err := CrossRatio(PointFromAffine(0, 0), PointFromAffine(1, 0),
		PointFromAffine(0, 1), PointFromAffine(1, 1), tol); err != ErrNotCollinear {
		t.Errorf("expected ErrNotCollinear, got %v", err)
	}
}

func TestHarmonicConjugate(t *testing.T) {
	a := PointFromAffine(0, 0)
	b := PointFromAffine(3, 0)
	c := PointFromAffine(1, 0)
	d, err := HarmonicConjugate(a, b, c, tol)
	if err != nil {
		t.Fatal(err)
	}
	x, _ := mustAffine(t, d)
	if !approx(x, -3) {
		t.Errorf("harmonic conjugate x = %v want -3", x)
	}
	if !AreHarmonic(a, b, c, d, 1e-6) {
		t.Errorf("AreHarmonic false for constructed conjugate")
	}
	// Midpoint has its conjugate at infinity.
	m := PointFromAffine(1.5, 0)
	dm, err := HarmonicConjugate(a, b, m, tol)
	if err != nil {
		t.Fatal(err)
	}
	if !dm.IsAtInfinity(1e-6) {
		t.Errorf("conjugate of midpoint should be at infinity, got %v", dm.V)
	}
}

func TestCrossRatioSixValues(t *testing.T) {
	vals := CrossRatioSixValues(2)
	want := []float64{2, 0.5, -1, -1, 2, 0.5}
	for i := range want {
		if !approx(vals[i], want[i]) {
			t.Errorf("six values[%d] = %v want %v", i, vals[i], want[i])
		}
	}
}

func TestHomographyFromPoints(t *testing.T) {
	src := [4]Point{
		PointFromAffine(0, 0), PointFromAffine(1, 0),
		PointFromAffine(1, 1), PointFromAffine(0, 1),
	}
	dst := [4]Point{
		PointFromAffine(2, 1), PointFromAffine(5, 0),
		PointFromAffine(6, 3), PointFromAffine(3, 4),
	}
	h, err := HomographyFromPoints(src, dst)
	if err != nil {
		t.Fatal(err)
	}
	for i := range src {
		got := h.Apply(src[i])
		if !got.Equal(dst[i], 1e-6) {
			gx, gy := mustAffine(t, got)
			dx, dy := mustAffine(t, dst[i])
			t.Errorf("map[%d] = (%v,%v) want (%v,%v)", i, gx, gy, dx, dy)
		}
	}
	// Inverse round-trips.
	inv, err := h.Inverse()
	if err != nil {
		t.Fatal(err)
	}
	rt := Compose(inv, h)
	if !rt.ApproxEqual(IdentityHomography(), 1e-6) {
		t.Errorf("h^-1 * h not identity: %v", rt.M)
	}
}

func TestHomographyLineTransform(t *testing.T) {
	h := Compose(Translation(1, 2), Rotation(0.3))
	l := NewLine(1, -1, 0) // y = x
	pOn := PointFromAffine(2, 2)
	if !OnLine(pOn, l, tol) {
		t.Fatal("setup wrong")
	}
	l2, err := h.ApplyLine(l)
	if err != nil {
		t.Fatal(err)
	}
	if !OnLine(h.Apply(pOn), l2, 1e-6) {
		t.Errorf("incidence not preserved under line transform")
	}
}

func TestHomographyClassification(t *testing.T) {
	if !Translation(1, 2).IsAffine(tol) {
		t.Errorf("translation should be affine")
	}
	if !Similarity(2, 0.5, 1, 1).IsSimilarity(tol) {
		t.Errorf("similarity should be similarity")
	}
	proj := NewHomography(NewMat3(1, 0, 0, 0, 1, 0, 0.5, 0.2, 1))
	if proj.IsAffine(tol) {
		t.Errorf("projective map wrongly reported affine")
	}
	if !ShearX(2).IsAffine(tol) || ShearX(2).IsSimilarity(tol) {
		t.Errorf("shear affine/similarity classification wrong")
	}
}

func TestConicFromFivePoints(t *testing.T) {
	// Unit circle.
	pts := [5]Point{
		PointFromAffine(1, 0), PointFromAffine(0, 1),
		PointFromAffine(-1, 0), PointFromAffine(0, -1),
		PointFromAffine(math.Sqrt2/2, math.Sqrt2/2),
	}
	q, err := ConicFromFivePoints(pts)
	if err != nil {
		t.Fatal(err)
	}
	// Another circle point must lie on it.
	if !q.OnConic(PointFromAffine(0.6, 0.8), 1e-7) {
		t.Errorf("expected (0.6,0.8) on unit circle conic, val=%v", q.Evaluate(PointFromAffine(0.6, 0.8)))
	}
	if q.Classify(1e-9) != ConicEllipse {
		t.Errorf("unit circle should classify as ellipse, got %v", q.Classify(1e-9))
	}
	c, err := q.Center()
	if err != nil {
		t.Fatal(err)
	}
	cx, cy := mustAffine(t, c)
	if !approx(cx, 0) || !approx(cy, 0) {
		t.Errorf("center = (%v,%v) want (0,0)", cx, cy)
	}
}

func TestConicClassify(t *testing.T) {
	tests := []struct {
		q    Conic
		want ConicType
	}{
		{NewConic(1, 0, 1, 0, 0, -1), ConicEllipse},    // circle
		{NewConic(1, 0, -1, 0, 0, -1), ConicHyperbola}, // x^2 - y^2 = 1
		{NewConic(1, 0, 0, 0, -1, 0), ConicParabola},   // y = x^2  -> x^2 - yz = 0
		{NewConic(4, 0, 9, 0, 0, -36), ConicEllipse},   // ellipse
	}
	for i, tc := range tests {
		if got := tc.q.Classify(1e-9); got != tc.want {
			t.Errorf("case %d classify = %v want %v", i, got, tc.want)
		}
	}
}

func TestPolePolar(t *testing.T) {
	circle := UnitCircle()
	// Polar of external point (2,0) is x = 1/2.
	polar := circle.Polar(PointFromAffine(2, 0))
	if !polar.Equal(NewLine(1, 0, -0.5), 1e-9) {
		t.Errorf("polar = %v want x=1/2", polar.V)
	}
	// Pole of that line is back to (2,0).
	pole, err := circle.Pole(polar)
	if err != nil {
		t.Fatal(err)
	}
	if !pole.Equal(PointFromAffine(2, 0), 1e-7) {
		px, py := mustAffine(t, pole)
		t.Errorf("pole = (%v,%v) want (2,0)", px, py)
	}
	// Point on circle: polar is tangent.
	tangent, err := circle.TangentAt(PointFromAffine(1, 0))
	if err != nil {
		t.Fatal(err)
	}
	if !tangent.Equal(NewLine(1, 0, -1), 1e-9) {
		t.Errorf("tangent at (1,0) = %v want x=1", tangent.V)
	}
}

func TestConicIntersectAndTangents(t *testing.T) {
	circle := UnitCircle()
	pts := circle.IntersectLine(NewLine(1, 0, 0)) // x = 0 meets circle at (0,+-1)
	if len(pts) != 2 {
		t.Fatalf("expected 2 intersections, got %d", len(pts))
	}
	// Tangent from external point (2,0): two contacts at x=1/2, y=+-sqrt(3)/2.
	tans, err := circle.TangentsFrom(PointFromAffine(2, 0))
	if err != nil {
		t.Fatal(err)
	}
	for _, l := range tans {
		if !OnLine(PointFromAffine(2, 0), l, 1e-7) {
			t.Errorf("tangent line %v does not pass through (2,0)", l.V)
		}
	}
	// A secant line meeting the circle in two known points.
	sec := circle.IntersectLine(NewLine(0, 1, 0)) // y = 0 -> (+-1,0)
	if len(sec) != 2 {
		t.Fatalf("expected 2 secant points, got %d", len(sec))
	}
	xs := []float64{}
	for _, p := range sec {
		x, _ := mustAffine(t, p)
		xs = append(xs, x)
	}
	if !((approx(xs[0], 1) && approx(xs[1], -1)) || (approx(xs[0], -1) && approx(xs[1], 1))) {
		t.Errorf("secant points x = %v want {-1,1}", xs)
	}
}

func TestConicTransform(t *testing.T) {
	circle := UnitCircle()
	h := Compose(Translation(3, 0), Scaling(2)) // scale by 2 then translate
	q2, err := circle.Transform(h)
	if err != nil {
		t.Fatal(err)
	}
	// Image of (1,0) is (5,0), must lie on transformed conic.
	if !q2.OnConic(PointFromAffine(5, 0), 1e-7) {
		t.Errorf("transformed conic missing image point")
	}
	center, err := q2.Center()
	if err != nil {
		t.Fatal(err)
	}
	cx, cy := mustAffine(t, center)
	if !approx(cx, 3) || !approx(cy, 0) {
		t.Errorf("transformed center = (%v,%v) want (3,0)", cx, cy)
	}
}

func TestConjugatePoints(t *testing.T) {
	circle := UnitCircle()
	// (2,0) and its polar foot (0.5,0) are conjugate.
	if !ConjugatePoints(circle, PointFromAffine(2, 0), PointFromAffine(0.5, 0), 1e-9) {
		t.Errorf("expected conjugate points")
	}
	if ConjugatePoints(circle, PointFromAffine(2, 0), PointFromAffine(2, 2), 1e-9) {
		t.Errorf("unexpected conjugacy")
	}
}

func TestPappus(t *testing.T) {
	cfg := PappusConfig{
		A:  PointFromAffine(0, 0),
		B:  PointFromAffine(1, 0),
		C:  PointFromAffine(2, 0),
		A2: PointFromAffine(0, 1),
		B2: PointFromAffine(2, 1),
		C2: PointFromAffine(4, 1),
	}
	if !cfg.Holds(1e-7) {
		t.Errorf("Pappus configuration should hold")
	}
	if _, err := cfg.PappusLine(); err != nil {
		t.Errorf("PappusLine error: %v", err)
	}
}

func TestDesargues(t *testing.T) {
	tp := TrianglePair{
		T1: Triangle{PointFromAffine(1, 0), PointFromAffine(0, 1), PointFromAffine(1, 1)},
		T2: Triangle{PointFromAffine(3, 0), PointFromAffine(0, 4), PointFromAffine(5, 5)},
	}
	if !tp.IsPerspectiveFromPoint(1e-7) {
		t.Errorf("triangles should be perspective from a point")
	}
	center, err := tp.PerspectiveCenter()
	if err != nil {
		t.Fatal(err)
	}
	cx, cy := mustAffine(t, center)
	if !approx(cx, 0) || !approx(cy, 0) {
		t.Errorf("perspective center = (%v,%v) want (0,0)", cx, cy)
	}
	if !tp.IsPerspectiveFromLine(1e-6) {
		t.Errorf("Desargues: should also be perspective from a line")
	}
}

func TestFixedPoints(t *testing.T) {
	// Pure rotation about origin has the origin (finite) and the two circular
	// points at infinity as fixed points; the real fixed point is the origin.
	h := RotationAbout(0.5, 0, 0)
	fps := h.FixedPointsReal(1e-7)
	foundOrigin := false
	for _, p := range fps {
		if p.Equal(PointFromAffine(0, 0), 1e-6) {
			foundOrigin = true
		}
	}
	if !foundOrigin {
		t.Errorf("rotation about origin should fix the origin, got %v", fps)
	}
}

func TestRP3PlanePoint(t *testing.T) {
	pl, err := PlaneThrough3Points(
		SPointFromAffine(0, 0, 0), SPointFromAffine(1, 0, 0), SPointFromAffine(0, 1, 0))
	if err != nil {
		t.Fatal(err)
	}
	if !pl.Equal(NewSPlane(0, 0, 1, 0), 1e-9) {
		t.Errorf("plane = %v want z=0", pl.V)
	}
	// Three coordinate planes meet at the origin.
	p, err := PointOf3Planes(NewSPlane(1, 0, 0, 0), NewSPlane(0, 1, 0, 0), NewSPlane(0, 0, 1, 0))
	if err != nil {
		t.Fatal(err)
	}
	x, y, z, err := p.Affine()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(x, 0) || !approx(y, 0) || !approx(z, 0) {
		t.Errorf("meet of planes = (%v,%v,%v) want origin", x, y, z)
	}
	if !Coplanar(SPointFromAffine(0, 0, 0), SPointFromAffine(1, 0, 0),
		SPointFromAffine(0, 1, 0), SPointFromAffine(2, 3, 0), 1e-7) {
		t.Errorf("coplanar points reported non-coplanar")
	}
	if Coplanar(SPointFromAffine(0, 0, 0), SPointFromAffine(1, 0, 0),
		SPointFromAffine(0, 1, 0), SPointFromAffine(0, 0, 1), 1e-7) {
		t.Errorf("tetrahedron reported coplanar")
	}
}

func TestPluecker(t *testing.T) {
	l, err := LineFromPoints(SPointFromAffine(0, 0, 0), SPointFromAffine(1, 0, 0))
	if err != nil {
		t.Fatal(err)
	}
	if !l.IsLine(1e-9) {
		t.Errorf("valid line failed Grassmann relation")
	}
	if !l.ContainsPoint(SPointFromAffine(5, 0, 0), 1e-7) {
		t.Errorf("x-axis should contain (5,0,0)")
	}
	if l.ContainsPoint(SPointFromAffine(0, 1, 0), 1e-7) {
		t.Errorf("x-axis should not contain (0,1,0)")
	}
	// Same line as intersection of y=0 and z=0.
	l2, err := LineFromPlanes(NewSPlane(0, 1, 0, 0), NewSPlane(0, 0, 1, 0))
	if err != nil {
		t.Fatal(err)
	}
	if !l.Equal(l2, 1e-7) {
		t.Errorf("line from points %v != line from planes %v", l, l2)
	}
	// Two skew lines do not meet; two intersecting lines do.
	xAxis, _ := LineFromPoints(SPointFromAffine(0, 0, 0), SPointFromAffine(1, 0, 0))
	yAxis, _ := LineFromPoints(SPointFromAffine(0, 0, 0), SPointFromAffine(0, 1, 0))
	skew, _ := LineFromPoints(SPointFromAffine(0, 0, 1), SPointFromAffine(1, 1, 1))
	if !LinesMeet(xAxis, yAxis, 1e-7) {
		t.Errorf("axes should meet at origin")
	}
	if LinesMeet(xAxis, skew, 1e-7) {
		t.Errorf("skew lines should not meet")
	}
}

func TestCamera(t *testing.T) {
	cam := CanonicalCamera()
	img := cam.Project(SPointFromAffine(1, 2, 3))
	u, v := mustAffine(t, img)
	if !approx(u, 1.0/3.0) || !approx(v, 2.0/3.0) {
		t.Errorf("projection = (%v,%v) want (1/3,2/3)", u, v)
	}
	center := cam.Center()
	if !center.Equal(SPointFromAffine(0, 0, 0), 1e-9) {
		t.Errorf("canonical camera center = %v want origin", center.V)
	}
	// Camera from K,R,t with identity K,R and translation.
	cam2 := CameraFromKRt(Identity3(), Identity3(), NewVec3(0, 0, 5))
	img2 := cam2.Project(SPointFromAffine(1, 1, 5))
	u2, v2 := mustAffine(t, img2)
	if !approx(u2, 0.1) || !approx(v2, 0.1) {
		t.Errorf("translated projection = (%v,%v) want (0.1,0.1)", u2, v2)
	}
}

func TestMat3InverseDet(t *testing.T) {
	m := NewMat3(2, 0, 1, 1, 3, 2, 1, 0, 3)
	inv, ok := m.Inverse()
	if !ok {
		t.Fatal("expected invertible")
	}
	if !m.Mul(inv).ApproxEqual(Identity3(), 1e-9) {
		t.Errorf("m*inv != I")
	}
	if !approx(m.Det(), m.Mul(Identity3()).Det()) {
		t.Errorf("det consistency")
	}
	// Singular matrix.
	if _, ok := NewMat3(1, 2, 3, 2, 4, 6, 0, 1, 0).Inverse(); ok {
		t.Errorf("singular matrix reported invertible")
	}
}

func TestMat4Inverse(t *testing.T) {
	m := Mat4{{1, 2, 0, 0}, {0, 1, 0, 0}, {0, 0, 2, 1}, {0, 0, 0, 1}}
	inv, ok := m.Inverse()
	if !ok {
		t.Fatal("expected invertible")
	}
	if !m.Mul(inv).ApproxEqual(Identity4(), 1e-9) {
		t.Errorf("mat4 inverse wrong")
	}
}

func TestDuality(t *testing.T) {
	p := NewPoint(2, 3, 1)
	l := DualOfPoint(p)
	if l.V != p.V {
		t.Errorf("dual should share coordinates")
	}
	if DualOfLine(l).V != p.V {
		t.Errorf("double dual mismatch")
	}
}

func TestRandomDeterminism(t *testing.T) {
	r1 := rand.New(rand.NewSource(42))
	r2 := rand.New(rand.NewSource(42))
	p1 := RandomPoint(r1)
	p2 := RandomPoint(r2)
	if p1.V != p2.V {
		t.Errorf("same seed should give same point: %v vs %v", p1.V, p2.V)
	}
	// Random homography preserves cross ratio of random collinear points.
	r := rand.New(rand.NewSource(7))
	h, ok := RandomHomography(r)
	if !ok {
		t.Skip("degenerate random homography")
	}
	a := PointFromAffine(0, 0)
	b := PointFromAffine(1, 0)
	c := PointFromAffine(2, 0)
	d := PointFromAffine(5, 0)
	cr0, _ := CrossRatio(a, b, c, d, tol)
	cr1, err := CrossRatio(h.Apply(a), h.Apply(b), h.Apply(c), h.Apply(d), 1e-6)
	if err != nil {
		t.Fatal(err)
	}
	if !approx(cr0, cr1) {
		t.Errorf("random homography changed cross ratio: %v vs %v", cr0, cr1)
	}
}

func ExampleCrossRatio() {
	a := PointFromAffine(0, 0)
	b := PointFromAffine(1, 0)
	c := PointFromAffine(2, 0)
	d := PointFromAffine(3, 0)
	cr, _ := CrossRatio(a, b, c, d, 1e-9)
	fmt.Printf("%.4f\n", cr)
	// Output: 1.3333
}

func ExampleMeet() {
	// Intersect the lines y = 0 and x = 2.
	p, _ := Meet(NewLine(0, 1, 0), NewLine(1, 0, -2))
	x, y, _ := p.Affine()
	fmt.Printf("(%.0f, %.0f)\n", x+0, y+0)
	// Output: (2, 0)
}
