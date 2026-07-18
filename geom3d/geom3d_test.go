package geom3d

import (
	"math"
	"testing"
)

func TestPlaneFromPoints(t *testing.T) {
	pl, ok := PlaneFromPoints(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0})
	if !ok {
		t.Fatal("expected valid plane")
	}
	if !pl.Normal.Equal(Vec3{0, 0, 1}, 1e-9) {
		t.Errorf("normal = %v", pl.Normal)
	}
	if math.Abs(pl.D) > testEps {
		t.Errorf("D = %v want 0", pl.D)
	}
	// Collinear -> false.
	if _, ok := PlaneFromPoints(Vec3{0, 0, 0}, Vec3{1, 1, 1}, Vec3{2, 2, 2}); ok {
		t.Errorf("expected collinear failure")
	}
}

func TestPlaneDistanceAndProject(t *testing.T) {
	pl := PlaneFromPointNormal(Vec3{0, 0, 5}, Vec3{0, 0, 2}) // z = 5 plane
	if got := pl.SignedDistance(Vec3{1, 2, 8}); math.Abs(got-3) > testEps {
		t.Errorf("signed dist = %v want 3", got)
	}
	if got := pl.SignedDistance(Vec3{0, 0, 1}); math.Abs(got+4) > testEps {
		t.Errorf("signed dist = %v want -4", got)
	}
	if got := pl.Distance(Vec3{0, 0, 1}); math.Abs(got-4) > testEps {
		t.Errorf("dist = %v want 4", got)
	}
	if got := pl.ClosestPoint(Vec3{1, 2, 8}); !got.Equal(Vec3{1, 2, 5}, 1e-9) {
		t.Errorf("closest = %v", got)
	}
	if pl.Side(Vec3{0, 0, 9}, geom3dEps) != 1 || pl.Side(Vec3{0, 0, 1}, geom3dEps) != -1 || pl.Side(Vec3{0, 0, 5}, geom3dEps) != 0 {
		t.Errorf("Side wrong")
	}
	if !pl.Contains(Vec3{7, -3, 5}, 1e-9) {
		t.Errorf("Contains failed")
	}
}

func TestLine3ClosestAndDistance(t *testing.T) {
	l := LineThrough(Vec3{0, 0, 0}, Vec3{0, 0, 1}) // z-axis
	if got := l.Distance(Vec3{3, 4, 10}); math.Abs(got-5) > testEps {
		t.Errorf("line dist = %v want 5", got)
	}
	if got := l.ClosestPoint(Vec3{3, 4, 10}); !got.Equal(Vec3{0, 0, 10}, 1e-9) {
		t.Errorf("closest = %v", got)
	}
	if got := l.At(2); !got.Equal(Vec3{0, 0, 2}, 1e-9) {
		t.Errorf("At = %v", got)
	}
}

func TestClosestPointsBetweenLines(t *testing.T) {
	// x-axis through z=0 and y-axis through z=1: closest points on the z gap.
	l1 := Line3{Point: Vec3{0, 0, 0}, Dir: Vec3{1, 0, 0}}
	l2 := Line3{Point: Vec3{0, 0, 1}, Dir: Vec3{0, 1, 0}}
	p1, p2, ok := ClosestPointsBetweenLines(l1, l2)
	if !ok {
		t.Fatal("expected non-parallel")
	}
	if !p1.Equal(Vec3{0, 0, 0}, 1e-9) || !p2.Equal(Vec3{0, 0, 1}, 1e-9) {
		t.Errorf("closest pts = %v %v", p1, p2)
	}
	if math.Abs(p1.Distance(p2)-1) > testEps {
		t.Errorf("gap = %v want 1", p1.Distance(p2))
	}
	// Parallel lines flagged.
	if _, _, ok := ClosestPointsBetweenLines(l1, Line3{Point: Vec3{0, 5, 0}, Dir: Vec3{2, 0, 0}}); ok {
		t.Errorf("expected parallel = false")
	}
}

func TestRayPlaneIntersect(t *testing.T) {
	pl := PlaneFromPointNormal(Vec3{0, 0, 0}, Vec3{0, 0, 1})
	r := NewRay(Vec3{0, 0, 5}, Vec3{0, 0, -1})
	tt, ok := r.IntersectPlane(pl)
	if !ok || math.Abs(tt-5) > testEps {
		t.Errorf("t = %v ok=%v want 5", tt, ok)
	}
	if !r.At(tt).Equal(Vec3{0, 0, 0}, 1e-9) {
		t.Errorf("hit = %v", r.At(tt))
	}
	// Parallel ray.
	if _, ok := NewRay(Vec3{0, 0, 5}, Vec3{1, 0, 0}).IntersectPlane(pl); ok {
		t.Errorf("expected parallel miss")
	}
}

func TestBarycentric(t *testing.T) {
	a := Vec3{0, 0, 0}
	b := Vec3{1, 0, 0}
	c := Vec3{0, 1, 0}
	// Centroid.
	cen := a.Add(b).Add(c).Scale(1.0 / 3)
	u, v, w, ok := Barycentric(cen, a, b, c)
	if !ok || math.Abs(u-1.0/3) > 1e-9 || math.Abs(v-1.0/3) > 1e-9 || math.Abs(w-1.0/3) > 1e-9 {
		t.Errorf("centroid bary = %v %v %v", u, v, w)
	}
	// Vertex a.
	u, v, w, _ = Barycentric(a, a, b, c)
	if math.Abs(u-1) > 1e-9 || math.Abs(v) > 1e-9 || math.Abs(w) > 1e-9 {
		t.Errorf("vertex bary = %v %v %v", u, v, w)
	}
	if !PointInTriangle(cen, a, b, c) {
		t.Errorf("centroid should be inside")
	}
	if PointInTriangle(Vec3{2, 2, 0}, a, b, c) {
		t.Errorf("outside point reported inside")
	}
}

func TestRayTriangle(t *testing.T) {
	a := Vec3{0, 0, 0}
	b := Vec3{1, 0, 0}
	c := Vec3{0, 1, 0}
	// Ray straight down onto the triangle interior.
	r := NewRay(Vec3{0.25, 0.25, 1}, Vec3{0, 0, -1})
	tt, u, v, ok := RayTriangle(r, a, b, c)
	if !ok {
		t.Fatal("expected hit")
	}
	if math.Abs(tt-1) > 1e-9 {
		t.Errorf("t = %v want 1", tt)
	}
	hit := r.At(tt)
	if !hit.Equal(Vec3{0.25, 0.25, 0}, 1e-9) {
		t.Errorf("hit = %v", hit)
	}
	if math.Abs(u-0.25) > 1e-9 || math.Abs(v-0.25) > 1e-9 {
		t.Errorf("u,v = %v %v", u, v)
	}
	// Miss (outside triangle).
	if _, _, _, ok := RayTriangle(NewRay(Vec3{2, 2, 1}, Vec3{0, 0, -1}), a, b, c); ok {
		t.Errorf("expected miss")
	}
	// Parallel ray.
	if _, _, _, ok := RayTriangle(NewRay(Vec3{0.25, 0.25, 0}, Vec3{1, 0, 0}), a, b, c); ok {
		t.Errorf("expected parallel miss")
	}
	// Behind the origin.
	if _, _, _, ok := RayTriangle(NewRay(Vec3{0.25, 0.25, 1}, Vec3{0, 0, 1}), a, b, c); ok {
		t.Errorf("expected behind miss")
	}
}

func TestAABB(t *testing.T) {
	bx := NewAABB(Vec3{2, 2, 2}, Vec3{-1, -1, -1})
	if !bx.Min.Equal(Vec3{-1, -1, -1}, 1e-12) || !bx.Max.Equal(Vec3{2, 2, 2}, 1e-12) {
		t.Errorf("NewAABB did not normalize: %v", bx)
	}
	if !bx.Center().Equal(Vec3{0.5, 0.5, 0.5}, 1e-12) {
		t.Errorf("center = %v", bx.Center())
	}
	if !bx.Size().Equal(Vec3{3, 3, 3}, 1e-12) {
		t.Errorf("size = %v", bx.Size())
	}
	if !bx.Contains(Vec3{0, 0, 0}) || bx.Contains(Vec3{3, 0, 0}) {
		t.Errorf("Contains wrong")
	}
	if got := bx.ClosestPoint(Vec3{5, 0, -5}); !got.Equal(Vec3{2, 0, -1}, 1e-12) {
		t.Errorf("closest = %v", got)
	}
	other := NewAABB(Vec3{1, 1, 1}, Vec3{3, 3, 3})
	if !bx.IntersectsAABB(other) {
		t.Errorf("should intersect")
	}
	far := NewAABB(Vec3{10, 10, 10}, Vec3{11, 11, 11})
	if bx.IntersectsAABB(far) {
		t.Errorf("should not intersect")
	}
	u := bx.Union(far)
	if !u.Min.Equal(Vec3{-1, -1, -1}, 1e-12) || !u.Max.Equal(Vec3{11, 11, 11}, 1e-12) {
		t.Errorf("union = %v", u)
	}
	fromPts, ok := AABBFromPoints(Vec3{1, 2, 3}, Vec3{-1, 5, 0}, Vec3{4, -2, 1})
	if !ok || !fromPts.Min.Equal(Vec3{-1, -2, 0}, 1e-12) || !fromPts.Max.Equal(Vec3{4, 5, 3}, 1e-12) {
		t.Errorf("fromPts = %v", fromPts)
	}
	if _, ok := AABBFromPoints(); ok {
		t.Errorf("empty points should fail")
	}
}

func TestAABBRay(t *testing.T) {
	bx := NewAABB(Vec3{-1, -1, -1}, Vec3{1, 1, 1})
	r := NewRay(Vec3{-5, 0, 0}, Vec3{1, 0, 0})
	tmin, tmax, ok := bx.IntersectRay(r)
	if !ok || math.Abs(tmin-4) > 1e-9 || math.Abs(tmax-6) > 1e-9 {
		t.Errorf("ray box = %v %v %v", tmin, tmax, ok)
	}
	// Origin inside -> tmin clamped to 0.
	tmin, _, ok = bx.IntersectRay(NewRay(Vec3{0, 0, 0}, Vec3{1, 0, 0}))
	if !ok || tmin != 0 {
		t.Errorf("inside tmin = %v", tmin)
	}
	// Miss.
	if _, _, ok := bx.IntersectRay(NewRay(Vec3{-5, 5, 0}, Vec3{1, 0, 0})); ok {
		t.Errorf("expected miss")
	}
	// Pointing away.
	if _, _, ok := bx.IntersectRay(NewRay(Vec3{-5, 0, 0}, Vec3{-1, 0, 0})); ok {
		t.Errorf("expected behind miss")
	}
}

func TestSphere(t *testing.T) {
	s := Sphere{Center: Vec3{0, 0, 0}, Radius: 2}
	if !s.Contains(Vec3{1, 1, 1}) || s.Contains(Vec3{2, 2, 2}) {
		t.Errorf("Contains wrong")
	}
	if !s.IntersectsSphere(Sphere{Center: Vec3{3, 0, 0}, Radius: 2}) {
		t.Errorf("spheres should touch/intersect")
	}
	if s.IntersectsSphere(Sphere{Center: Vec3{5, 0, 0}, Radius: 2}) {
		t.Errorf("spheres should not intersect")
	}
	// Sphere vs AABB.
	bx := NewAABB(Vec3{2, -1, -1}, Vec3{4, 1, 1})
	if !s.IntersectsAABB(bx) {
		t.Errorf("sphere should touch box at x=2")
	}
	if !bx.IntersectsSphere(s) {
		t.Errorf("box-sphere symmetric failed")
	}
	if s.IntersectsAABB(NewAABB(Vec3{10, 10, 10}, Vec3{11, 11, 11})) {
		t.Errorf("far box should not intersect")
	}
	// Ray.
	tmin, tmax, ok := s.IntersectRay(NewRay(Vec3{-5, 0, 0}, Vec3{1, 0, 0}))
	if !ok || math.Abs(tmin-3) > 1e-9 || math.Abs(tmax-7) > 1e-9 {
		t.Errorf("sphere ray = %v %v %v", tmin, tmax, ok)
	}
	if _, _, ok := s.IntersectRay(NewRay(Vec3{-5, 5, 0}, Vec3{1, 0, 0})); ok {
		t.Errorf("expected ray miss")
	}
}

func BenchmarkRayTriangle(b *testing.B) {
	r := NewRay(Vec3{0.25, 0.25, 1}, Vec3{0, 0, -1})
	a := Vec3{0, 0, 0}
	v1 := Vec3{1, 0, 0}
	v2 := Vec3{0, 1, 0}
	b.ReportAllocs()
	var acc float64
	for i := 0; i < b.N; i++ {
		t, _, _, ok := RayTriangle(r, a, v1, v2)
		if ok {
			acc += t
		}
	}
	_ = acc
}
