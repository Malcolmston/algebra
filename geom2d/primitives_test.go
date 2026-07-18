package geom2d

import (
	"math"
	"testing"
)

const testTol = 1e-9

func approx(a, b float64) bool { return math.Abs(a-b) <= testTol }

func pointApprox(p, q Point2) bool {
	return approx(p.X, q.X) && approx(p.Y, q.Y)
}

func vecApprox(u, v Vec2) bool {
	return approx(u.X, v.X) && approx(u.Y, v.Y)
}

func TestVecArithmetic(t *testing.T) {
	u := NewVec2(3, 4)
	v := NewVec2(1, -2)
	if got := u.Add(v); !vecApprox(got, Vec2{4, 2}) {
		t.Errorf("Add = %v", got)
	}
	if got := u.Sub(v); !vecApprox(got, Vec2{2, 6}) {
		t.Errorf("Sub = %v", got)
	}
	if got := u.Scale(2); !vecApprox(got, Vec2{6, 8}) {
		t.Errorf("Scale = %v", got)
	}
	if got := u.Neg(); !vecApprox(got, Vec2{-3, -4}) {
		t.Errorf("Neg = %v", got)
	}
	if got := u.Dot(v); !approx(got, 3*1+4*-2) {
		t.Errorf("Dot = %v", got)
	}
	if got := u.Cross(v); !approx(got, 3*-2-4*1) {
		t.Errorf("Cross = %v", got)
	}
	if got := u.Hadamard(v); !vecApprox(got, Vec2{3, -8}) {
		t.Errorf("Hadamard = %v", got)
	}
}

func TestNorms(t *testing.T) {
	u := NewVec2(3, 4)
	if !approx(u.Norm(), 5) {
		t.Errorf("Norm = %v", u.Norm())
	}
	if !approx(u.NormSq(), 25) {
		t.Errorf("NormSq = %v", u.NormSq())
	}
	unit, n := u.Normalize()
	if !approx(n, 5) || !vecApprox(unit, Vec2{0.6, 0.8}) {
		t.Errorf("Normalize = %v, %v", unit, n)
	}
	if !vecApprox(u.Unit(), Vec2{0.6, 0.8}) {
		t.Errorf("Unit = %v", u.Unit())
	}
	if z, n := (Vec2{}).Normalize(); !vecApprox(z, Vec2{}) || n != 0 {
		t.Errorf("zero Normalize = %v, %v", z, n)
	}
	if !(Vec2{1e-12, 0}).IsZero(1e-9) {
		t.Errorf("IsZero should be true")
	}
}

func TestPerpAndAngle(t *testing.T) {
	u := NewVec2(1, 0)
	if !vecApprox(u.Perp(), Vec2{0, 1}) {
		t.Errorf("Perp = %v", u.Perp())
	}
	if !vecApprox(u.PerpCW(), Vec2{0, -1}) {
		t.Errorf("PerpCW = %v", u.PerpCW())
	}
	if !approx(NewVec2(0, 1).Angle(), math.Pi/2) {
		t.Errorf("Angle = %v", NewVec2(0, 1).Angle())
	}
	if !approx(NewVec2(-1, 0).Angle(), math.Pi) {
		t.Errorf("Angle(-1,0) = %v", NewVec2(-1, 0).Angle())
	}
}

func TestRotate(t *testing.T) {
	u := NewVec2(1, 0)
	got := u.Rotate(math.Pi / 2)
	if !vecApprox(got, Vec2{0, 1}) {
		t.Errorf("Rotate 90 = %v", got)
	}
	got = u.Rotate(math.Pi)
	if !vecApprox(got, Vec2{-1, 0}) {
		t.Errorf("Rotate 180 = %v", got)
	}
	// Rotation preserves length.
	w := NewVec2(3, -4).Rotate(0.7)
	if !approx(w.Norm(), 5) {
		t.Errorf("Rotate length = %v", w.Norm())
	}
}

func TestAngleBetween(t *testing.T) {
	cases := []struct {
		u, v Vec2
		want float64
	}{
		{Vec2{1, 0}, Vec2{0, 1}, math.Pi / 2},
		{Vec2{1, 0}, Vec2{1, 0}, 0},
		{Vec2{1, 0}, Vec2{-1, 0}, math.Pi},
		{Vec2{1, 0}, Vec2{1, 1}, math.Pi / 4},
	}
	for _, c := range cases {
		if got := c.u.AngleBetween(c.v); !approx(got, c.want) {
			t.Errorf("AngleBetween(%v,%v) = %v want %v", c.u, c.v, got, c.want)
		}
	}
	if got := NewVec2(1, 0).SignedAngle(NewVec2(0, 1)); !approx(got, math.Pi/2) {
		t.Errorf("SignedAngle ccw = %v", got)
	}
	if got := NewVec2(1, 0).SignedAngle(NewVec2(0, -1)); !approx(got, -math.Pi/2) {
		t.Errorf("SignedAngle cw = %v", got)
	}
}

func TestProjectReject(t *testing.T) {
	u := NewVec2(2, 3)
	v := NewVec2(1, 0)
	if !vecApprox(u.Project(v), Vec2{2, 0}) {
		t.Errorf("Project = %v", u.Project(v))
	}
	if !vecApprox(u.Reject(v), Vec2{0, 3}) {
		t.Errorf("Reject = %v", u.Reject(v))
	}
	if !approx(u.ScalarProject(v), 2) {
		t.Errorf("ScalarProject = %v", u.ScalarProject(v))
	}
	// Project + Reject reconstructs u.
	if !vecApprox(u.Project(v).Add(u.Reject(v)), u) {
		t.Errorf("Project+Reject != u")
	}
}

func TestReflect(t *testing.T) {
	// Reflect across X axis.
	u := NewVec2(2, 3)
	if !vecApprox(u.Reflect(NewVec2(1, 0)), Vec2{2, -3}) {
		t.Errorf("Reflect axis = %v", u.Reflect(NewVec2(1, 0)))
	}
	// Bounce off surface with upward normal reverses Y.
	got := NewVec2(1, -1).ReflectAcrossNormal(NewVec2(0, 1))
	if !vecApprox(got, Vec2{1, 1}) {
		t.Errorf("ReflectAcrossNormal = %v", got)
	}
}

func TestLerpAndMidpoint(t *testing.T) {
	if !vecApprox(NewVec2(0, 0).Lerp(NewVec2(10, 20), 0.5), Vec2{5, 10}) {
		t.Errorf("Vec Lerp wrong")
	}
	p, q := NewPoint2(0, 0), NewPoint2(4, 8)
	if !pointApprox(Lerp(p, q, 0.25), Point2{1, 2}) {
		t.Errorf("Lerp wrong")
	}
	if !pointApprox(Midpoint(p, q), Point2{2, 4}) {
		t.Errorf("Midpoint wrong")
	}
	if !pointApprox(Centroid(p, q, NewPoint2(2, 1)), Point2{2, 3}) {
		t.Errorf("Centroid wrong")
	}
}

func TestDistance(t *testing.T) {
	p, q := NewPoint2(0, 0), NewPoint2(3, 4)
	if !approx(Distance(p, q), 5) {
		t.Errorf("Distance = %v", Distance(p, q))
	}
	if !approx(DistanceSq(p, q), 25) {
		t.Errorf("DistanceSq = %v", DistanceSq(p, q))
	}
	if !approx(p.DistanceTo(q), 5) || !approx(p.DistanceSqTo(q), 25) {
		t.Errorf("Point distance methods wrong")
	}
}

func TestOrientation(t *testing.T) {
	a, b := NewPoint2(0, 0), NewPoint2(1, 0)
	if Orientation(a, b, NewPoint2(0, 1), testTol) != 1 {
		t.Errorf("expected CCW")
	}
	if Orientation(a, b, NewPoint2(0, -1), testTol) != -1 {
		t.Errorf("expected CW")
	}
	if Orientation(a, b, NewPoint2(2, 0), testTol) != 0 {
		t.Errorf("expected collinear")
	}
	if !CCW(a, b, NewPoint2(0, 1), testTol) {
		t.Errorf("CCW true expected")
	}
	if !Collinear(a, b, NewPoint2(5, 0), testTol) {
		t.Errorf("Collinear expected")
	}
}

func TestAreas(t *testing.T) {
	a, b, c := NewPoint2(0, 0), NewPoint2(4, 0), NewPoint2(0, 3)
	if !approx(SignedArea2(a, b, c), 12) {
		t.Errorf("SignedArea2 = %v", SignedArea2(a, b, c))
	}
	if !approx(TriangleArea(a, b, c), 6) {
		t.Errorf("TriangleArea = %v", TriangleArea(a, b, c))
	}
	// Clockwise ordering flips sign.
	if !approx(SignedArea2(a, c, b), -12) {
		t.Errorf("SignedArea2 cw = %v", SignedArea2(a, c, b))
	}
}

func TestAngleAt(t *testing.T) {
	a, b, c := NewPoint2(1, 0), NewPoint2(0, 0), NewPoint2(0, 1)
	if !approx(AngleAt(a, b, c), math.Pi/2) {
		t.Errorf("AngleAt = %v", AngleAt(a, b, c))
	}
}

func TestSlope(t *testing.T) {
	if m, ok := Slope(NewPoint2(0, 0), NewPoint2(2, 4)); !ok || !approx(m, 2) {
		t.Errorf("Slope = %v, %v", m, ok)
	}
	if _, ok := Slope(NewPoint2(1, 0), NewPoint2(1, 5)); ok {
		t.Errorf("vertical Slope should report false")
	}
}

func TestRotateAround(t *testing.T) {
	p := NewPoint2(2, 0)
	c := NewPoint2(1, 0)
	got := RotateAround(p, c, math.Pi/2)
	if !pointApprox(got, Point2{1, 1}) {
		t.Errorf("RotateAround = %v", got)
	}
}

func TestClosestPointAndDistances(t *testing.T) {
	l := Line{A: NewPoint2(0, 0), B: NewPoint2(10, 0)}
	p := NewPoint2(3, 5)
	if !pointApprox(ClosestPointOnLine(p, l), Point2{3, 0}) {
		t.Errorf("ClosestPointOnLine = %v", ClosestPointOnLine(p, l))
	}
	if !approx(PointLineDistance(p, l), 5) {
		t.Errorf("PointLineDistance = %v", PointLineDistance(p, l))
	}
	if !approx(SignedPointLineDistance(p, l), 5) {
		t.Errorf("SignedPointLineDistance = %v", SignedPointLineDistance(p, l))
	}
	if !approx(SignedPointLineDistance(NewPoint2(3, -5), l), -5) {
		t.Errorf("SignedPointLineDistance below = %v", SignedPointLineDistance(NewPoint2(3, -5), l))
	}

	s := Segment{A: NewPoint2(0, 0), B: NewPoint2(10, 0)}
	// Foot beyond B clamps to B.
	if !pointApprox(ClosestPointOnSegment(NewPoint2(20, 5), s), Point2{10, 0}) {
		t.Errorf("clamp to B wrong")
	}
	// Foot before A clamps to A.
	if !pointApprox(ClosestPointOnSegment(NewPoint2(-5, 5), s), Point2{0, 0}) {
		t.Errorf("clamp to A wrong")
	}
	// Interior projection.
	if !pointApprox(ClosestPointOnSegment(p, s), Point2{3, 0}) {
		t.Errorf("interior projection wrong")
	}
	if !approx(PointSegmentDistance(NewPoint2(20, 0), s), 10) {
		t.Errorf("PointSegmentDistance beyond = %v", PointSegmentDistance(NewPoint2(20, 0), s))
	}
	if !approx(PointSegmentDistance(p, s), 5) {
		t.Errorf("PointSegmentDistance interior = %v", PointSegmentDistance(p, s))
	}
}

func TestProjectAndReflectPoint(t *testing.T) {
	l := Line{A: NewPoint2(0, 0), B: NewPoint2(1, 1)}
	// Reflect (1,0) across y=x gives (0,1).
	if !pointApprox(ReflectPointAcrossLine(NewPoint2(1, 0), l), Point2{0, 1}) {
		t.Errorf("ReflectPointAcrossLine = %v", ReflectPointAcrossLine(NewPoint2(1, 0), l))
	}
	if !pointApprox(ProjectPointOntoLine(NewPoint2(2, 0), l), Point2{1, 1}) {
		t.Errorf("ProjectPointOntoLine = %v", ProjectPointOntoLine(NewPoint2(2, 0), l))
	}
}

func TestPerpendicularBisector(t *testing.T) {
	pb := PerpendicularBisector(NewPoint2(0, 0), NewPoint2(4, 0))
	if !pointApprox(pb.A, Point2{2, 0}) {
		t.Errorf("bisector midpoint = %v", pb.A)
	}
	// Direction must be perpendicular to pq (horizontal), so vertical.
	if !approx(pb.Vec().Dot(NewVec2(4, 0)), 0) {
		t.Errorf("bisector not perpendicular")
	}
}

func TestLineIntersection(t *testing.T) {
	l1 := Line{A: NewPoint2(0, 0), B: NewPoint2(1, 1)}
	l2 := Line{A: NewPoint2(0, 2), B: NewPoint2(2, 0)}
	pt, ok := LineIntersection(l1, l2, testTol)
	if !ok || !pointApprox(pt, Point2{1, 1}) {
		t.Errorf("LineIntersection = %v, %v", pt, ok)
	}
	// Parallel lines.
	l3 := Line{A: NewPoint2(0, 1), B: NewPoint2(1, 2)}
	if _, ok := LineIntersection(l1, l3, testTol); ok {
		t.Errorf("parallel lines should not intersect")
	}
}

func TestSegmentIntersection(t *testing.T) {
	// Crossing X.
	s1 := Segment{A: NewPoint2(0, 0), B: NewPoint2(2, 2)}
	s2 := Segment{A: NewPoint2(0, 2), B: NewPoint2(2, 0)}
	pt, ok := SegmentIntersection(s1, s2, testTol)
	if !ok || !pointApprox(pt, Point2{1, 1}) {
		t.Errorf("SegmentIntersection = %v, %v", pt, ok)
	}
	if !SegmentsIntersect(s1, s2, testTol) {
		t.Errorf("SegmentsIntersect should be true")
	}
	// Disjoint.
	s3 := Segment{A: NewPoint2(3, 3), B: NewPoint2(4, 4)}
	if _, ok := SegmentIntersection(s1, s3, testTol); ok {
		t.Errorf("disjoint should not intersect")
	}
	if SegmentsIntersect(s1, s3, testTol) {
		t.Errorf("disjoint SegmentsIntersect false")
	}
	// Shared endpoint (touching).
	s4 := Segment{A: NewPoint2(2, 2), B: NewPoint2(5, 0)}
	if !SegmentsIntersect(s1, s4, testTol) {
		t.Errorf("touching segments should intersect")
	}
	// Collinear overlap: no single point but SegmentsIntersect true.
	s5 := Segment{A: NewPoint2(1, 1), B: NewPoint2(3, 3)}
	if !SegmentsIntersect(s1, s5, testTol) {
		t.Errorf("collinear overlap should intersect")
	}
	if _, ok := SegmentIntersection(s1, s5, testTol); ok {
		t.Errorf("collinear overlap has no unique point")
	}
}

func TestLineSegmentIntersection(t *testing.T) {
	l := Line{A: NewPoint2(0, 1), B: NewPoint2(10, 1)}
	s := Segment{A: NewPoint2(5, 0), B: NewPoint2(5, 5)}
	pt, ok := LineSegmentIntersection(l, s, testTol)
	if !ok || !pointApprox(pt, Point2{5, 1}) {
		t.Errorf("LineSegmentIntersection = %v, %v", pt, ok)
	}
	// Segment entirely below the line.
	s2 := Segment{A: NewPoint2(5, 2), B: NewPoint2(5, 5)}
	if _, ok := LineSegmentIntersection(l, s2, testTol); ok {
		t.Errorf("no intersection expected")
	}
}

func TestSegmentSegmentDistance(t *testing.T) {
	s1 := Segment{A: NewPoint2(0, 0), B: NewPoint2(1, 0)}
	s2 := Segment{A: NewPoint2(0, 3), B: NewPoint2(1, 3)}
	if !approx(SegmentSegmentDistance(s1, s2, testTol), 3) {
		t.Errorf("parallel segment distance = %v", SegmentSegmentDistance(s1, s2, testTol))
	}
	// Intersecting -> 0.
	s3 := Segment{A: NewPoint2(0, 0), B: NewPoint2(2, 2)}
	s4 := Segment{A: NewPoint2(0, 2), B: NewPoint2(2, 0)}
	if !approx(SegmentSegmentDistance(s3, s4, testTol), 0) {
		t.Errorf("intersecting distance = %v", SegmentSegmentDistance(s3, s4, testTol))
	}
}

func TestOnSegmentAndContains(t *testing.T) {
	s := Segment{A: NewPoint2(0, 0), B: NewPoint2(4, 4)}
	if !OnSegment(NewPoint2(2, 2), s, testTol) {
		t.Errorf("midpoint should be on segment")
	}
	if OnSegment(NewPoint2(5, 5), s, testTol) {
		t.Errorf("beyond endpoint not on segment")
	}
	if !s.Contains(NewPoint2(2, 2), testTol) {
		t.Errorf("Contains midpoint")
	}
	l := Line{A: NewPoint2(0, 0), B: NewPoint2(1, 1)}
	if !l.Contains(NewPoint2(5, 5), testTol) {
		t.Errorf("line should contain collinear far point")
	}
}

func TestPolarAndConversions(t *testing.T) {
	p := PointFromPolar(2, math.Pi/2)
	if !pointApprox(p, Point2{0, 2}) {
		t.Errorf("PointFromPolar = %v", p)
	}
	v := VecFromPolar(5, 0)
	if !vecApprox(v, Vec2{5, 0}) {
		t.Errorf("VecFromPolar = %v", v)
	}
	if !vecApprox(NewPoint2(3, 4).Vec(), Vec2{3, 4}) {
		t.Errorf("Point.Vec wrong")
	}
	if !pointApprox(NewVec2(3, 4).AsPoint(), Point2{3, 4}) {
		t.Errorf("Vec.AsPoint wrong")
	}
}

func TestFreeFunctionWrappers(t *testing.T) {
	u, v := NewVec2(1, 2), NewVec2(3, 4)
	if !vecApprox(Add(u, v), u.Add(v)) || !vecApprox(Sub(u, v), u.Sub(v)) {
		t.Errorf("Add/Sub wrappers wrong")
	}
	if !vecApprox(Scale(u, 3), u.Scale(3)) {
		t.Errorf("Scale wrapper wrong")
	}
	if !approx(Dot(u, v), u.Dot(v)) || !approx(Cross(u, v), u.Cross(v)) {
		t.Errorf("Dot/Cross wrappers wrong")
	}
}

func TestSegmentAndLineHelpers(t *testing.T) {
	s := Segment{A: NewPoint2(0, 0), B: NewPoint2(3, 4)}
	if !approx(s.Length(), 5) || !approx(s.LengthSq(), 25) {
		t.Errorf("segment length wrong")
	}
	if !vecApprox(s.Direction(), Vec2{0.6, 0.8}) {
		t.Errorf("segment direction wrong")
	}
	if !pointApprox(s.Midpoint(), Point2{1.5, 2}) {
		t.Errorf("segment midpoint wrong")
	}
	if !pointApprox(s.PointAt(0.5), Point2{1.5, 2}) {
		t.Errorf("segment PointAt wrong")
	}
	l := s.AsLine()
	if !pointApprox(l.PointAt(2), Point2{6, 8}) {
		t.Errorf("line PointAt wrong")
	}
	if !vecApprox(l.Vec(), Vec2{3, 4}) {
		t.Errorf("line Vec wrong")
	}
}

func TestEqualityPredicates(t *testing.T) {
	if !NewPoint2(1, 2).Equal(NewPoint2(1, 2)) {
		t.Errorf("Point Equal wrong")
	}
	if !NewPoint2(1, 2).ApproxEqual(NewPoint2(1+1e-12, 2), testTol) {
		t.Errorf("Point ApproxEqual wrong")
	}
	if !NewVec2(1, 2).Equal(NewVec2(1, 2)) {
		t.Errorf("Vec Equal wrong")
	}
	if !NewVec2(1, 2).ApproxEqual(NewVec2(1, 2+1e-12), testTol) {
		t.Errorf("Vec ApproxEqual wrong")
	}
}

func BenchmarkSegmentIntersection(b *testing.B) {
	s1 := Segment{A: NewPoint2(0, 0), B: NewPoint2(2, 2)}
	s2 := Segment{A: NewPoint2(0, 2), B: NewPoint2(2, 0)}
	var acc float64
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		pt, ok := SegmentIntersection(s1, s2, geom2dEps)
		if ok {
			acc += pt.X + pt.Y
		}
	}
	_ = acc
}

func BenchmarkSegmentSegmentDistance(b *testing.B) {
	s1 := Segment{A: NewPoint2(0, 0), B: NewPoint2(1, 0)}
	s2 := Segment{A: NewPoint2(0, 3), B: NewPoint2(1, 3)}
	var acc float64
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		acc += SegmentSegmentDistance(s1, s2, geom2dEps)
	}
	_ = acc
}
