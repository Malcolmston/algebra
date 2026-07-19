package hull3d

import (
	"errors"
	"math"
)

// ConvexShape is any convex body that can report its support point: the point
// of the body farthest in a given direction. Support functions are the only
// interface the GJK and EPA algorithms need, so arbitrary convex shapes
// (polytopes, spheres, capsules, sums) can be handled uniformly.
type ConvexShape interface {
	// Support returns the point of the shape maximising the dot product with d.
	Support(d Vec3) Vec3
}

// PointShape is a single point treated as a degenerate convex body.
type PointShape struct{ P Vec3 }

// Support implements [ConvexShape].
func (s PointShape) Support(Vec3) Vec3 { return s.P }

// SegmentShape is the line segment between A and B.
type SegmentShape struct{ A, B Vec3 }

// Support implements [ConvexShape].
func (s SegmentShape) Support(d Vec3) Vec3 {
	if s.A.Dot(d) >= s.B.Dot(d) {
		return s.A
	}
	return s.B
}

// SphereShape is a ball with the given centre and radius.
type SphereShape struct {
	Center Vec3
	Radius float64
}

// Support implements [ConvexShape].
func (s SphereShape) Support(d Vec3) Vec3 {
	u, err := d.Normalize()
	if err != nil {
		return s.Center
	}
	return s.Center.Add(u.Scale(s.Radius))
}

// BoxShape is an axis-aligned box with the given centre and half-extents.
type BoxShape struct {
	Center Vec3
	Half   Vec3
}

// Support implements [ConvexShape].
func (s BoxShape) Support(d Vec3) Vec3 {
	sgn := func(x, h float64) float64 {
		if x >= 0 {
			return h
		}
		return -h
	}
	return s.Center.Add(Vec3{sgn(d.X, s.Half.X), sgn(d.Y, s.Half.Y), sgn(d.Z, s.Half.Z)})
}

// PointCloudShape is the convex hull of a finite set of points, given only by
// the points themselves (no explicit hull needed for support queries).
type PointCloudShape struct{ Points []Vec3 }

// Support implements [ConvexShape].
func (s PointCloudShape) Support(d Vec3) Vec3 {
	_, v, err := FarthestPoint(s.Points, d)
	if err != nil {
		return Vec3{}
	}
	return v
}

// PolytopeShape adapts a [Polytope] to [ConvexShape] using its vertices.
type PolytopeShape struct{ P *Polytope }

// Support implements [ConvexShape].
func (s PolytopeShape) Support(d Vec3) Vec3 {
	v, err := s.P.Support(d)
	if err != nil {
		return Vec3{}
	}
	return v
}

// TranslatedShape is a convex shape shifted by Offset.
type TranslatedShape struct {
	Shape  ConvexShape
	Offset Vec3
}

// Support implements [ConvexShape].
func (s TranslatedShape) Support(d Vec3) Vec3 {
	return s.Shape.Support(d).Add(s.Offset)
}

// MinkowskiSumShape is the Minkowski sum of two convex shapes, whose support in
// direction d is the sum of the two supports in d.
type MinkowskiSumShape struct{ A, B ConvexShape }

// Support implements [ConvexShape].
func (s MinkowskiSumShape) Support(d Vec3) Vec3 {
	return s.A.Support(d).Add(s.B.Support(d))
}

// mvert is a vertex of the Minkowski difference A-B, retaining the individual
// supports so witness points on A and B can be recovered by interpolation.
type mvert struct {
	p, sa, sb Vec3
}

func mdSupport(a, b ConvexShape, d Vec3) mvert {
	sa := a.Support(d)
	sb := b.Support(d.Neg())
	return mvert{p: sa.Sub(sb), sa: sa, sb: sb}
}

// GJKIntersect reports whether the two convex shapes overlap, using the
// Gilbert–Johnson–Keerthi algorithm to test whether the origin lies in their
// Minkowski difference.
func GJKIntersect(a, b ConvexShape) bool {
	ok, _ := gjkSimplex(a, b)
	return ok
}

// gjkSimplex runs GJK to determine intersection. On intersection it returns a
// tetrahedron simplex (4 mverts) enclosing the origin, suitable for EPA.
func gjkSimplex(a, b ConvexShape) (bool, []mvert) {
	d := Vec3{1, 0, 0}
	s := mdSupport(a, b, d)
	simplex := []mvert{s}
	d = s.p.Neg()
	if d.IsZero() {
		d = Vec3{1, 0, 0}
	}
	for iter := 0; iter < 64; iter++ {
		if d.LengthSq() < 1e-30 {
			return true, ensureTetra(a, b, simplex)
		}
		w := mdSupport(a, b, d)
		if w.p.Dot(d) < 0 {
			return false, nil
		}
		simplex = append(simplex, w)
		var contains bool
		simplex, d, contains = doSimplex(simplex)
		if contains {
			return true, ensureTetra(a, b, simplex)
		}
	}
	return false, nil
}

// ensureTetra returns a tetrahedron of Minkowski-difference support points that
// strictly encloses the origin. It uses the GJK terminal simplex when that
// already qualifies and otherwise rebuilds one with [expandToTetra].
func ensureTetra(a, b ConvexShape, simplex []mvert) []mvert {
	if len(simplex) == 4 {
		vol := Orient3D(simplex[0].p, simplex[1].p, simplex[2].p, simplex[3].p)
		if math.Abs(vol) > 1e-14 &&
			pointStrictlyInTetra(Vec3{}, simplex[0].p, simplex[1].p, simplex[2].p, simplex[3].p, vol) {
			return simplex
		}
	}
	return expandToTetra(a, b, simplex)
}

// expandToTetra builds a non-degenerate tetrahedron of Minkowski-difference
// support points that strictly encloses the origin, suitable for seeding EPA. It
// collects a spread of support points (the incoming simplex, the six axis
// extremes and several diagonal directions) and searches four of them that
// strictly enclose the origin with the largest volume. When the two bodies
// overlap the origin is strictly interior to their Minkowski difference, so such
// a quadruple exists.
func expandToTetra(a, b ConvexShape, simplex []mvert) []mvert {
	dirs := []Vec3{
		{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
		{1, 1, 1}, {-1, -1, -1}, {1, -1, 1}, {-1, 1, -1},
		{1, 1, -1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, 1},
	}
	var cand []mvert
	add := func(w mvert) {
		for _, v := range cand {
			if v.p.ApproxEqual(w.p, 1e-12) {
				return
			}
		}
		cand = append(cand, w)
	}
	for _, w := range simplex {
		add(w)
	}
	for _, d := range dirs {
		add(mdSupport(a, b, d))
	}

	best := -1.0
	var bestTet []mvert
	n := len(cand)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			for k := j + 1; k < n; k++ {
				for l := k + 1; l < n; l++ {
					vol := Orient3D(cand[i].p, cand[j].p, cand[k].p, cand[l].p)
					av := math.Abs(vol)
					if av <= best {
						continue
					}
					if pointStrictlyInTetra(Vec3{}, cand[i].p, cand[j].p, cand[k].p, cand[l].p, vol) {
						best = av
						bestTet = []mvert{cand[i], cand[j], cand[k], cand[l]}
					}
				}
			}
		}
	}
	if bestTet != nil {
		return bestTet
	}
	return expandToTetraSearch(a, b, simplex)
}

// expandToTetraSearch is the fallback incremental builder used when the
// candidate search fails to find a strictly enclosing tetrahedron.
func expandToTetraSearch(a, b ConvexShape, _ []mvert) []mvert {
	// 1. Any first support, then search back toward the origin for the second.
	s0 := mdSupport(a, b, Vec3{0, 1, 0})
	dir := s0.p.Neg()
	if dir.LengthSq() < 1e-20 {
		dir = Vec3{1, 0, 0}
	}
	s1 := mdSupport(a, b, dir)

	// 2. Perpendicular to segment s0-s1, oriented toward the origin.
	ab := s1.p.Sub(s0.p)
	dir = ab.Cross(s0.p.Neg()).Cross(ab)
	if dir.LengthSq() < 1e-20 {
		dir = ab.Cross(Vec3{1, 0, 0})
		if dir.LengthSq() < 1e-20 {
			dir = ab.Cross(Vec3{0, 0, 1})
		}
	}
	s2 := mdSupport(a, b, dir)

	// 3. Triangle normal, oriented toward the origin, gives the apex.
	n := s1.p.Sub(s0.p).Cross(s2.p.Sub(s0.p))
	if n.Dot(s0.p) > 0 {
		n = n.Neg()
	}
	s3 := mdSupport(a, b, n)

	verts := []mvert{s0, s1, s2, s3}
	// Ensure the origin is strictly enclosed. Each iteration finds the face
	// relative to which the origin is most outside (or on) and pushes the
	// opposite vertex past it with a fresh support query.
	faces := [4][3]int{{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2}}
	opp := [4]int{3, 1, 2, 0}
	for iter := 0; iter < 16; iter++ {
		worst := -1
		worstVal := math.Inf(-1)
		var worstDir Vec3
		for fi, f := range faces {
			fa, fb, fc := verts[f[0]].p, verts[f[1]].p, verts[f[2]].p
			nn := fb.Sub(fa).Cross(fc.Sub(fa))
			if l := nn.Length(); l > 1e-300 {
				nn = nn.Scale(1 / l)
			}
			// Orient outward (away from the opposite vertex).
			if nn.Dot(verts[opp[fi]].p.Sub(fa)) > 0 {
				nn = nn.Neg()
			}
			val := nn.Dot(fa.Neg()) // signed distance of origin outside this face
			if val > worstVal {
				worstVal, worst, worstDir = val, fi, nn
			}
		}
		if worstVal < -1e-12 {
			break // origin strictly inside
		}
		w := mdSupport(a, b, worstDir)
		if w.p.Dot(worstDir)-verts[opp[worst]].p.Dot(worstDir) < 1e-12 {
			break // cannot push the face any farther; origin on boundary
		}
		verts[opp[worst]] = w
	}
	return verts
}

// doSimplex processes the current GJK simplex, reduces it toward the origin, and
// returns the new simplex, the next search direction, and whether the origin is
// enclosed.
func doSimplex(s []mvert) ([]mvert, Vec3, bool) {
	switch len(s) {
	case 2:
		return doLine(s)
	case 3:
		return doTriangle(s)
	case 4:
		return doTetra(s)
	}
	return s, s[len(s)-1].p.Neg(), false
}

func doLine(s []mvert) ([]mvert, Vec3, bool) {
	a := s[1]
	b := s[0]
	ab := b.p.Sub(a.p)
	ao := a.p.Neg()
	if ab.Dot(ao) > 0 {
		d := ab.Cross(ao).Cross(ab)
		if d.LengthSq() < 1e-30 {
			return []mvert{a, b}, arbitraryPerp(ab), false
		}
		return []mvert{b, a}, d, false
	}
	return []mvert{a}, ao, false
}

func doTriangle(s []mvert) ([]mvert, Vec3, bool) {
	a := s[2]
	b := s[1]
	c := s[0]
	ab := b.p.Sub(a.p)
	ac := c.p.Sub(a.p)
	ao := a.p.Neg()
	abc := ab.Cross(ac)

	if abc.Cross(ac).Dot(ao) > 0 {
		if ac.Dot(ao) > 0 {
			return []mvert{c, a}, ac.Cross(ao).Cross(ac), false
		}
		return doLine([]mvert{b, a})
	}
	if ab.Cross(abc).Dot(ao) > 0 {
		return doLine([]mvert{b, a})
	}
	if abc.Dot(ao) > 0 {
		return []mvert{c, b, a}, abc, false
	}
	return []mvert{b, c, a}, abc.Neg(), false
}

func doTetra(s []mvert) ([]mvert, Vec3, bool) {
	a := s[3]
	b := s[2]
	c := s[1]
	d := s[0]
	ao := a.p.Neg()
	abc := b.p.Sub(a.p).Cross(c.p.Sub(a.p))
	acd := c.p.Sub(a.p).Cross(d.p.Sub(a.p))
	adb := d.p.Sub(a.p).Cross(b.p.Sub(a.p))

	// Orient so each normal points away from the opposite vertex.
	if abc.Dot(d.p.Sub(a.p)) > 0 {
		abc = abc.Neg()
	}
	if acd.Dot(b.p.Sub(a.p)) > 0 {
		acd = acd.Neg()
	}
	if adb.Dot(c.p.Sub(a.p)) > 0 {
		adb = adb.Neg()
	}

	if abc.Dot(ao) > 0 {
		return doTriangle([]mvert{c, b, a})
	}
	if acd.Dot(ao) > 0 {
		return doTriangle([]mvert{d, c, a})
	}
	if adb.Dot(ao) > 0 {
		return doTriangle([]mvert{b, d, a})
	}
	return []mvert{d, c, b, a}, Vec3{}, true
}

func arbitraryPerp(v Vec3) Vec3 {
	p, err := v.AnyPerpendicular()
	if err != nil {
		return Vec3{1, 0, 0}
	}
	return p
}

// GJKDistance returns the Euclidean distance between two convex shapes together
// with the closest points on each. When the shapes overlap the distance is zero
// and the returned witness points coincide at an arbitrary contact point; use
// [EPAPenetration] to measure penetration depth in that case.
func GJKDistance(a, b ConvexShape) (dist float64, ptA, ptB Vec3, err error) {
	if a == nil || b == nil {
		return 0, Vec3{}, Vec3{}, errors.New("hull3d: nil shape")
	}
	d := Vec3{1, 0, 0}
	simplex := []mvert{mdSupport(a, b, d)}
	d = simplex[0].p.Neg()
	prev := math.Inf(1)
	for iter := 0; iter < 128; iter++ {
		if d.LengthSq() < 1e-30 {
			// Origin in Minkowski difference: shapes overlap.
			return 0, simplex[0].sa, simplex[0].sa, nil
		}
		w := mdSupport(a, b, d)
		duplicate := false
		for _, v := range simplex {
			if v.p.ApproxEqual(w.p, 1e-15) {
				duplicate = true
				break
			}
		}
		simplex = append(simplex, w)
		closest, weights, reduced := closestOnSimplex(simplex)
		simplex = reduced
		dist = closest.Length()
		if closest.LengthSq() < 1e-30 {
			return 0, simplex[0].sa, simplex[0].sa, nil
		}
		if duplicate || math.Abs(prev-dist) < 1e-12*(1+dist) {
			pa, pb := witnessPoints(simplex, weights)
			return dist, pa, pb, nil
		}
		prev = dist
		d = closest.Neg()
	}
	closest, weights, _ := closestOnSimplex(simplex)
	pa, pb := witnessPoints(simplex, weights)
	return closest.Length(), pa, pb, nil
}

func witnessPoints(s []mvert, w []float64) (Vec3, Vec3) {
	var pa, pb Vec3
	for i := range s {
		pa = pa.Add(s[i].sa.Scale(w[i]))
		pb = pb.Add(s[i].sb.Scale(w[i]))
	}
	return pa, pb
}

// closestOnSimplex returns the point of the simplex closest to the origin, the
// barycentric weights of that point over the returned reduced simplex, and the
// reduced simplex (only the vertices with non-zero weight).
func closestOnSimplex(s []mvert) (Vec3, []float64, []mvert) {
	switch len(s) {
	case 1:
		return s[0].p, []float64{1}, s
	case 2:
		return closestSeg(s)
	case 3:
		return closestTri(s)
	default:
		return closestTetra(s)
	}
}

func closestSeg(s []mvert) (Vec3, []float64, []mvert) {
	a, b := s[0].p, s[1].p
	ab := b.Sub(a)
	t := a.Neg().Dot(ab)
	den := ab.LengthSq()
	if den < 1e-300 || t <= 0 {
		return a, []float64{1}, s[:1]
	}
	if t >= den {
		return b, []float64{1}, s[1:]
	}
	t /= den
	return a.Add(ab.Scale(t)), []float64{1 - t, t}, s
}

func closestTri(s []mvert) (Vec3, []float64, []mvert) {
	a, b, c := s[0].p, s[1].p, s[2].p
	// Ericson closest point on triangle to origin (P = 0).
	ab := b.Sub(a)
	ac := c.Sub(a)
	ap := a.Neg()
	d1 := ab.Dot(ap)
	d2 := ac.Dot(ap)
	if d1 <= 0 && d2 <= 0 {
		return a, []float64{1}, s[:1]
	}
	bp := b.Neg()
	d3 := ab.Dot(bp)
	d4 := ac.Dot(bp)
	if d3 >= 0 && d4 <= d3 {
		return b, []float64{1}, []mvert{s[1]}
	}
	vc := d1*d4 - d3*d2
	if vc <= 0 && d1 >= 0 && d3 <= 0 {
		v := d1 / (d1 - d3)
		return a.Add(ab.Scale(v)), []float64{1 - v, v}, []mvert{s[0], s[1]}
	}
	cp := c.Neg()
	d5 := ab.Dot(cp)
	d6 := ac.Dot(cp)
	if d6 >= 0 && d5 <= d6 {
		return c, []float64{1}, []mvert{s[2]}
	}
	vb := d5*d2 - d1*d6
	if vb <= 0 && d2 >= 0 && d6 <= 0 {
		w := d2 / (d2 - d6)
		return a.Add(ac.Scale(w)), []float64{1 - w, w}, []mvert{s[0], s[2]}
	}
	va := d3*d6 - d5*d4
	if va <= 0 && (d4-d3) >= 0 && (d5-d6) >= 0 {
		w := (d4 - d3) / ((d4 - d3) + (d5 - d6))
		return b.Add(c.Sub(b).Scale(w)), []float64{1 - w, w}, []mvert{s[1], s[2]}
	}
	denom := 1 / (va + vb + vc)
	v := vb * denom
	w := vc * denom
	pt := a.Add(ab.Scale(v)).Add(ac.Scale(w))
	return pt, []float64{1 - v - w, v, w}, s
}

func closestTetra(s []mvert) (Vec3, []float64, []mvert) {
	// If the origin is strictly inside a non-degenerate tetrahedron, distance is
	// zero. The volume gate avoids treating a flat simplex as enclosing.
	a, b, c, d := s[0].p, s[1].p, s[2].p, s[3].p
	vol := Orient3D(a, b, c, d)
	if math.Abs(vol) > 1e-14 && pointStrictlyInTetra(Vec3{}, a, b, c, d, vol) {
		return Vec3{}, []float64{0.25, 0.25, 0.25, 0.25}, s
	}
	bestD := math.Inf(1)
	var bestP Vec3
	var bestW []float64
	var bestS []mvert
	tris := [][3]int{{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2}}
	for _, t := range tris {
		sub := []mvert{s[t[0]], s[t[1]], s[t[2]]}
		p, w, rs := closestTri(sub)
		if l := p.LengthSq(); l < bestD {
			bestD, bestP, bestW, bestS = l, p, w, rs
		}
	}
	return bestP, bestW, bestS
}

// pointStrictlyInTetra reports whether p is strictly inside the tetrahedron
// (a,b,c,d) whose signed volume determinant is vol (assumed non-zero). A point on
// any face (zero sub-determinant) is not strictly inside.
func pointStrictlyInTetra(p, a, b, c, d Vec3, vol float64) bool {
	s := math.Copysign(1, vol)
	o0 := Orient3D(p, b, c, d)
	o1 := Orient3D(a, p, c, d)
	o2 := Orient3D(a, b, p, d)
	o3 := Orient3D(a, b, c, p)
	return o0 != 0 && math.Copysign(1, o0) == s &&
		o1 != 0 && math.Copysign(1, o1) == s &&
		o2 != 0 && math.Copysign(1, o2) == s &&
		o3 != 0 && math.Copysign(1, o3) == s
}
