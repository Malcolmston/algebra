package hull3d

import (
	"errors"
	"math"
)

// SupportingPlaneOf returns the supporting plane of the polytope in direction d:
// the plane with outward normal d touching the polytope at its farthest vertex,
// so the whole polytope lies in the plane's negative half-space. It returns an
// error if the polytope has no vertices.
func (p *Polytope) SupportingPlaneOf(d Vec3) (Plane, error) {
	v, err := p.Support(d)
	if err != nil {
		return Plane{}, err
	}
	return Plane{Normal: d, Offset: d.Dot(v)}, nil
}

// Width returns the extent of the polytope along direction d: the difference
// between its maximum and minimum support in d. It returns an error for an empty
// polytope.
func (p *Polytope) Width(d Vec3) (float64, error) {
	if len(p.Vertices) == 0 {
		return 0, errors.New("hull3d: empty polytope")
	}
	hi := math.Inf(-1)
	lo := math.Inf(1)
	for _, v := range p.Vertices {
		x := v.Dot(d)
		hi = math.Max(hi, x)
		lo = math.Min(lo, x)
	}
	nd := d.Length()
	if nd < 1e-300 {
		return 0, errors.New("hull3d: zero direction")
	}
	return (hi - lo) / nd, nil
}

// SeparatingPlane returns a plane strictly separating two disjoint convex
// polytopes, with a and b on the negative and positive sides respectively. It
// uses the closest points reported by GJK to build the plane. It returns an
// error if the polytopes overlap or touch (no strict separator exists).
func SeparatingPlane(a, b *Polytope, eps float64) (Plane, error) {
	if a == nil || b == nil {
		return Plane{}, errors.New("hull3d: nil polytope")
	}
	dist, pa, pb, err := GJKDistance(PolytopeShape{a}, PolytopeShape{b})
	if err != nil {
		return Plane{}, err
	}
	if dist <= eps {
		return Plane{}, errors.New("hull3d: polytopes are not strictly separated")
	}
	normal, err := pb.Sub(pa).Normalize()
	if err != nil {
		return Plane{}, err
	}
	mid := pa.Midpoint(pb)
	// Oriented so a (near pa) is on the negative side.
	return Plane{Normal: normal, Offset: normal.Dot(mid)}, nil
}

// AreDisjoint reports whether two convex polytopes are disjoint (do not overlap)
// within tolerance eps, via their GJK distance.
func AreDisjoint(a, b *Polytope, eps float64) bool {
	dist, _, _, err := GJKDistance(PolytopeShape{a}, PolytopeShape{b})
	if err != nil {
		return false
	}
	return dist > eps
}

// ClosestPointInPolytope returns the point of the (convex) polytope closest to
// q and the distance to it, using GJK between the polytope and the single point
// q. It returns an error for an empty polytope.
func ClosestPointInPolytope(p *Polytope, q Vec3) (Vec3, float64, error) {
	if p == nil || len(p.Vertices) == 0 {
		return Vec3{}, 0, errors.New("hull3d: empty polytope")
	}
	dist, pt, _, err := GJKDistance(PolytopeShape{p}, PointShape{q})
	if err != nil {
		return Vec3{}, 0, err
	}
	if dist == 0 {
		return q, 0, nil
	}
	return pt, dist, nil
}

// DistancePointToPolytope returns the Euclidean distance from q to the polytope
// (zero if q is inside). It returns an error for an empty polytope.
func DistancePointToPolytope(p *Polytope, q Vec3) (float64, error) {
	_, d, err := ClosestPointInPolytope(p, q)
	return d, err
}

// RayPolytopeIntersection computes the entry and exit parameters of the ray
// origin + t*dir with the convex polytope, using the slab method against every
// face half-space. It returns (tEnter, tExit, true) with 0 <= tEnter <= tExit
// when the ray meets the polytope, or ok=false when it misses. A ray starting
// inside yields tEnter = 0.
func RayPolytopeIntersection(p *Polytope, origin, dir Vec3) (tEnter, tExit float64, ok bool) {
	tEnter = 0
	tExit = math.Inf(1)
	for _, f := range p.Faces {
		pl, err := p.FacePlane(f)
		if err != nil {
			continue
		}
		// Face plane: n·x <= offset is the inside.
		denom := pl.Normal.Dot(dir)
		num := pl.Offset - pl.Normal.Dot(origin) // >= 0 inside currently
		if math.Abs(denom) < 1e-300 {
			if num < 0 {
				return 0, 0, false // parallel and outside this slab
			}
			continue
		}
		t := num / denom
		if denom < 0 {
			// entering this half-space as t increases
			if t > tEnter {
				tEnter = t
			}
		} else {
			if t < tExit {
				tExit = t
			}
		}
		if tEnter > tExit {
			return 0, 0, false
		}
	}
	if tExit < 0 {
		return 0, 0, false
	}
	return tEnter, tExit, true
}

// SegmentPolytopeIntersect reports whether the segment from s0 to s1 meets the
// convex polytope.
func SegmentPolytopeIntersect(p *Polytope, s0, s1 Vec3) bool {
	dir := s1.Sub(s0)
	tEnter, _, ok := RayPolytopeIntersection(p, s0, dir)
	if !ok {
		return false
	}
	return tEnter <= 1
}

// LineSegmentClosest returns the closest points between two segments (p0,p1) and
// (q0,q1) and the distance between them.
func LineSegmentClosest(p0, p1, q0, q1 Vec3) (cp, cq Vec3, dist float64) {
	d1 := p1.Sub(p0)
	d2 := q1.Sub(q0)
	r := p0.Sub(q0)
	a := d1.Dot(d1)
	e := d2.Dot(d2)
	f := d2.Dot(r)
	var s, t float64
	const eps = 1e-300
	if a <= eps && e <= eps {
		return p0, q0, p0.Distance(q0)
	}
	if a <= eps {
		t = clamp01(f / e)
	} else {
		c := d1.Dot(r)
		if e <= eps {
			s = clamp01(-c / a)
		} else {
			b := d1.Dot(d2)
			den := a*e - b*b
			if den > eps {
				s = clamp01((b*f - c*e) / den)
			}
			t = (b*s + f) / e
			if t < 0 {
				t = 0
				s = clamp01(-c / a)
			} else if t > 1 {
				t = 1
				s = clamp01((b - c) / a)
			}
		}
	}
	cp = p0.Add(d1.Scale(s))
	cq = q0.Add(d2.Scale(t))
	return cp, cq, cp.Distance(cq)
}

func clamp01(x float64) float64 {
	if x < 0 {
		return 0
	}
	if x > 1 {
		return 1
	}
	return x
}

// ClosestPointOnTriangle returns the point on triangle (a,b,c) closest to q.
func ClosestPointOnTriangle(q, a, b, c Vec3) Vec3 {
	ab := b.Sub(a)
	ac := c.Sub(a)
	ap := q.Sub(a)
	d1 := ab.Dot(ap)
	d2 := ac.Dot(ap)
	if d1 <= 0 && d2 <= 0 {
		return a
	}
	bp := q.Sub(b)
	d3 := ab.Dot(bp)
	d4 := ac.Dot(bp)
	if d3 >= 0 && d4 <= d3 {
		return b
	}
	vc := d1*d4 - d3*d2
	if vc <= 0 && d1 >= 0 && d3 <= 0 {
		v := d1 / (d1 - d3)
		return a.Add(ab.Scale(v))
	}
	cp := q.Sub(c)
	d5 := ab.Dot(cp)
	d6 := ac.Dot(cp)
	if d6 >= 0 && d5 <= d6 {
		return c
	}
	vb := d5*d2 - d1*d6
	if vb <= 0 && d2 >= 0 && d6 <= 0 {
		w := d2 / (d2 - d6)
		return a.Add(ac.Scale(w))
	}
	va := d3*d6 - d5*d4
	if va <= 0 && (d4-d3) >= 0 && (d5-d6) >= 0 {
		w := (d4 - d3) / ((d4 - d3) + (d5 - d6))
		return b.Add(c.Sub(b).Scale(w))
	}
	denom := 1 / (va + vb + vc)
	v := vb * denom
	w := vc * denom
	return a.Add(ab.Scale(v)).Add(ac.Scale(w))
}

// PointInTetrahedron reports whether q lies inside (or on) the tetrahedron with
// the given corners, within tolerance eps on the orientation determinant.
func PointInTetrahedron(q, a, b, c, d Vec3, eps float64) bool {
	ref := Orient3D(a, b, c, d)
	if math.Abs(ref) < eps {
		return false
	}
	s := math.Copysign(1, ref)
	tests := [][4]Vec3{
		{q, b, c, d},
		{a, q, c, d},
		{a, b, q, d},
		{a, b, c, q},
	}
	for _, t := range tests {
		o := Orient3D(t[0], t[1], t[2], t[3])
		if math.Copysign(1, o) != s && math.Abs(o) > eps {
			return false
		}
	}
	return true
}

// BarycentricTetrahedron returns the four barycentric coordinates of q with
// respect to the tetrahedron (a,b,c,d). They sum to one; q is inside iff all are
// non-negative. It returns an error for a degenerate (flat) tetrahedron.
func BarycentricTetrahedron(q, a, b, c, d Vec3) ([4]float64, error) {
	vol := Orient3D(a, b, c, d)
	if math.Abs(vol) < 1e-300 {
		return [4]float64{}, errors.New("hull3d: degenerate tetrahedron")
	}
	wa := Orient3D(q, b, c, d) / vol
	wb := Orient3D(a, q, c, d) / vol
	wc := Orient3D(a, b, q, d) / vol
	wd := Orient3D(a, b, c, q) / vol
	return [4]float64{wa, wb, wc, wd}, nil
}
