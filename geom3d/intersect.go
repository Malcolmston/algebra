package geom3d

import "math"

// Barycentric returns the barycentric coordinates (u, v, w) of point p with
// respect to triangle (a, b, c), so that p ≈ u*a + v*b + w*c and u+v+w = 1. The
// point is projected onto the triangle's plane if it does not lie exactly on
// it. It returns false if the triangle is degenerate (its vertices are
// collinear).
func Barycentric(p, a, b, c Vec3) (u, v, w float64, ok bool) {
	v0 := b.Sub(a)
	v1 := c.Sub(a)
	v2 := p.Sub(a)
	d00 := v0.Dot(v0)
	d01 := v0.Dot(v1)
	d11 := v1.Dot(v1)
	d20 := v2.Dot(v0)
	d21 := v2.Dot(v1)
	denom := d00*d11 - d01*d01
	if math.Abs(denom) <= geom3dEps {
		return 0, 0, 0, false
	}
	v = (d11*d20 - d01*d21) / denom
	w = (d00*d21 - d01*d20) / denom
	u = 1 - v - w
	return u, v, w, true
}

// PointInTriangle reports whether point p lies inside (or on the edge of)
// triangle (a, b, c). The point is assumed to lie in the triangle's plane;
// callers testing an arbitrary point should first project it. A degenerate
// triangle yields false.
func PointInTriangle(p, a, b, c Vec3) bool {
	u, v, w, ok := Barycentric(p, a, b, c)
	if !ok {
		return false
	}
	return u >= -geom3dEps && v >= -geom3dEps && w >= -geom3dEps
}

// RayTriangle intersects a ray with triangle (a, b, c) using the
// Moeller-Trumbore algorithm. On a hit it returns the ray parameter t (distance
// along Dir to the intersection), the barycentric coordinates u and v within
// the triangle (with the a-weight equal to 1-u-v), and true. It returns false
// if the ray is parallel to the triangle, misses it, or would intersect behind
// the origin. Both faces of the triangle are considered.
func RayTriangle(r Ray, a, b, c Vec3) (t, u, v float64, ok bool) {
	edge1 := b.Sub(a)
	edge2 := c.Sub(a)
	pvec := r.Dir.Cross(edge2)
	det := edge1.Dot(pvec)
	if math.Abs(det) <= geom3dEps {
		return 0, 0, 0, false
	}
	invDet := 1 / det
	tvec := r.Origin.Sub(a)
	u = tvec.Dot(pvec) * invDet
	if u < -geom3dEps || u > 1+geom3dEps {
		return 0, 0, 0, false
	}
	qvec := tvec.Cross(edge1)
	v = r.Dir.Dot(qvec) * invDet
	if v < -geom3dEps || u+v > 1+geom3dEps {
		return 0, 0, 0, false
	}
	t = edge2.Dot(qvec) * invDet
	if t < 0 {
		return 0, 0, 0, false
	}
	return t, u, v, true
}

// AABB is an axis-aligned bounding box defined by its minimum and maximum
// corners. A valid box has Min[i] <= Max[i] for every axis i.
type AABB struct {
	Min, Max Vec3
}

// NewAABB returns the axis-aligned box with the given corners, normalizing them
// so that Min holds the component-wise minimum and Max the component-wise
// maximum.
func NewAABB(a, b Vec3) AABB {
	return AABB{Min: a.Min(b), Max: a.Max(b)}
}

// AABBFromPoints returns the smallest axis-aligned box containing all the given
// points, and true. It returns the zero box and false if no points are given.
func AABBFromPoints(pts ...Vec3) (AABB, bool) {
	if len(pts) == 0 {
		return AABB{}, false
	}
	mn, mx := pts[0], pts[0]
	for _, p := range pts[1:] {
		mn = mn.Min(p)
		mx = mx.Max(p)
	}
	return AABB{Min: mn, Max: mx}, true
}

// Center returns the geometric center of the box.
func (bx AABB) Center() Vec3 {
	return bx.Min.Add(bx.Max).Scale(0.5)
}

// Size returns the full extent of the box along each axis (Max - Min).
func (bx AABB) Size() Vec3 {
	return bx.Max.Sub(bx.Min)
}

// Extents returns the half-size of the box along each axis.
func (bx AABB) Extents() Vec3 {
	return bx.Size().Scale(0.5)
}

// Contains reports whether point p lies inside or on the boundary of the box.
func (bx AABB) Contains(p Vec3) bool {
	return p.X >= bx.Min.X && p.X <= bx.Max.X &&
		p.Y >= bx.Min.Y && p.Y <= bx.Max.Y &&
		p.Z >= bx.Min.Z && p.Z <= bx.Max.Z
}

// ClosestPoint returns the point within the box (on its surface or interior)
// nearest to p, obtained by clamping p to the box on each axis.
func (bx AABB) ClosestPoint(p Vec3) Vec3 {
	return Vec3{
		geom3dclamp(p.X, bx.Min.X, bx.Max.X),
		geom3dclamp(p.Y, bx.Min.Y, bx.Max.Y),
		geom3dclamp(p.Z, bx.Min.Z, bx.Max.Z),
	}
}

// IntersectsAABB reports whether the box overlaps the other box, touching
// boundaries included.
func (bx AABB) IntersectsAABB(other AABB) bool {
	return bx.Min.X <= other.Max.X && bx.Max.X >= other.Min.X &&
		bx.Min.Y <= other.Max.Y && bx.Max.Y >= other.Min.Y &&
		bx.Min.Z <= other.Max.Z && bx.Max.Z >= other.Min.Z
}

// Union returns the smallest axis-aligned box containing both bx and other.
func (bx AABB) Union(other AABB) AABB {
	return AABB{Min: bx.Min.Min(other.Min), Max: bx.Max.Max(other.Max)}
}

// Expand returns the box grown outward by the vector v on every side; each
// component of v is subtracted from Min and added to Max.
func (bx AABB) Expand(v Vec3) AABB {
	return AABB{Min: bx.Min.Sub(v), Max: bx.Max.Add(v)}
}

// IntersectsSphere reports whether the box overlaps the sphere s.
func (bx AABB) IntersectsSphere(s Sphere) bool {
	return s.IntersectsAABB(bx)
}

// IntersectRay intersects a ray with the box using the slab method. On a hit it
// returns the entry parameter tmin (distance along the ray to the first
// intersection, clamped to 0 if the origin is inside), the exit parameter tmax,
// and true. It returns false if the ray misses the box or the box lies entirely
// behind the origin.
func (bx AABB) IntersectRay(r Ray) (tmin, tmax float64, ok bool) {
	tmin = math.Inf(-1)
	tmax = math.Inf(1)
	origin := [3]float64{r.Origin.X, r.Origin.Y, r.Origin.Z}
	dir := [3]float64{r.Dir.X, r.Dir.Y, r.Dir.Z}
	mn := [3]float64{bx.Min.X, bx.Min.Y, bx.Min.Z}
	mx := [3]float64{bx.Max.X, bx.Max.Y, bx.Max.Z}
	for i := 0; i < 3; i++ {
		if math.Abs(dir[i]) <= geom3dEps {
			// Ray parallel to this slab: origin must lie within it.
			if origin[i] < mn[i] || origin[i] > mx[i] {
				return 0, 0, false
			}
			continue
		}
		inv := 1 / dir[i]
		t1 := (mn[i] - origin[i]) * inv
		t2 := (mx[i] - origin[i]) * inv
		if t1 > t2 {
			t1, t2 = t2, t1
		}
		if t1 > tmin {
			tmin = t1
		}
		if t2 < tmax {
			tmax = t2
		}
		if tmin > tmax {
			return 0, 0, false
		}
	}
	if tmax < 0 {
		return 0, 0, false
	}
	if tmin < 0 {
		tmin = 0
	}
	return tmin, tmax, true
}

// Sphere is a solid ball defined by its Center and Radius.
type Sphere struct {
	Center Vec3
	Radius float64
}

// Contains reports whether point p lies inside or on the surface of the sphere.
func (s Sphere) Contains(p Vec3) bool {
	return s.Center.DistanceSq(p) <= s.Radius*s.Radius
}

// IntersectsSphere reports whether the sphere overlaps the other sphere,
// touching surfaces included.
func (s Sphere) IntersectsSphere(other Sphere) bool {
	rr := s.Radius + other.Radius
	return s.Center.DistanceSq(other.Center) <= rr*rr
}

// IntersectsAABB reports whether the sphere overlaps the axis-aligned box bx,
// by comparing the box's closest point to the sphere's radius.
func (s Sphere) IntersectsAABB(bx AABB) bool {
	closest := bx.ClosestPoint(s.Center)
	return closest.DistanceSq(s.Center) <= s.Radius*s.Radius
}

// IntersectRay intersects a ray with the sphere. On a hit it returns the entry
// parameter tmin and exit parameter tmax along the ray (both may coincide for a
// tangent ray) and true. If the origin is inside the sphere tmin is negative.
// It returns false when the ray misses the sphere or only meets it behind the
// origin. The ray direction is assumed to be unit length.
func (s Sphere) IntersectRay(r Ray) (tmin, tmax float64, ok bool) {
	oc := r.Origin.Sub(s.Center)
	a := r.Dir.Dot(r.Dir)
	b := 2 * oc.Dot(r.Dir)
	c := oc.Dot(oc) - s.Radius*s.Radius
	disc := b*b - 4*a*c
	if disc < 0 {
		return 0, 0, false
	}
	sq := math.Sqrt(disc)
	t0 := (-b - sq) / (2 * a)
	t1 := (-b + sq) / (2 * a)
	if t1 < 0 {
		return 0, 0, false
	}
	return t0, t1, true
}
