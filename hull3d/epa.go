package hull3d

import (
	"errors"
	"math"
)

// Penetration describes the result of an [EPAPenetration] query for two
// overlapping convex shapes.
type Penetration struct {
	// Depth is the minimum translation distance to separate the shapes.
	Depth float64
	// Normal is the unit contact normal (the direction to push shape A out of
	// shape B).
	Normal Vec3
	// PointA and PointB are witness points on the boundary of A and B.
	PointA, PointB Vec3
}

// epaFace is a triangular face of the expanding polytope over the Minkowski
// difference, with a cached outward plane through the origin.
type epaFace struct {
	v      [3]int
	normal Vec3
	dist   float64
}

// EPAPenetration computes the penetration depth and contact normal of two
// overlapping convex shapes with the Expanding Polytope Algorithm. It first runs
// GJK to obtain a tetrahedron enclosing the origin of the Minkowski difference,
// then repeatedly expands the closest face toward the boundary. It returns an
// error if the shapes do not actually overlap.
func EPAPenetration(a, b ConvexShape) (Penetration, error) {
	ok, simplex := gjkSimplex(a, b)
	if !ok || len(simplex) < 4 {
		return Penetration{}, errors.New("hull3d: shapes do not overlap")
	}

	verts := append([]mvert(nil), simplex[:4]...)
	// Ensure positive orientation of the seed tetrahedron.
	if Orient3D(verts[0].p, verts[1].p, verts[2].p, verts[3].p) > 0 {
		verts[2], verts[3] = verts[3], verts[2]
	}
	faces := []epaFace{
		{v: [3]int{0, 1, 2}},
		{v: [3]int{0, 2, 3}},
		{v: [3]int{0, 3, 1}},
		{v: [3]int{1, 3, 2}},
	}
	for i := range faces {
		computeFace(&faces[i], verts)
	}

	const maxIter = 128
	const tol = 1e-9
	for iter := 0; iter < maxIter; iter++ {
		// Closest face to the origin.
		ci := 0
		for i := 1; i < len(faces); i++ {
			if faces[i].dist < faces[ci].dist {
				ci = i
			}
		}
		fc := faces[ci]
		w := mdSupport(a, b, fc.normal)
		d := w.p.Dot(fc.normal)
		if d-fc.dist < tol {
			return finishEPA(fc, verts, d), nil
		}

		// Remove faces visible from w and stitch the horizon to the new vertex.
		newIdx := len(verts)
		verts = append(verts, w)
		type edge struct{ a, b int }
		horizon := map[edge]int{}
		var kept []epaFace
		for _, f := range faces {
			if f.normal.Dot(w.p.Sub(verts[f.v[0]].p)) > 0 {
				es := [3]edge{{f.v[0], f.v[1]}, {f.v[1], f.v[2]}, {f.v[2], f.v[0]}}
				for _, e := range es {
					horizon[e]++
				}
				continue
			}
			kept = append(kept, f)
		}
		faces = kept
		for e, cnt := range horizon {
			if cnt != 1 {
				continue // interior edge shared by two removed faces
			}
			if _, twin := horizon[edge{e.b, e.a}]; twin {
				continue
			}
			nf := epaFace{v: [3]int{e.a, e.b, newIdx}}
			computeFace(&nf, verts)
			faces = append(faces, nf)
		}
		if len(faces) == 0 {
			return Penetration{}, errors.New("hull3d: EPA degenerated")
		}
	}
	// Return best-so-far on non-convergence.
	ci := 0
	for i := 1; i < len(faces); i++ {
		if faces[i].dist < faces[ci].dist {
			ci = i
		}
	}
	fc := faces[ci]
	return finishEPA(fc, verts, fc.dist), nil
}

func computeFace(f *epaFace, verts []mvert) {
	a := verts[f.v[0]].p
	b := verts[f.v[1]].p
	c := verts[f.v[2]].p
	n := b.Sub(a).Cross(c.Sub(a))
	if l := n.Length(); l > 1e-300 {
		n = n.Scale(1 / l)
	}
	// Orient the normal to point away from the origin.
	if n.Dot(a) < 0 {
		n = n.Neg()
		f.v[1], f.v[2] = f.v[2], f.v[1]
	}
	f.normal = n
	f.dist = n.Dot(a)
	if f.dist < 0 {
		f.dist = 0
	}
}

func finishEPA(fc epaFace, verts []mvert, depth float64) Penetration {
	a := verts[fc.v[0]]
	b := verts[fc.v[1]]
	c := verts[fc.v[2]]
	// Barycentric coordinates of the projection of the origin onto the face.
	u, v, w := baryOnTriangle(fc.normal.Scale(fc.dist), a.p, b.p, c.p)
	pa := a.sa.Scale(u).Add(b.sa.Scale(v)).Add(c.sa.Scale(w))
	pb := a.sb.Scale(u).Add(b.sb.Scale(v)).Add(c.sb.Scale(w))
	return Penetration{
		Depth:  math.Max(depth, fc.dist),
		Normal: fc.normal,
		PointA: pa,
		PointB: pb,
	}
}

// baryOnTriangle returns the barycentric coordinates of point p with respect to
// triangle (a,b,c). p is assumed to lie in (or near) the triangle's plane.
func baryOnTriangle(p, a, b, c Vec3) (u, v, w float64) {
	v0 := b.Sub(a)
	v1 := c.Sub(a)
	v2 := p.Sub(a)
	d00 := v0.Dot(v0)
	d01 := v0.Dot(v1)
	d11 := v1.Dot(v1)
	d20 := v2.Dot(v0)
	d21 := v2.Dot(v1)
	den := d00*d11 - d01*d01
	if math.Abs(den) < 1e-300 {
		return 1, 0, 0
	}
	v = (d11*d20 - d01*d21) / den
	w = (d00*d21 - d01*d20) / den
	u = 1 - v - w
	return u, v, w
}
