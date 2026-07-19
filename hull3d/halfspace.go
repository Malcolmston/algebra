package hull3d

import (
	"errors"
	"math"
)

// HalfSpace is a closed half-space, the feasible side of a [Plane]. The
// convention here is that the feasible region is {x : Plane.Eval(x) <= 0}, i.e.
// Normal·x <= Offset; the outward normal points away from the feasible side.
type HalfSpace struct {
	Plane Plane
}

// NewHalfSpace returns the half-space {x : normal·x <= offset}.
func NewHalfSpace(normal Vec3, offset float64) HalfSpace {
	return HalfSpace{Plane{normal, offset}}
}

// HalfSpaceFromPointNormal returns the half-space with the given boundary point
// and outward normal: {x : normal·(x-p) <= 0}.
func HalfSpaceFromPointNormal(p, normal Vec3) HalfSpace {
	return HalfSpace{PlaneFromPointNormal(p, normal)}
}

// Contains reports whether q lies in the (closed, eps-relaxed) half-space.
func (h HalfSpace) Contains(q Vec3, eps float64) bool {
	return h.Plane.SignedDistance(q) <= eps
}

// HalfSpaces returns one outward half-space per face of the polytope, so that
// the polytope equals the intersection of the returned half-spaces (for a convex
// outward-oriented polytope).
func (p *Polytope) HalfSpaces() []HalfSpace {
	var out []HalfSpace
	for _, f := range p.Faces {
		if pl, err := p.FacePlane(f); err == nil {
			out = append(out, HalfSpace{pl})
		}
	}
	return out
}

// HalfSpaceIntersection computes the bounded convex polytope that is the
// intersection of the given half-spaces, using the polar-dual construction:
// with the supplied strictly interior point moved to the origin, each half-space
// Normal·y <= b (b>0) maps to the dual point Normal/b, the intersection is the
// polar dual of the convex hull of those dual points, and each facet of that
// hull yields a vertex of the intersection.
//
// The interior point must satisfy every constraint strictly. The routine returns
// an error if the interior point is infeasible or if the intersection is
// unbounded (the dual hull does not contain the origin in its interior).
func HalfSpaceIntersection(halfspaces []HalfSpace, interior Vec3) (*Polytope, error) {
	if len(halfspaces) < 4 {
		return nil, errors.New("hull3d: need at least four half-spaces for a bounded body")
	}
	dual := make([]Vec3, 0, len(halfspaces))
	for _, hs := range halfspaces {
		b := hs.Plane.Offset - hs.Plane.Normal.Dot(interior)
		if b <= 1e-12 {
			return nil, errors.New("hull3d: interior point is not strictly feasible")
		}
		dual = append(dual, hs.Plane.Normal.Scale(1/b))
	}

	q, err := ConvexHull(dual)
	if err != nil {
		return nil, errors.New("hull3d: half-space intersection is unbounded or degenerate")
	}

	var verts []Vec3
	for _, f := range q.Faces {
		pl, err := q.FacePlane(f)
		if err != nil {
			continue
		}
		// Face plane is Normal·z = Offset (outward). Origin must be interior,
		// so Offset > 0; the corresponding primal vertex is Normal/Offset.
		if pl.Offset <= 1e-12 {
			return nil, errors.New("hull3d: half-space intersection is unbounded")
		}
		verts = append(verts, pl.Normal.Scale(1/pl.Offset).Add(interior))
	}

	verts = DedupPoints(verts, 1e-9)
	return ConvexHull(verts)
}

// IntersectHalfSpacesWithBox intersects the given half-spaces together with the
// six faces of a large axis-aligned bounding box centred at interior with the
// given half-extent, guaranteeing a bounded result. It is a convenience wrapper
// that clips an otherwise unbounded region to a finite box.
func IntersectHalfSpacesWithBox(halfspaces []HalfSpace, interior Vec3, extent float64) (*Polytope, error) {
	if extent <= 0 {
		return nil, errors.New("hull3d: box extent must be positive")
	}
	all := append([]HalfSpace(nil), halfspaces...)
	dirs := []Vec3{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}
	for _, d := range dirs {
		all = append(all, HalfSpaceFromPointNormal(interior.Add(d.Scale(extent)), d))
	}
	return HalfSpaceIntersection(all, interior)
}

// FeasiblePoint returns a strictly interior point of the intersection of the
// half-spaces, if one is found by the Chebyshev-centre heuristic of averaging a
// small random-free search over vertex-like directions. It returns an error if
// no strictly feasible point is found among the tried candidates. The candidate
// guess is a good seed for [HalfSpaceIntersection].
func FeasiblePoint(halfspaces []HalfSpace) (Vec3, error) {
	if len(halfspaces) == 0 {
		return Vec3{}, errors.New("hull3d: no half-spaces")
	}
	// Try the average of boundary base points as an initial guess, then a small
	// deterministic grid of offsets.
	var base Vec3
	m := 0
	for _, hs := range halfspaces {
		if pt, err := hs.Plane.PointOnPlane(); err == nil {
			base = base.Add(pt)
			m++
		}
	}
	if m > 0 {
		base = base.Scale(1 / float64(m))
	}
	feasible := func(q Vec3) bool {
		for _, hs := range halfspaces {
			if hs.Plane.SignedDistance(q) > -1e-9 {
				return false
			}
		}
		return true
	}
	if feasible(base) {
		return base, nil
	}
	// Deterministic search over a lattice around the base.
	for _, s := range []float64{0.1, 0.5, 1, 2, 5, 10} {
		for dx := -1; dx <= 1; dx++ {
			for dy := -1; dy <= 1; dy++ {
				for dz := -1; dz <= 1; dz++ {
					q := base.Add(Vec3{float64(dx), float64(dy), float64(dz)}.Scale(s))
					if feasible(q) {
						return q, nil
					}
				}
			}
		}
	}
	return Vec3{}, errors.New("hull3d: no strictly feasible point found")
}

// ChebyshevCenterEstimate returns an estimate of the deepest interior point and
// its clearance for the intersection of the half-spaces, by hill-climbing the
// minimum signed clearance from [FeasiblePoint]. The normals need not be unit
// length; distances use the true Euclidean clearance. It returns an error if no
// feasible seed exists.
func ChebyshevCenterEstimate(halfspaces []HalfSpace) (center Vec3, radius float64, err error) {
	seed, err := FeasiblePoint(halfspaces)
	if err != nil {
		return Vec3{}, 0, err
	}
	clearance := func(q Vec3) float64 {
		r := math.Inf(1)
		for _, hs := range halfspaces {
			if d := -hs.Plane.SignedDistance(q); d < r {
				r = d
			}
		}
		return r
	}
	center = seed
	best := clearance(center)
	step := 1.0
	dirs := []Vec3{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}
	for i := 0; i < 200 && step > 1e-9; i++ {
		improved := false
		for _, d := range dirs {
			cand := center.Add(d.Scale(step))
			if c := clearance(cand); c > best {
				best, center, improved = c, cand, true
			}
		}
		if !improved {
			step *= 0.5
		}
	}
	return center, best, nil
}
