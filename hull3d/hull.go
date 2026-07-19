package hull3d

import (
	"errors"
	"math"
	"math/rand"
	"sort"
)

// hullFace is the internal mutable face used during incremental construction.
type hullFace struct {
	v      [3]int
	plane  Plane
	alive  bool
	points []int // indices of not-yet-processed points on the outer side
}

// ConvexHull computes the convex hull of the given points and returns it as a
// closed, outward-oriented triangulated [Polytope]. Duplicate and interior
// points are discarded; the returned polytope's vertices are a subset of the
// input (re-indexed). It returns an error if fewer than four points are supplied
// or the points are coplanar (do not span three dimensions).
//
// The algorithm is the incremental / QuickHull method: an initial tetrahedron of
// extreme points is built, then remaining points are added one at a time, each
// time removing the faces it can see and stitching new faces across the horizon.
func ConvexHull(points []Vec3) (*Polytope, error) {
	return convexHull(points, 0)
}

// QuickHull is an alias for [ConvexHull]; the incremental construction used here
// selects, for each face, the farthest outside point to add next, which is the
// distinguishing step of the QuickHull algorithm.
func QuickHull(points []Vec3) (*Polytope, error) {
	return convexHull(points, 0)
}

// ConvexHullSeeded behaves like [ConvexHull] but shuffles the input point order
// using the supplied seed before construction. Shuffling can improve numerical
// behaviour on adversarial inputs; the resulting hull is identical as a set.
func ConvexHullSeeded(points []Vec3, seed int64) (*Polytope, error) {
	pts := append([]Vec3(nil), points...)
	r := rand.New(rand.NewSource(seed))
	r.Shuffle(len(pts), func(i, j int) { pts[i], pts[j] = pts[j], pts[i] })
	return convexHull(pts, 0)
}

func convexHull(points []Vec3, epsScale float64) (*Polytope, error) {
	pts := DedupPoints(points, 1e-12)
	if len(pts) < 4 {
		return nil, errors.New("hull3d: need at least four distinct points")
	}

	// Scale-aware epsilon based on the bounding box extent.
	min, max, _ := BoundingBox(pts)
	extent := max.Sub(min).LInf()
	if extent < 1e-300 {
		return nil, errNoInterior
	}
	eps := extent * 1e-9
	if epsScale > 0 {
		eps = extent * epsScale
	}

	i0, i1, i2, i3, err := initialSimplex(pts, eps)
	if err != nil {
		return nil, err
	}

	interior := pts[i0].Add(pts[i1]).Add(pts[i2]).Add(pts[i3]).Scale(0.25)

	var faces []*hullFace
	newFace := func(a, b, c int) *hullFace {
		pl, _ := SupportingPlane(pts[a], pts[b], pts[c], interior)
		// Order the triple so (B-A)x(C-A) matches the outward plane normal.
		va, vb, vc := pts[a], pts[b], pts[c]
		if vb.Sub(va).Cross(vc.Sub(va)).Dot(pl.Normal) < 0 {
			b, c = c, b
		}
		return &hullFace{v: [3]int{a, b, c}, plane: pl, alive: true}
	}
	faces = append(faces,
		newFace(i0, i1, i2),
		newFace(i0, i1, i3),
		newFace(i0, i2, i3),
		newFace(i1, i2, i3),
	)

	// Assign each remaining point to the outside set of the first face it sees.
	assigned := map[int]bool{i0: true, i1: true, i2: true, i3: true}
	for idx := range pts {
		if assigned[idx] {
			continue
		}
		for _, f := range faces {
			if f.plane.SignedDistance(pts[idx]) > eps {
				f.points = append(f.points, idx)
				break
			}
		}
	}

	// Process faces with outside points until none remain.
	for {
		var work *hullFace
		for _, f := range faces {
			if f.alive && len(f.points) > 0 {
				work = f
				break
			}
		}
		if work == nil {
			break
		}

		// Farthest outside point of this face becomes the new apex.
		apex := work.points[0]
		best := work.plane.SignedDistance(pts[apex])
		for _, pi := range work.points[1:] {
			if d := work.plane.SignedDistance(pts[pi]); d > best {
				best, apex = d, pi
			}
		}

		// Find all faces visible from the apex.
		var visible []*hullFace
		for _, f := range faces {
			if f.alive && f.plane.SignedDistance(pts[apex]) > eps {
				visible = append(visible, f)
			}
		}
		if len(visible) == 0 {
			// Numerical: point already inside; drop it.
			work.points = removeInt(work.points, apex)
			continue
		}

		visibleSet := make(map[*hullFace]bool, len(visible))
		for _, f := range visible {
			visibleSet[f] = true
		}

		// Collect horizon: directed edges of visible faces whose twin lies in a
		// non-visible face.
		type dedge struct{ a, b int }
		var horizon []dedge
		for _, f := range visible {
			edges := [3]dedge{{f.v[0], f.v[1]}, {f.v[1], f.v[2]}, {f.v[2], f.v[0]}}
			for _, e := range edges {
				if !edgeInVisible(faces, visibleSet, e.b, e.a) {
					horizon = append(horizon, e)
				}
			}
		}

		// Gather orphaned outside points from the visible faces.
		var orphans []int
		for _, f := range visible {
			orphans = append(orphans, f.points...)
			f.alive = false
			f.points = nil
		}

		// Build new faces from horizon edges to the apex.
		var created []*hullFace
		for _, e := range horizon {
			nf := newFace(e.a, e.b, apex)
			faces = append(faces, nf)
			created = append(created, nf)
		}

		// Redistribute orphaned points to the new faces.
		seen := make(map[int]bool)
		for _, pi := range orphans {
			if pi == apex || seen[pi] {
				continue
			}
			seen[pi] = true
			for _, nf := range created {
				if nf.plane.SignedDistance(pts[pi]) > eps {
					nf.points = append(nf.points, pi)
					break
				}
			}
		}
	}

	return assembleHull(pts, faces)
}

// edgeInVisible reports whether the directed edge (a,b) belongs to some alive
// face that is in the visible set.
func edgeInVisible(faces []*hullFace, visible map[*hullFace]bool, a, b int) bool {
	for _, f := range faces {
		if !f.alive || !visible[f] {
			continue
		}
		if (f.v[0] == a && f.v[1] == b) ||
			(f.v[1] == a && f.v[2] == b) ||
			(f.v[2] == a && f.v[0] == b) {
			return true
		}
	}
	return false
}

func removeInt(s []int, x int) []int {
	out := s[:0]
	for _, v := range s {
		if v != x {
			out = append(out, v)
		}
	}
	return out
}

// assembleHull compacts alive faces into a Polytope with re-indexed vertices.
func assembleHull(pts []Vec3, faces []*hullFace) (*Polytope, error) {
	remap := make(map[int]int)
	var verts []Vec3
	idxOf := func(i int) int {
		if j, ok := remap[i]; ok {
			return j
		}
		j := len(verts)
		remap[i] = j
		verts = append(verts, pts[i])
		return j
	}
	var outFaces []Face
	for _, f := range faces {
		if !f.alive {
			continue
		}
		outFaces = append(outFaces, Face{idxOf(f.v[0]), idxOf(f.v[1]), idxOf(f.v[2])})
	}
	if len(outFaces) < 4 || len(verts) < 4 {
		return nil, errNoInterior
	}
	return &Polytope{Vertices: verts, Faces: outFaces}, nil
}

// initialSimplex chooses four affinely independent points forming a
// non-degenerate tetrahedron and returns their indices.
func initialSimplex(pts []Vec3, eps float64) (i0, i1, i2, i3 int, err error) {
	// Extreme points along each axis.
	var minI, maxI [3]int
	for a := 0; a < 3; a++ {
		for i := 1; i < len(pts); i++ {
			if pts[i].Get(a) < pts[minI[a]].Get(a) {
				minI[a] = i
			}
			if pts[i].Get(a) > pts[maxI[a]].Get(a) {
				maxI[a] = i
			}
		}
	}
	// Pick the two extreme points farthest apart.
	i0, i1 = 0, 1
	best := -1.0
	cand := []int{minI[0], maxI[0], minI[1], maxI[1], minI[2], maxI[2]}
	for a := 0; a < len(cand); a++ {
		for b := a + 1; b < len(cand); b++ {
			if d := pts[cand[a]].DistanceSq(pts[cand[b]]); d > best {
				best, i0, i1 = d, cand[a], cand[b]
			}
		}
	}
	if best <= eps*eps {
		return 0, 0, 0, 0, errNoInterior
	}

	// Third point maximises distance to the line i0-i1.
	dir := pts[i1].Sub(pts[i0])
	best = -1.0
	i2 = -1
	for i := range pts {
		d := dir.Cross(pts[i].Sub(pts[i0])).LengthSq()
		if d > best {
			best, i2 = d, i
		}
	}
	if i2 < 0 || best <= eps*eps {
		return 0, 0, 0, 0, errNoInterior
	}

	// Fourth point maximises absolute distance to the plane i0-i1-i2.
	pl, err2 := PlaneFromPoints(pts[i0], pts[i1], pts[i2])
	if err2 != nil {
		return 0, 0, 0, 0, errNoInterior
	}
	best = -1.0
	i3 = -1
	for i := range pts {
		d := math.Abs(pl.SignedDistance(pts[i]))
		if d > best {
			best, i3 = d, i
		}
	}
	if i3 < 0 || best <= eps {
		return 0, 0, 0, 0, errNoInterior
	}
	return i0, i1, i2, i3, nil
}

// GiftWrapHull computes the convex hull by enumerating supporting facets and
// returns it as a closed, outward-oriented triangulated [Polytope]. It is an
// O(n^4) brute-force method whose logic is entirely independent of the
// incremental [ConvexHull], and is intended as a cross-check on small inputs.
// Coplanar boundary points are merged into a single planar facet and then
// fan-triangulated, so the surface is manifold even for boxes and prisms. It
// returns an error under the same degenerate conditions as [ConvexHull].
func GiftWrapHull(points []Vec3) (*Polytope, error) {
	pts := DedupPoints(points, 1e-12)
	if len(pts) < 4 {
		return nil, errors.New("hull3d: need at least four distinct points")
	}
	min, max, _ := BoundingBox(pts)
	extent := max.Sub(min).LInf()
	if extent < 1e-300 {
		return nil, errNoInterior
	}
	eps := extent * 1e-9
	i0, i1, i2, i3, err := initialSimplex(pts, eps)
	if err != nil {
		return nil, err
	}
	interior := pts[i0].Add(pts[i1]).Add(pts[i2]).Add(pts[i3]).Scale(0.25)

	// Enumerate distinct supporting planes: for each non-collinear triple whose
	// outward plane keeps every point on the closed inner side, record the plane
	// keyed by a rounded (unit-normal, offset) signature so coplanar triples map
	// together.
	type facet struct {
		plane Plane
		verts map[int]bool
	}
	facets := make(map[[4]int64]*facet)
	round := func(x float64) int64 { return int64(math.Round(x / (eps * 10))) }
	n := len(pts)
	for a := 0; a < n; a++ {
		for b := a + 1; b < n; b++ {
			for c := b + 1; c < n; c++ {
				pl, e := PlaneFromPoints(pts[a], pts[b], pts[c])
				if e != nil {
					continue
				}
				npl, e2 := pl.Normalized()
				if e2 != nil {
					continue
				}
				if npl.Eval(interior) > 0 {
					npl = npl.Flip()
				}
				allInside := true
				for k := 0; k < n; k++ {
					if npl.SignedDistance(pts[k]) > eps {
						allInside = false
						break
					}
				}
				if !allInside {
					continue
				}
				key := [4]int64{round(npl.Normal.X), round(npl.Normal.Y), round(npl.Normal.Z), round(npl.Offset)}
				f := facets[key]
				if f == nil {
					f = &facet{plane: npl, verts: make(map[int]bool)}
					facets[key] = f
				}
				for k := 0; k < n; k++ {
					if math.Abs(npl.SignedDistance(pts[k])) <= eps {
						f.verts[k] = true
					}
				}
			}
		}
	}
	if len(facets) < 4 {
		return nil, errNoInterior
	}

	remap := make(map[int]int)
	var verts []Vec3
	idxOf := func(i int) int {
		if j, ok := remap[i]; ok {
			return j
		}
		j := len(verts)
		remap[i] = j
		verts = append(verts, pts[i])
		return j
	}

	var faces []Face
	for _, f := range facets {
		ring := orderCoplanar(pts, f.verts, f.plane)
		if len(ring) < 3 {
			continue
		}
		// Fan-triangulate the convex facet polygon.
		a := idxOf(ring[0])
		for i := 1; i+1 < len(ring); i++ {
			faces = append(faces, Face{a, idxOf(ring[i]), idxOf(ring[i+1])})
		}
	}
	if len(faces) < 4 {
		return nil, errNoInterior
	}
	return &Polytope{Vertices: verts, Faces: faces}, nil
}

// orderCoplanar returns the vertex indices in set, ordered counter-clockwise
// (as seen from the outward side of pl) around their centroid within the plane.
func orderCoplanar(pts []Vec3, set map[int]bool, pl Plane) []int {
	idx := make([]int, 0, len(set))
	for i := range set {
		idx = append(idx, i)
	}
	if len(idx) < 3 {
		return idx
	}
	var c Vec3
	for _, i := range idx {
		c = c.Add(pts[i])
	}
	c = c.Scale(1 / float64(len(idx)))
	u, v, err := pl.Normal.OrthonormalBasis()
	if err != nil {
		return idx
	}
	sort.Slice(idx, func(i, j int) bool {
		pi := pts[idx[i]].Sub(c)
		pj := pts[idx[j]].Sub(c)
		ai := math.Atan2(pi.Dot(v), pi.Dot(u))
		aj := math.Atan2(pj.Dot(v), pj.Dot(u))
		return ai < aj
	})
	return idx
}

// HullVolume computes the volume of the convex hull of the given points.
func HullVolume(points []Vec3) (float64, error) {
	h, err := ConvexHull(points)
	if err != nil {
		return 0, err
	}
	return h.Volume(), nil
}

// HullSurfaceArea computes the surface area of the convex hull of the given
// points.
func HullSurfaceArea(points []Vec3) (float64, error) {
	h, err := ConvexHull(points)
	if err != nil {
		return 0, err
	}
	return h.SurfaceArea(), nil
}

// ExtremePoints returns the subset of input points that are vertices of the
// convex hull.
func ExtremePoints(points []Vec3) ([]Vec3, error) {
	h, err := ConvexHull(points)
	if err != nil {
		return nil, err
	}
	return h.Vertices, nil
}

// IsConvexPosition reports whether p is a vertex (extreme point) of the convex
// hull of pts, i.e. it is not contained in the hull of the remaining points.
func IsConvexPosition(p Vec3, pts []Vec3, eps float64) bool {
	var rest []Vec3
	for _, q := range pts {
		if !q.ApproxEqual(p, eps) {
			rest = append(rest, q)
		}
	}
	h, err := ConvexHull(rest)
	if err != nil {
		return true
	}
	return !h.Contains(p, eps)
}
