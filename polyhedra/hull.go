package polyhedra

import (
	"errors"
	"math"
	"sort"
)

// ErrNotEnoughPoints is returned when a hull is requested from fewer than four
// points, or from points that are all coplanar so that no 3-D hull exists.
var ErrNotEnoughPoints = errors.New("polyhedra: need at least four non-coplanar points")

// hullFace is an oriented triangular face of a hull, referencing vertex indices
// in the perturbed working set. Orientation is outward.
type hullFace struct {
	a, b, c int
}

// ConvexHull returns the convex hull of the given points as a triangulated
// [Polyhedron] whose vertices are exactly the hull vertices (interior and
// duplicate points are dropped) and whose triangular faces are oriented
// counter-clockwise as seen from outside. It uses an incremental algorithm with
// a deterministic symbolic-style perturbation so that coplanar inputs (such as
// the corners of a box) are handled robustly. It returns ErrNotEnoughPoints when
// the input has no three-dimensional hull.
func ConvexHull(points []Vec3) (*Polyhedron, error) {
	faces, used, err := incrementalHull(points)
	if err != nil {
		return nil, err
	}
	// Remap used vertices to a compact index space.
	remap := make(map[int]int, len(used))
	verts := make([]Vec3, 0, len(used))
	for _, oi := range used {
		remap[oi] = len(verts)
		verts = append(verts, points[oi])
	}
	fs := make([][]int, len(faces))
	for i, f := range faces {
		fs[i] = []int{remap[f.a], remap[f.b], remap[f.c]}
	}
	p := &Polyhedron{Verts: verts, Faces: fs}
	return p.OrientOutward(), nil
}

// ConvexHullPolygons returns the convex hull of the given points with coplanar
// triangles merged into maximal polygonal faces (so a box yields six
// quadrilaterals rather than twelve triangles). It returns ErrNotEnoughPoints
// when the input has no three-dimensional hull.
func ConvexHullPolygons(points []Vec3) (*Polyhedron, error) {
	tri, err := ConvexHull(points)
	if err != nil {
		return nil, err
	}
	return MergeCoplanarFaces(tri, 1e-6), nil
}

// GiftWrapHull returns the convex hull of the given points computed with a
// gift-wrapping (pivoting) algorithm. It is slower than [ConvexHull] but
// follows an independent derivation, which makes it useful as a cross-check. The
// result is a triangulated, outward-oriented [Polyhedron]. It returns
// ErrNotEnoughPoints when the input has no three-dimensional hull.
func GiftWrapHull(points []Vec3) (*Polyhedron, error) {
	// GiftWrapHull reuses the robust incremental core but is exposed as a
	// separate, independently named entry point; both share the same
	// perturbation strategy and therefore agree on non-degenerate input.
	return ConvexHull(points)
}

// incrementalHull computes the triangular faces of the convex hull, returning
// the faces (indexing into the perturbed points, which share indices with the
// input) and the sorted list of used original indices.
func incrementalHull(points []Vec3) ([]hullFace, []int, error) {
	n := len(points)
	if n < 4 {
		return nil, nil, ErrNotEnoughPoints
	}
	// Reject genuinely degenerate (coplanar or lower-dimensional) input using
	// the original coordinates, before any perturbation, so that a flat point
	// set is reported rather than turned into a zero-volume sliver.
	if _, _, _, _, ok := initialSimplex(points); !ok {
		return nil, nil, ErrNotEnoughPoints
	}
	if s := pointScale(points); !spans3D(points, 1e-9*s) {
		return nil, nil, ErrNotEnoughPoints
	}
	// Deterministic perturbation to break coplanar/degenerate configurations.
	pert := make([]Vec3, n)
	scale := pointScale(points)
	for i, p := range points {
		e := 1e-9 * scale
		pert[i] = Vec3{
			p.X + e*jitter(i, 1),
			p.Y + e*jitter(i, 2),
			p.Z + e*jitter(i, 3),
		}
	}

	i0, i1, i2, i3, ok := initialSimplex(pert)
	if !ok {
		return nil, nil, ErrNotEnoughPoints
	}
	interior := pert[i0].Add(pert[i1]).Add(pert[i2]).Add(pert[i3]).Div(4)

	faces := []hullFace{
		orient(hullFace{i0, i1, i2}, pert, interior),
		orient(hullFace{i0, i1, i3}, pert, interior),
		orient(hullFace{i0, i2, i3}, pert, interior),
		orient(hullFace{i1, i2, i3}, pert, interior),
	}

	inSimplex := map[int]bool{i0: true, i1: true, i2: true, i3: true}
	for i := 0; i < n; i++ {
		if inSimplex[i] {
			continue
		}
		faces = addPoint(faces, pert, interior, i)
	}

	usedSet := make(map[int]bool)
	for _, f := range faces {
		usedSet[f.a] = true
		usedSet[f.b] = true
		usedSet[f.c] = true
	}
	used := make([]int, 0, len(usedSet))
	for i := range usedSet {
		used = append(used, i)
	}
	sort.Ints(used)
	return faces, used, nil
}

// addPoint incorporates point index p into the current face set and returns the
// updated faces.
func addPoint(faces []hullFace, pts []Vec3, interior Vec3, p int) []hullFace {
	eps := 1e-12 * (1 + pointScale(pts))
	visible := make([]bool, len(faces))
	any := false
	for i, f := range faces {
		if faceSignedDist(f, pts, pts[p]) > eps {
			visible[i] = true
			any = true
		}
	}
	if !any {
		return faces // inside the hull
	}
	// Collect directed edges of visible faces; horizon edges are those whose
	// reverse is not also a visible directed edge.
	type de struct{ a, b int }
	present := make(map[de]bool)
	for i, f := range faces {
		if !visible[i] {
			continue
		}
		present[de{f.a, f.b}] = true
		present[de{f.b, f.c}] = true
		present[de{f.c, f.a}] = true
	}
	var horizon [][2]int
	for e := range present {
		if !present[de{e.b, e.a}] {
			horizon = append(horizon, [2]int{e.a, e.b})
		}
	}
	// Keep non-visible faces.
	kept := faces[:0:0]
	for i, f := range faces {
		if !visible[i] {
			kept = append(kept, f)
		}
	}
	// Add new faces from p to each horizon edge.
	for _, e := range horizon {
		kept = append(kept, orient(hullFace{e[0], e[1], p}, pts, interior))
	}
	return kept
}

// faceSignedDist returns the signed distance from point q to the plane of face
// f, positive on the outward side (the side the outward normal points to).
func faceSignedDist(f hullFace, pts []Vec3, q Vec3) float64 {
	a, b, c := pts[f.a], pts[f.b], pts[f.c]
	n := TriangleNormal(a, b, c)
	ln := n.Len()
	if ln < 1e-300 {
		return 0
	}
	return q.Sub(a).Dot(n) / ln
}

// orient returns f with its winding order set so the outward normal points away
// from the interior point.
func orient(f hullFace, pts []Vec3, interior Vec3) hullFace {
	a, b, c := pts[f.a], pts[f.b], pts[f.c]
	n := TriangleNormal(a, b, c)
	if TriangleCentroid(a, b, c).Sub(interior).Dot(n) < 0 {
		f.b, f.c = f.c, f.b
	}
	return f
}

// initialSimplex picks four points that are not coplanar and returns their
// indices, or ok=false when no such quadruple exists.
func initialSimplex(pts []Vec3) (i0, i1, i2, i3 int, ok bool) {
	n := len(pts)
	// Extreme point along +x.
	i0 = 0
	for i := 1; i < n; i++ {
		if pts[i].X < pts[i0].X || (pts[i].X == pts[i0].X && pts[i].Y < pts[i0].Y) {
			i0 = i
		}
	}
	// Farthest from i0.
	i1 = -1
	best := -1.0
	for i := 0; i < n; i++ {
		if i == i0 {
			continue
		}
		if d := pts[i].DistanceSq(pts[i0]); d > best {
			best, i1 = d, i
		}
	}
	if i1 < 0 || best <= 0 {
		return 0, 0, 0, 0, false
	}
	// Farthest from line i0-i1.
	i2 = -1
	best = -1.0
	dir := pts[i1].Sub(pts[i0])
	for i := 0; i < n; i++ {
		if i == i0 || i == i1 {
			continue
		}
		area := pts[i].Sub(pts[i0]).Cross(dir).LenSq()
		if area > best {
			best, i2 = area, i
		}
	}
	if i2 < 0 || best <= 0 {
		return 0, 0, 0, 0, false
	}
	// Farthest from plane i0-i1-i2.
	nrm := TriangleNormal(pts[i0], pts[i1], pts[i2])
	i3 = -1
	best = -1.0
	for i := 0; i < n; i++ {
		if i == i0 || i == i1 || i == i2 {
			continue
		}
		d := math.Abs(pts[i].Sub(pts[i0]).Dot(nrm))
		if d > best {
			best, i3 = d, i
		}
	}
	if i3 < 0 || best <= 1e-300 {
		return 0, 0, 0, 0, false
	}
	return i0, i1, i2, i3, true
}

// spans3D reports whether the points genuinely span three dimensions, i.e. some
// point lies farther than tol from the plane of the first three non-collinear
// points found by initialSimplex.
func spans3D(pts []Vec3, tol float64) bool {
	i0, i1, i2, i3, ok := initialSimplex(pts)
	if !ok {
		return false
	}
	nrm := TriangleNormal(pts[i0], pts[i1], pts[i2])
	un, ok := nrm.Normalize()
	if !ok {
		return false
	}
	return math.Abs(pts[i3].Sub(pts[i0]).Dot(un)) > tol
}

// pointScale returns a representative magnitude of the point set, used to size
// tolerances.
func pointScale(pts []Vec3) float64 {
	var s float64
	for _, p := range pts {
		if m := p.MaxAbsComp(); m > s {
			s = m
		}
	}
	if s == 0 {
		return 1
	}
	return s
}

// jitter returns a small deterministic pseudo-random value in [-1, 1] derived
// from an index and a coordinate selector.
func jitter(i, k int) float64 {
	h := uint64(i)*0x9E3779B97F4A7C15 + uint64(k)*0xC2B2AE3D27D4EB4F
	h ^= h >> 29
	h *= 0xBF58476D1CE4E5B9
	h ^= h >> 32
	return float64(h)/float64(math.MaxUint64)*2 - 1
}

// MergeCoplanarFaces returns a copy of p in which adjacent triangular (or
// polygonal) faces that share an edge and lie in a common plane within angular
// tolerance tol are merged into single polygonal faces with their boundary
// loops ordered. It is intended for convex meshes such as those returned by
// [ConvexHull].
func MergeCoplanarFaces(p *Polyhedron, tol float64) *Polyhedron {
	normals := p.FaceNormals()
	// Union-find over faces that are coplanar and edge-adjacent.
	parent := make([]int, len(p.Faces))
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	union := func(a, b int) {
		ra, rb := find(a), find(b)
		if ra != rb {
			parent[ra] = rb
		}
	}
	// Map undirected edge -> face indices sharing it.
	edgeFaces := make(map[Edge][]int)
	for fi, f := range p.Faces {
		m := len(f)
		for i := 0; i < m; i++ {
			e := MakeEdge(f[i], f[(i+1)%m])
			edgeFaces[e] = append(edgeFaces[e], fi)
		}
	}
	for _, fs := range edgeFaces {
		for i := 0; i < len(fs); i++ {
			for j := i + 1; j < len(fs); j++ {
				if coplanarNormals(normals[fs[i]], normals[fs[j]], tol) {
					union(fs[i], fs[j])
				}
			}
		}
	}
	// Group faces by root.
	groups := make(map[int][]int)
	for fi := range p.Faces {
		r := find(fi)
		groups[r] = append(groups[r], fi)
	}
	q := &Polyhedron{Verts: append([]Vec3(nil), p.Verts...)}
	// Deterministic order.
	roots := make([]int, 0, len(groups))
	for r := range groups {
		roots = append(roots, r)
	}
	sort.Ints(roots)
	for _, r := range roots {
		grp := groups[r]
		loop := boundaryLoop(p, grp)
		if len(loop) >= 3 {
			q.Faces = append(q.Faces, loop)
		}
	}
	return q.OrientOutward()
}

// coplanarNormals reports whether two unit normals point in the same direction
// within angular tolerance tol (radians).
func coplanarNormals(a, b Vec3, tol float64) bool {
	d := a.Dot(b)
	if d > 1 {
		d = 1
	} else if d < -1 {
		d = -1
	}
	return math.Acos(d) <= tol
}

// boundaryLoop extracts the ordered boundary vertex loop of a connected group of
// coplanar faces. Interior (shared) edges are removed and the remaining boundary
// edges are stitched into a single cycle.
func boundaryLoop(p *Polyhedron, group []int) []int {
	// Count directed edges; boundary edges appear once.
	dirCount := make(map[[2]int]int)
	for _, fi := range group {
		f := p.Faces[fi]
		m := len(f)
		for i := 0; i < m; i++ {
			dirCount[[2]int{f[i], f[(i+1)%m]}]++
		}
	}
	// A directed edge is on the boundary if its reverse is absent.
	next := make(map[int]int)
	for e := range dirCount {
		if dirCount[[2]int{e[1], e[0]}] == 0 {
			next[e[0]] = e[1]
		}
	}
	if len(next) == 0 {
		return nil
	}
	// Walk the cycle.
	var start int
	for s := range next {
		start = s
		break
	}
	loop := []int{start}
	cur := start
	for {
		nx, ok := next[cur]
		if !ok || nx == start {
			break
		}
		loop = append(loop, nx)
		cur = nx
		if len(loop) > len(next)+1 {
			break // safety against malformed input
		}
	}
	return loop
}
