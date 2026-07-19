package meshgen

import (
	"errors"
	"math"
	"sort"
)

// ErrNotEnoughPoints is returned when fewer than three non-degenerate points
// are supplied to a triangulation routine.
var ErrNotEnoughPoints = errors.New("meshgen: need at least three non-collinear points")

// Triangulation is a mutable planar triangulation supporting incremental
// point insertion (Bowyer-Watson), constraint recovery and refinement. The
// points at indices [0, NumInput) are the caller's input points; internal
// helper vertices (super-triangle corners and refinement points) follow.
type Triangulation struct {
	pts      []Vec2
	tris     []Tri
	numInput int
	super    [3]int // indices of super-triangle corners in pts
}

// Triangulate returns the Delaunay triangulation of the given planar points
// using the incremental Bowyer-Watson algorithm. Duplicate points are ignored.
func Triangulate(points []Vec2) (*Triangulation, error) {
	if len(points) < 3 {
		return nil, ErrNotEnoughPoints
	}
	t := &Triangulation{numInput: len(points)}
	t.pts = append(t.pts, points...)
	if err := t.buildSuper(points); err != nil {
		return nil, err
	}
	order := make([]int, len(points))
	for i := range order {
		order[i] = i
	}
	for _, i := range order {
		if t.duplicateOf(i) >= 0 {
			continue
		}
		t.insertPoint(i)
	}
	if !t.hasRealTriangle() {
		return nil, ErrNotEnoughPoints
	}
	return t, nil
}

// buildSuper appends three super-triangle corners enclosing all points.
func (t *Triangulation) buildSuper(points []Vec2) error {
	min, max := BoundingBox2(points)
	if !min.IsFinite() || !max.IsFinite() {
		return ErrDegenerate
	}
	dx := max.X - min.X
	dy := max.Y - min.Y
	d := math.Max(dx, dy)
	if d == 0 {
		d = 1
	}
	mid := min.Add(max).Div(2)
	m := 20 * d
	s0 := Vec2{mid.X - m, mid.Y - d}
	s1 := Vec2{mid.X, mid.Y + m}
	s2 := Vec2{mid.X + m, mid.Y - d}
	base := len(t.pts)
	t.pts = append(t.pts, s0, s1, s2)
	t.super = [3]int{base, base + 1, base + 2}
	t.tris = []Tri{{base, base + 1, base + 2}}
	return nil
}

func (t *Triangulation) duplicateOf(i int) int {
	for j := 0; j < i; j++ {
		if t.pts[j].Equal(t.pts[i]) {
			return j
		}
	}
	return -1
}

// insertPoint inserts vertex index p into the current triangulation.
func (t *Triangulation) insertPoint(p int) {
	pt := t.pts[p]
	var bad []int
	for i, tr := range t.tris {
		a := t.orientedCCW(tr)
		if InCircle(t.pts[a.A], t.pts[a.B], t.pts[a.C], pt) > 0 {
			bad = append(bad, i)
		}
	}
	if len(bad) == 0 {
		return
	}
	// Collect boundary of the cavity: edges belonging to exactly one bad tri.
	type de struct{ a, b int } // directed edge
	count := make(map[Edge]int)
	dirEdges := make([]de, 0, len(bad)*3)
	badSet := make(map[int]bool, len(bad))
	for _, bi := range bad {
		badSet[bi] = true
		tr := t.tris[bi]
		vs := [3]int{tr.A, tr.B, tr.C}
		for k := 0; k < 3; k++ {
			a, b := vs[k], vs[(k+1)%3]
			dirEdges = append(dirEdges, de{a, b})
			count[NewEdge(a, b)]++
		}
	}
	// Remove bad triangles.
	kept := t.tris[:0:0]
	for i, tr := range t.tris {
		if !badSet[i] {
			kept = append(kept, tr)
		}
	}
	t.tris = kept
	// Retriangulate the cavity using boundary directed edges.
	for _, e := range dirEdges {
		if count[NewEdge(e.a, e.b)] == 1 {
			t.tris = append(t.tris, Tri{e.a, e.b, p})
		}
	}
}

// orientedCCW returns the triangle with vertices ordered counterclockwise.
func (t *Triangulation) orientedCCW(tr Tri) Tri {
	if Orient2D(t.pts[tr.A], t.pts[tr.B], t.pts[tr.C]) < 0 {
		return Tri{tr.A, tr.C, tr.B}
	}
	return tr
}

func (t *Triangulation) isSuper(v int) bool {
	return v == t.super[0] || v == t.super[1] || v == t.super[2]
}

func (t *Triangulation) triHasSuper(tr Tri) bool {
	return t.isSuper(tr.A) || t.isSuper(tr.B) || t.isSuper(tr.C)
}

func (t *Triangulation) hasRealTriangle() bool {
	for _, tr := range t.tris {
		if !t.triHasSuper(tr) {
			return true
		}
	}
	return false
}

// Points returns the input points (excluding internal helper vertices). The
// slice may be longer than the original input if refinement added Steiner
// points; those are included and appear after the input points.
func (t *Triangulation) Points() []Vec2 {
	out := make([]Vec2, 0, len(t.pts))
	for i, p := range t.pts {
		if !t.isSuper(i) {
			out = append(out, p)
		}
	}
	return out
}

// NumInput returns the number of original input points.
func (t *Triangulation) NumInput() int { return t.numInput }

// Mesh returns the triangulation as a Mesh, dropping triangles that touch the
// super-triangle and renumbering vertices to a compact 0-based index set.
func (t *Triangulation) Mesh() *Mesh {
	remap := make(map[int]int)
	var verts []Vec2
	idx := func(v int) int {
		if r, ok := remap[v]; ok {
			return r
		}
		r := len(verts)
		remap[v] = r
		verts = append(verts, t.pts[v])
		return r
	}
	var tris []Tri
	for _, tr := range t.tris {
		if t.triHasSuper(tr) {
			continue
		}
		ct := t.orientedCCW(tr)
		tris = append(tris, Tri{idx(ct.A), idx(ct.B), idx(ct.C)})
	}
	return &Mesh{Vertices: verts, Triangles: tris}
}

// NumTriangles returns the number of real (non-super) triangles.
func (t *Triangulation) NumTriangles() int {
	n := 0
	for _, tr := range t.tris {
		if !t.triHasSuper(tr) {
			n++
		}
	}
	return n
}

// IsDelaunay reports whether every real triangle satisfies the empty-circle
// property with respect to all input vertices, within tolerance eps.
func (t *Triangulation) IsDelaunay(eps float64) bool {
	for _, tr := range t.tris {
		if t.triHasSuper(tr) {
			continue
		}
		ct := t.orientedCCW(tr)
		for i := 0; i < len(t.pts); i++ {
			if t.isSuper(i) || i == ct.A || i == ct.B || i == ct.C {
				continue
			}
			if InCircle(t.pts[ct.A], t.pts[ct.B], t.pts[ct.C], t.pts[i]) > eps {
				return false
			}
		}
	}
	return true
}

// DelaunayMesh is a convenience wrapper returning the Delaunay triangulation of
// the points directly as a Mesh.
func DelaunayMesh(points []Vec2) (*Mesh, error) {
	t, err := Triangulate(points)
	if err != nil {
		return nil, err
	}
	return t.Mesh(), nil
}

// ConvexHull returns the counterclockwise convex hull of the points using
// Andrew's monotone chain algorithm. Duplicate and collinear interior points
// are removed.
func ConvexHull(points []Vec2) []Vec2 {
	n := len(points)
	if n < 3 {
		out := make([]Vec2, n)
		copy(out, points)
		return out
	}
	pts := make([]Vec2, n)
	copy(pts, points)
	sortVec2(pts)
	// deduplicate
	dedup := pts[:1]
	for _, p := range pts[1:] {
		if !p.Equal(dedup[len(dedup)-1]) {
			dedup = append(dedup, p)
		}
	}
	pts = dedup
	if len(pts) < 3 {
		return pts
	}
	var hull []Vec2
	// lower
	for _, p := range pts {
		for len(hull) >= 2 && Orient2D(hull[len(hull)-2], hull[len(hull)-1], p) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, p)
	}
	// upper
	lower := len(hull) + 1
	for i := len(pts) - 2; i >= 0; i-- {
		p := pts[i]
		for len(hull) >= lower && Orient2D(hull[len(hull)-2], hull[len(hull)-1], p) <= 0 {
			hull = hull[:len(hull)-1]
		}
		hull = append(hull, p)
	}
	return hull[:len(hull)-1]
}

func sortVec2(pts []Vec2) {
	sort.Slice(pts, func(i, j int) bool {
		if pts[i].X != pts[j].X {
			return pts[i].X < pts[j].X
		}
		return pts[i].Y < pts[j].Y
	})
}
