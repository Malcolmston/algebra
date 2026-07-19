package meshgen

import "errors"

// ErrConstraintFailed is returned when a constraint edge cannot be recovered,
// for example because another input vertex lies on the segment.
var ErrConstraintFailed = errors.New("meshgen: unable to recover constraint edge")

// hasEdge reports whether the undirected edge (a,b) is present in some triangle.
func (t *Triangulation) hasEdge(a, b int) bool {
	e := NewEdge(a, b)
	for _, tr := range t.tris {
		for _, te := range tr.EdgesOf() {
			if te == e {
				return true
			}
		}
	}
	return false
}

// apexes returns the two apex vertices of the (up to two) triangles sharing the
// edge (u,v), and reports whether exactly two triangles share it.
func (t *Triangulation) apexes(u, v int) (w1, w2 int, ok bool) {
	e := NewEdge(u, v)
	found := 0
	for _, tr := range t.tris {
		has := false
		for _, te := range tr.EdgesOf() {
			if te == e {
				has = true
				break
			}
		}
		if !has {
			continue
		}
		w, o := tr.Opposite(u, v)
		if !o {
			continue
		}
		switch found {
		case 0:
			w1 = w
		case 1:
			w2 = w
		default:
			return 0, 0, false
		}
		found++
	}
	return w1, w2, found == 2
}

// doFlip replaces the two triangles sharing edge (u,v) - with apexes w1,w2 - by
// the two triangles sharing the flipped diagonal (w1,w2).
func (t *Triangulation) doFlip(u, v, w1, w2 int) {
	e := NewEdge(u, v)
	kept := t.tris[:0:0]
	for _, tr := range t.tris {
		isShared := false
		for _, te := range tr.EdgesOf() {
			if te == e {
				isShared = true
				break
			}
		}
		if !isShared {
			kept = append(kept, tr)
		}
	}
	kept = append(kept, Tri{u, w1, w2}, Tri{w1, v, w2})
	t.tris = kept
}

// crossingEdges returns the edges whose segment properly crosses segment (a,b).
func (t *Triangulation) crossingEdges(a, b int) []Edge {
	pa, pb := t.pts[a], t.pts[b]
	set := make(map[Edge]struct{})
	var out []Edge
	for _, tr := range t.tris {
		for _, e := range tr.EdgesOf() {
			if e.U == a || e.U == b || e.V == a || e.V == b {
				continue
			}
			if _, seen := set[e]; seen {
				continue
			}
			if SegmentsProperlyIntersect(pa, pb, t.pts[e.U], t.pts[e.V]) {
				set[e] = struct{}{}
				out = append(out, e)
			}
		}
	}
	return out
}

// InsertConstraint forces the edge between input vertices a and b to appear in
// the triangulation by flipping the diagonals that cross it. The indices refer
// to the original input points. It returns ErrConstraintFailed if the edge
// cannot be recovered.
func (t *Triangulation) InsertConstraint(a, b int) error {
	if a == b || a < 0 || b < 0 || a >= len(t.pts) || b >= len(t.pts) {
		return ErrDegenerate
	}
	const maxIter = 1 << 20
	for iter := 0; iter < maxIter; iter++ {
		if t.hasEdge(a, b) {
			return nil
		}
		crossings := t.crossingEdges(a, b)
		if len(crossings) == 0 {
			return ErrConstraintFailed
		}
		flipped := false
		for _, e := range crossings {
			w1, w2, ok := t.apexes(e.U, e.V)
			if !ok {
				continue
			}
			if SegmentsProperlyIntersect(t.pts[e.U], t.pts[e.V], t.pts[w1], t.pts[w2]) {
				t.doFlip(e.U, e.V, w1, w2)
				flipped = true
				break
			}
		}
		if !flipped {
			return ErrConstraintFailed
		}
	}
	return ErrConstraintFailed
}

// TriangulateConstrained returns the constrained Delaunay triangulation of the
// points in which every segment (a pair of input point indices) is guaranteed
// to appear as an edge.
func TriangulateConstrained(points []Vec2, segments [][2]int) (*Triangulation, error) {
	t, err := Triangulate(points)
	if err != nil {
		return nil, err
	}
	for _, s := range segments {
		if err := t.InsertConstraint(s[0], s[1]); err != nil {
			return nil, err
		}
	}
	return t, nil
}
