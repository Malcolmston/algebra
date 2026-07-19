package meshgen

import "math"

// TriMinAngle returns the smallest interior angle, in radians, of the triangle
// with vertices a, b, c. A degenerate triangle returns 0.
func TriMinAngle(a, b, c Vec2) float64 {
	angs := TriAngles(a, b, c)
	return math.Min(angs[0], math.Min(angs[1], angs[2]))
}

// TriMaxAngle returns the largest interior angle, in radians, of the triangle.
func TriMaxAngle(a, b, c Vec2) float64 {
	angs := TriAngles(a, b, c)
	return math.Max(angs[0], math.Max(angs[1], angs[2]))
}

// TriAngles returns the three interior angles of the triangle at vertices a, b
// and c respectively, in radians. A degenerate triangle returns zeros.
func TriAngles(a, b, c Vec2) [3]float64 {
	la := b.Distance(c)
	lb := a.Distance(c)
	lc := a.Distance(b)
	if la == 0 || lb == 0 || lc == 0 {
		return [3]float64{0, 0, 0}
	}
	angA := lawOfCosinesAngle(lb, lc, la)
	angB := lawOfCosinesAngle(la, lc, lb)
	angC := lawOfCosinesAngle(la, lb, lc)
	return [3]float64{angA, angB, angC}
}

// lawOfCosinesAngle returns the angle opposite side of length opp in a triangle
// with the other two sides of length x and y.
func lawOfCosinesAngle(x, y, opp float64) float64 {
	den := 2 * x * y
	if den == 0 {
		return 0
	}
	cos := (x*x + y*y - opp*opp) / den
	if cos > 1 {
		cos = 1
	} else if cos < -1 {
		cos = -1
	}
	return math.Acos(cos)
}

// locateReal reports whether p lies inside some real (non-super) triangle,
// returning that triangle index and true.
func (t *Triangulation) locateReal(p Vec2) (int, bool) {
	for i, tr := range t.tris {
		if t.triHasSuper(tr) {
			continue
		}
		ct := t.orientedCCW(tr)
		if PointInTriangle(p, t.pts[ct.A], t.pts[ct.B], t.pts[ct.C], 1e-12) {
			return i, true
		}
	}
	return -1, false
}

// worstTriangle returns the real triangle whose minimum angle is smallest and
// below minAngleRad, or ok=false when none is below the threshold.
func (t *Triangulation) worstTriangle(minAngleRad float64) (int, bool) {
	best := -1
	bestAngle := minAngleRad
	for i, tr := range t.tris {
		if t.triHasSuper(tr) {
			continue
		}
		ct := t.orientedCCW(tr)
		ang := TriMinAngle(t.pts[ct.A], t.pts[ct.B], t.pts[ct.C])
		if ang < bestAngle {
			bestAngle = ang
			best = i
		}
	}
	return best, best >= 0
}

// realBoundaryEdges returns the edges incident to exactly one real (non-super)
// triangle, i.e. the boundary of the meshed region.
func (t *Triangulation) realBoundaryEdges() []Edge {
	count := make(map[Edge]int)
	for _, tr := range t.tris {
		if t.triHasSuper(tr) {
			continue
		}
		for _, e := range tr.EdgesOf() {
			count[e]++
		}
	}
	var out []Edge
	for e, n := range count {
		if n == 1 {
			out = append(out, e)
		}
	}
	return out
}

// encroachedBoundary returns a boundary edge whose diametral circle strictly
// contains p, if any.
func (t *Triangulation) encroachedBoundary(p Vec2) (Edge, bool) {
	for _, e := range t.realBoundaryEdges() {
		a, b := t.pts[e.U], t.pts[e.V]
		center := a.Midpoint(b)
		radius := a.Distance(b) / 2
		if p.Distance(center) < radius-1e-12 {
			return e, true
		}
	}
	return Edge{}, false
}

// longestBoundaryEdgeOf returns the longest real-boundary edge of triangle ct.
func (t *Triangulation) longestBoundaryEdgeOf(ct Tri) (Edge, bool) {
	bset := make(map[Edge]struct{})
	for _, e := range t.realBoundaryEdges() {
		bset[e] = struct{}{}
	}
	best := Edge{}
	bestLen := -1.0
	found := false
	for _, e := range ct.EdgesOf() {
		if _, ok := bset[e]; !ok {
			continue
		}
		if l := t.pts[e.U].Distance(t.pts[e.V]); l > bestLen {
			bestLen = l
			best = e
			found = true
		}
	}
	return best, found
}

func (t *Triangulation) pointExists(p Vec2) bool {
	for _, q := range t.pts {
		if q.ApproxEqual(p, 1e-12) {
			return true
		}
	}
	return false
}

// globalWorstAngle returns the smallest interior angle (radians) over all real
// triangles, or +Inf when there are none.
func (t *Triangulation) globalWorstAngle() float64 {
	worst := math.Inf(1)
	for _, tr := range t.tris {
		if t.triHasSuper(tr) {
			continue
		}
		ct := t.orientedCCW(tr)
		if a := TriMinAngle(t.pts[ct.A], t.pts[ct.B], t.pts[ct.C]); a < worst {
			worst = a
		}
	}
	return worst
}

// worstTriangleExcluding returns the worst real triangle below minAngleRad
// whose centroid key is not in skip.
func (t *Triangulation) worstTriangleExcluding(minAngleRad float64, skip map[[2]int64]bool) (int, bool) {
	best := -1
	bestAngle := minAngleRad
	for i, tr := range t.tris {
		if t.triHasSuper(tr) {
			continue
		}
		ct := t.orientedCCW(tr)
		if skip[centroidKey(t.pts[ct.A], t.pts[ct.B], t.pts[ct.C])] {
			continue
		}
		ang := TriMinAngle(t.pts[ct.A], t.pts[ct.B], t.pts[ct.C])
		if ang < bestAngle {
			bestAngle = ang
			best = i
		}
	}
	return best, best >= 0
}

func centroidKey(a, b, c Vec2) [2]int64 {
	cx := (a.X + b.X + c.X) / 3
	cy := (a.Y + b.Y + c.Y) / 3
	return [2]int64{int64(cx * 1e6), int64(cy * 1e6)}
}

// Refine performs Ruppert-style Delaunay refinement guarded for robustness: it
// repeatedly picks the worst real triangle whose minimum angle is below
// minAngleDeg and inserts a Steiner point at its circumcenter, or at the
// midpoint of an encroached or exterior boundary edge when the circumcenter is
// unsuitable. Each candidate insertion is accepted only if it does not reduce
// the mesh's overall minimum angle; otherwise it is rolled back and that
// triangle is set aside. This guarantees the mesh quality never degrades.
// At most maxSteiner points are inserted. It returns the number added.
//
// A minimum angle up to about 20 degrees is achievable for well-posed inputs
// whose boundary contains no small angles; the guard may plateau earlier on
// hard domains such as long thin strips.
func (t *Triangulation) Refine(minAngleDeg float64, maxSteiner int) int {
	minAngle := minAngleDeg * math.Pi / 180
	added := 0
	skip := make(map[[2]int64]bool)
	for added < maxSteiner {
		bi, ok := t.worstTriangleExcluding(minAngle, skip)
		if !ok {
			break
		}
		ct := t.orientedCCW(t.tris[bi])
		a, b, c := t.pts[ct.A], t.pts[ct.B], t.pts[ct.C]
		cc, err := Circumcenter(a, b, c)
		if err != nil {
			skip[centroidKey(a, b, c)] = true
			continue
		}
		var target Vec2
		if _, inside := t.locateReal(cc); inside {
			if e, enc := t.encroachedBoundary(cc); enc {
				target = t.pts[e.U].Midpoint(t.pts[e.V])
			} else {
				target = cc
			}
		} else if e, ok := t.longestBoundaryEdgeOf(ct); ok {
			target = t.pts[e.U].Midpoint(t.pts[e.V])
		} else {
			skip[centroidKey(a, b, c)] = true
			continue
		}
		if t.pointExists(target) {
			skip[centroidKey(a, b, c)] = true
			continue
		}
		before := t.globalWorstAngle()
		savedPts := append([]Vec2(nil), t.pts...)
		savedTris := append([]Tri(nil), t.tris...)
		t.pts = append(t.pts, target)
		t.insertPoint(len(t.pts) - 1)
		if t.globalWorstAngle() <= before+1e-12 {
			// Insertion did not improve overall quality: roll back.
			t.pts = savedPts
			t.tris = savedTris
			skip[centroidKey(a, b, c)] = true
			continue
		}
		added++
		skip = make(map[[2]int64]bool)
	}
	return added
}

// RefineToMinAngle triangulates the points and applies Delaunay refinement to
// the given minimum angle, returning the refined mesh.
func RefineToMinAngle(points []Vec2, minAngleDeg float64, maxSteiner int) (*Mesh, error) {
	t, err := Triangulate(points)
	if err != nil {
		return nil, err
	}
	t.Refine(minAngleDeg, maxSteiner)
	return t.Mesh(), nil
}
