package designs

import "errors"

// ProjectivePlane is the finite Desarguesian projective plane PG(2,q) of
// prime-power order q. It has q^2+q+1 points and the same number of lines; every
// line contains q+1 points, every point lies on q+1 lines, two distinct points
// determine a unique line and two distinct lines meet in a unique point.
type ProjectivePlane struct {
	q      int
	f      *GaloisField
	points [][3]int
	lines  [][3]int
	pIndex map[[3]int]int
	lIndex map[[3]int]int
}

// normalizeTriple scales a non-zero triple so its first non-zero coordinate is
// 1, giving a canonical representative of a projective point or line.
func (pl *ProjectivePlane) normalizeTriple(t [3]int) ([3]int, bool) {
	for i := 0; i < 3; i++ {
		if t[i] != 0 {
			inv, _ := pl.f.Inv(t[i])
			var out [3]int
			for j := 0; j < 3; j++ {
				out[j] = pl.f.Mul(t[j], inv)
			}
			return out, true
		}
	}
	return [3]int{}, false
}

// enumerateTriples lists the canonical representatives of the projective points
// of PG(2,q): (x,y,1), then (x,1,0), then (1,0,0).
func (pl *ProjectivePlane) enumerateTriples() [][3]int {
	q := pl.q
	var out [][3]int
	for x := 0; x < q; x++ {
		for y := 0; y < q; y++ {
			out = append(out, [3]int{x, y, 1})
		}
	}
	for x := 0; x < q; x++ {
		out = append(out, [3]int{x, 1, 0})
	}
	out = append(out, [3]int{1, 0, 0})
	return out
}

// NewProjectivePlane constructs PG(2,q) for a prime power q>=2. It reports an
// error when q is not a prime power.
func NewProjectivePlane(q int) (*ProjectivePlane, error) {
	if !IsPrimePower(q) {
		return nil, errors.New("designs: order must be a prime power")
	}
	f, err := NewGaloisField(q)
	if err != nil {
		return nil, err
	}
	pl := &ProjectivePlane{q: q, f: f}
	pl.points = pl.enumerateTriples()
	pl.lines = pl.enumerateTriples() // points and lines share coordinates
	pl.pIndex = make(map[[3]int]int, len(pl.points))
	pl.lIndex = make(map[[3]int]int, len(pl.lines))
	for i, p := range pl.points {
		pl.pIndex[p] = i
	}
	for j, l := range pl.lines {
		pl.lIndex[l] = j
	}
	return pl, nil
}

// Order returns the order q of the plane.
func (pl *ProjectivePlane) Order() int { return pl.q }

// NumPoints returns the number of points q^2+q+1.
func (pl *ProjectivePlane) NumPoints() int { return len(pl.points) }

// NumLines returns the number of lines q^2+q+1.
func (pl *ProjectivePlane) NumLines() int { return len(pl.lines) }

// Point returns the homogeneous coordinates of point index i.
func (pl *ProjectivePlane) Point(i int) [3]int { return pl.points[i] }

// Line returns the homogeneous coordinates of line index j.
func (pl *ProjectivePlane) Line(j int) [3]int { return pl.lines[j] }

// dot evaluates the field inner product a·b over GF(q).
func (pl *ProjectivePlane) dot(a, b [3]int) int {
	s := 0
	for i := 0; i < 3; i++ {
		s = pl.f.Add(s, pl.f.Mul(a[i], b[i]))
	}
	return s
}

// cross returns the field cross product a×b over GF(q).
func (pl *ProjectivePlane) cross(a, b [3]int) [3]int {
	sub := pl.f.Sub
	mul := pl.f.Mul
	return [3]int{
		sub(mul(a[1], b[2]), mul(a[2], b[1])),
		sub(mul(a[2], b[0]), mul(a[0], b[2])),
		sub(mul(a[0], b[1]), mul(a[1], b[0])),
	}
}

// IsIncident reports whether point index p lies on line index l, i.e. the field
// inner product of their coordinates is zero.
func (pl *ProjectivePlane) IsIncident(p, l int) bool {
	return pl.dot(pl.points[p], pl.lines[l]) == 0
}

// PointsOnLine returns the indices of the q+1 points incident with line l.
func (pl *ProjectivePlane) PointsOnLine(l int) []int {
	var out []int
	for p := range pl.points {
		if pl.IsIncident(p, l) {
			out = append(out, p)
		}
	}
	return out
}

// LinesThroughPoint returns the indices of the q+1 lines incident with point p.
func (pl *ProjectivePlane) LinesThroughPoint(p int) []int {
	var out []int
	for l := range pl.lines {
		if pl.IsIncident(p, l) {
			out = append(out, l)
		}
	}
	return out
}

// LineThroughPoints returns the index of the unique line joining two distinct
// points. It reports an error when the point indices coincide.
func (pl *ProjectivePlane) LineThroughPoints(p1, p2 int) (int, error) {
	if p1 == p2 {
		return 0, errors.New("designs: points must be distinct")
	}
	c := pl.cross(pl.points[p1], pl.points[p2])
	n, ok := pl.normalizeTriple(c)
	if !ok {
		return 0, errors.New("designs: degenerate join")
	}
	return pl.lIndex[n], nil
}

// MeetOfLines returns the index of the unique point common to two distinct
// lines. It reports an error when the line indices coincide.
func (pl *ProjectivePlane) MeetOfLines(l1, l2 int) (int, error) {
	if l1 == l2 {
		return 0, errors.New("designs: lines must be distinct")
	}
	c := pl.cross(pl.lines[l1], pl.lines[l2])
	n, ok := pl.normalizeTriple(c)
	if !ok {
		return 0, errors.New("designs: degenerate meet")
	}
	return pl.pIndex[n], nil
}

// IncidenceDesign returns the points-and-lines incidence structure as a Design:
// one point per plane point and one block per line listing the points on it.
// This is a symmetric 2-(q^2+q+1, q+1, 1) design.
func (pl *ProjectivePlane) IncidenceDesign() *Design {
	blocks := make([][]int, len(pl.lines))
	for l := range pl.lines {
		blocks[l] = pl.PointsOnLine(l)
	}
	d, _ := NewDesign(len(pl.points), blocks)
	return d
}

// IsProjectivePlane verifies the incidence axioms: the plane has q^2+q+1 points
// and lines, every line has q+1 points, and every pair of points lies on
// exactly one common line.
func (pl *ProjectivePlane) IsProjectivePlane() bool {
	d := pl.IncidenceDesign()
	params, err := d.Parameters()
	if err != nil {
		return false
	}
	return params.K == pl.q+1 && params.Lambda == 1 && params.V == pl.q*pl.q+pl.q+1
}

// FanoPlane returns the Fano plane PG(2,2), the unique projective plane of
// order 2 with 7 points and 7 lines.
func FanoPlane() *ProjectivePlane {
	pl, _ := NewProjectivePlane(2)
	return pl
}

// AffinePlane is the finite affine plane AG(2,q) of prime-power order q, with
// q^2 points and q^2+q lines partitioned into q+1 parallel classes of q lines
// each. Every line has q points and two points determine a unique line.
type AffinePlane struct {
	q         int
	f         *GaloisField
	lineSets  [][]int // each line as sorted point indices
	parallels [][]int // parallel classes as line indices
}

// NewAffinePlane constructs AG(2,q) for a prime power q>=2. It reports an error
// when q is not a prime power.
func NewAffinePlane(q int) (*AffinePlane, error) {
	if !IsPrimePower(q) {
		return nil, errors.New("designs: order must be a prime power")
	}
	f, err := NewGaloisField(q)
	if err != nil {
		return nil, err
	}
	ap := &AffinePlane{q: q, f: f}
	idx := func(x, y int) int { return x*q + y }
	// Lines with finite slope m and intercept b: y = m*x + b.
	for m := 0; m < q; m++ {
		class := []int{}
		for b := 0; b < q; b++ {
			line := make([]int, 0, q)
			for x := 0; x < q; x++ {
				y := f.Add(f.Mul(m, x), b)
				line = append(line, idx(x, y))
			}
			class = append(class, len(ap.lineSets))
			ap.lineSets = append(ap.lineSets, line)
		}
		ap.parallels = append(ap.parallels, class)
	}
	// Vertical lines x = c.
	vclass := []int{}
	for c := 0; c < q; c++ {
		line := make([]int, 0, q)
		for y := 0; y < q; y++ {
			line = append(line, idx(c, y))
		}
		vclass = append(vclass, len(ap.lineSets))
		ap.lineSets = append(ap.lineSets, line)
	}
	ap.parallels = append(ap.parallels, vclass)
	return ap, nil
}

// Order returns the order q of the affine plane.
func (ap *AffinePlane) Order() int { return ap.q }

// NumPoints returns the number of points q^2.
func (ap *AffinePlane) NumPoints() int { return ap.q * ap.q }

// NumLines returns the number of lines q^2+q.
func (ap *AffinePlane) NumLines() int { return len(ap.lineSets) }

// NumParallelClasses returns the number of parallel classes q+1.
func (ap *AffinePlane) NumParallelClasses() int { return len(ap.parallels) }

// PointsOnLine returns the indices of the q points on line l.
func (ap *AffinePlane) PointsOnLine(l int) []int {
	return append([]int(nil), ap.lineSets[l]...)
}

// IsIncident reports whether point index p lies on line index l.
func (ap *AffinePlane) IsIncident(p, l int) bool {
	for _, x := range ap.lineSets[l] {
		if x == p {
			return true
		}
	}
	return false
}

// ParallelClasses returns the parallel classes, each as a slice of line
// indices; lines in the same class are pairwise disjoint and together cover
// every point exactly once.
func (ap *AffinePlane) ParallelClasses() [][]int {
	out := make([][]int, len(ap.parallels))
	for i, c := range ap.parallels {
		out[i] = append([]int(nil), c...)
	}
	return out
}

// AreParallel reports whether two distinct lines belong to the same parallel
// class (are disjoint). A line is considered parallel to itself.
func (ap *AffinePlane) AreParallel(l1, l2 int) bool {
	if l1 == l2 {
		return true
	}
	for _, c := range ap.parallels {
		in1, in2 := false, false
		for _, x := range c {
			if x == l1 {
				in1 = true
			}
			if x == l2 {
				in2 = true
			}
		}
		if in1 && in2 {
			return true
		}
	}
	return false
}

// IncidenceDesign returns the affine plane as a resolvable 2-(q^2, q, 1) design
// whose blocks are the lines.
func (ap *AffinePlane) IncidenceDesign() *Design {
	blocks := make([][]int, len(ap.lineSets))
	for l := range ap.lineSets {
		blocks[l] = append([]int(nil), ap.lineSets[l]...)
	}
	d, _ := NewDesign(ap.q*ap.q, blocks)
	return d
}

// IsAffinePlane verifies that the incidence structure is a 2-(q^2,q,1) design
// with q+1 parallel classes.
func (ap *AffinePlane) IsAffinePlane() bool {
	d := ap.IncidenceDesign()
	params, err := d.Parameters()
	if err != nil {
		return false
	}
	return params.K == ap.q && params.Lambda == 1 && params.V == ap.q*ap.q &&
		len(ap.parallels) == ap.q+1
}
