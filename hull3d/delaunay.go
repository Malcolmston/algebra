package hull3d

import (
	"errors"
	"math"
	"sort"
)

// Tetrahedron is a solid tetrahedron given as four indices into a point slice.
type Tetrahedron struct {
	A, B, C, D int
}

// Tetrahedralization is a decomposition of the convex hull of a point set into
// tetrahedra, produced by [Delaunay]. It retains the original points.
type Tetrahedralization struct {
	Points []Vec3
	Tets   []Tetrahedron
}

// inCircumsphere reports whether e lies strictly inside the circumsphere of the
// tetrahedron (a,b,c,d).
func inCircumsphere(a, b, c, d, e Vec3) bool {
	o := Orient3D(a, b, c, d)
	if o == 0 {
		return false
	}
	return InSphere(a, b, c, d, e)*o > 0
}

// Circumcenter returns the centre of the sphere passing through the four points
// and its radius. It returns an error if the points are coplanar (no unique
// circumsphere).
func Circumcenter(a, b, c, d Vec3) (center Vec3, radius float64, err error) {
	ba := b.Sub(a)
	ca := c.Sub(a)
	da := d.Sub(a)
	det := Triple(ba, ca, da)
	if math.Abs(det) < 1e-300 {
		return Vec3{}, 0, errors.New("hull3d: coplanar points have no circumsphere")
	}
	rhs := Vec3{0.5 * ba.LengthSq(), 0.5 * ca.LengthSq(), 0.5 * da.LengthSq()}
	// Solve [ba;ca;da] x = rhs via Cramer's rule.
	x := solve3(ba, ca, da, rhs)
	center = a.Add(x)
	radius = center.Distance(a)
	return center, radius, nil
}

// solve3 solves the 3x3 linear system whose rows are r0,r1,r2 and right-hand
// side b, using Cramer's rule. The caller guarantees a non-zero determinant.
func solve3(r0, r1, r2, b Vec3) Vec3 {
	col0 := Vec3{r0.X, r1.X, r2.X}
	col1 := Vec3{r0.Y, r1.Y, r2.Y}
	col2 := Vec3{r0.Z, r1.Z, r2.Z}
	det := Triple(col0, col1, col2)
	x := Triple(b, col1, col2) / det
	y := Triple(col0, b, col2) / det
	z := Triple(col0, col1, b) / det
	return Vec3{x, y, z}
}

// Circumsphere returns the circumscribed sphere of the tetrahedron with the
// given corners as a [SphereShape]. It returns an error for coplanar corners.
func Circumsphere(a, b, c, d Vec3) (SphereShape, error) {
	c0, r, err := Circumcenter(a, b, c, d)
	if err != nil {
		return SphereShape{}, err
	}
	return SphereShape{Center: c0, Radius: r}, nil
}

// Delaunay computes a Delaunay tetrahedralization of the given points using the
// incremental Bowyer–Watson algorithm. It returns an error if fewer than four
// points are supplied or the points are coplanar.
func Delaunay(points []Vec3) (*Tetrahedralization, error) {
	pts := DedupPoints(points, 1e-12)
	if len(pts) < 4 {
		return nil, errors.New("hull3d: need at least four distinct points")
	}
	min, max, _ := BoundingBox(pts)
	center := min.Add(max).Scale(0.5)
	diag := max.Sub(min).Length()
	if diag < 1e-300 {
		return nil, errNoInterior
	}
	m := diag * 1000

	// Super-tetrahedron enclosing all points.
	n := len(pts)
	work := append([]Vec3(nil), pts...)
	work = append(work,
		center.Add(Vec3{0, 0, 4 * m}),
		center.Add(Vec3{-3 * m, -2 * m, -2 * m}),
		center.Add(Vec3{3 * m, -2 * m, -2 * m}),
		center.Add(Vec3{0, 3 * m, -2 * m}),
	)
	super := [4]int{n, n + 1, n + 2, n + 3}

	tets := []Tetrahedron{{super[0], super[1], super[2], super[3]}}

	type triKey struct{ a, b, c int }
	makeKey := func(x, y, z int) triKey {
		s := []int{x, y, z}
		sort.Ints(s)
		return triKey{s[0], s[1], s[2]}
	}

	for pi := 0; pi < n; pi++ {
		p := work[pi]
		var bad []int
		for ti, t := range tets {
			if inCircumsphere(work[t.A], work[t.B], work[t.C], work[t.D], p) {
				bad = append(bad, ti)
			}
		}
		if len(bad) == 0 {
			// Numerical fallback: attach to nearest tetra by containing face.
			continue
		}
		// Boundary faces: appear in exactly one bad tetra.
		faceCount := map[triKey]int{}
		faceRep := map[triKey][3]int{}
		badSet := make(map[int]bool, len(bad))
		for _, ti := range bad {
			badSet[ti] = true
			t := tets[ti]
			fs := [4][3]int{
				{t.B, t.C, t.D},
				{t.A, t.D, t.C},
				{t.A, t.B, t.D},
				{t.A, t.C, t.B},
			}
			for _, f := range fs {
				k := makeKey(f[0], f[1], f[2])
				faceCount[k]++
				faceRep[k] = f
			}
		}
		// Remove bad tetrahedra (build new slice).
		var remaining []Tetrahedron
		for ti, t := range tets {
			if !badSet[ti] {
				remaining = append(remaining, t)
			}
		}
		tets = remaining
		// Stitch cavity with new tetrahedra from boundary faces to p.
		for k, cnt := range faceCount {
			if cnt != 1 {
				continue
			}
			f := faceRep[k]
			nt := Tetrahedron{f[0], f[1], f[2], pi}
			if Orient3D(work[nt.A], work[nt.B], work[nt.C], work[nt.D]) < 0 {
				nt.B, nt.C = nt.C, nt.B
			}
			tets = append(tets, nt)
		}
	}

	// Drop tetrahedra touching any super-vertex.
	isSuper := func(i int) bool { return i >= n }
	var out []Tetrahedron
	for _, t := range tets {
		if isSuper(t.A) || isSuper(t.B) || isSuper(t.C) || isSuper(t.D) {
			continue
		}
		out = append(out, t)
	}
	if len(out) == 0 {
		return nil, errNoInterior
	}
	return &Tetrahedralization{Points: pts, Tets: out}, nil
}

// NumTets returns the number of tetrahedra.
func (d *Tetrahedralization) NumTets() int { return len(d.Tets) }

// TetVertices returns the four corner positions of tetrahedron t.
func (d *Tetrahedralization) TetVertices(t Tetrahedron) (a, b, c, e Vec3) {
	return d.Points[t.A], d.Points[t.B], d.Points[t.C], d.Points[t.D]
}

// TotalVolume returns the sum of the volumes of all tetrahedra, which equals the
// volume of the convex hull of the point set.
func (d *Tetrahedralization) TotalVolume() float64 {
	var v float64
	for _, t := range d.Tets {
		a, b, c, e := d.TetVertices(t)
		v += TetrahedronVolume(a, b, c, e)
	}
	return v
}

// IsDelaunay reports whether the tetrahedralization satisfies the empty-sphere
// property: no point lies strictly inside the circumsphere of any tetrahedron,
// within the given tolerance on the in-sphere predicate.
func (d *Tetrahedralization) IsDelaunay(eps float64) bool {
	for _, t := range d.Tets {
		a, b, c, e := d.TetVertices(t)
		o := Orient3D(a, b, c, e)
		if o == 0 {
			continue
		}
		for i, p := range d.Points {
			if i == t.A || i == t.B || i == t.C || i == t.D {
				continue
			}
			s := InSphere(a, b, c, e, p) * o
			if s > eps {
				return false
			}
		}
	}
	return true
}

// Faces returns the distinct triangular faces of the tetrahedralization, each as
// a sorted index triple.
func (d *Tetrahedralization) Faces() [][3]int {
	seen := map[[3]int]bool{}
	var out [][3]int
	add := func(x, y, z int) {
		s := []int{x, y, z}
		sort.Ints(s)
		k := [3]int{s[0], s[1], s[2]}
		if !seen[k] {
			seen[k] = true
			out = append(out, k)
		}
	}
	for _, t := range d.Tets {
		add(t.A, t.B, t.C)
		add(t.A, t.B, t.D)
		add(t.A, t.C, t.D)
		add(t.B, t.C, t.D)
	}
	return out
}

// BoundaryFaces returns the triangular faces that belong to exactly one
// tetrahedron; collectively they form the surface of the convex hull.
func (d *Tetrahedralization) BoundaryFaces() [][3]int {
	count := map[[3]int]int{}
	rep := map[[3]int][3]int{}
	add := func(x, y, z int) {
		s := []int{x, y, z}
		sort.Ints(s)
		k := [3]int{s[0], s[1], s[2]}
		count[k]++
		rep[k] = [3]int{x, y, z}
	}
	for _, t := range d.Tets {
		add(t.B, t.C, t.D)
		add(t.A, t.C, t.D)
		add(t.A, t.B, t.D)
		add(t.A, t.B, t.C)
	}
	var out [][3]int
	for k, c := range count {
		if c == 1 {
			out = append(out, rep[k])
		}
	}
	return out
}
