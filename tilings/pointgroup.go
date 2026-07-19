package tilings

import (
	"math"
	"sort"
)

// PointGroup is a finite group of planar isometries that all fix a common
// point (taken to be the origin). Its elements are rotations about the origin
// and reflections through the origin.
type PointGroup struct {
	// Name is a short label such as "C4" or "D6".
	Name string
	// Elements are the isometries of the group.
	Elements []Isometry
}

// Order returns the number of elements in the point group.
func (g PointGroup) Order() int { return len(g.Elements) }

// CyclicGroup returns the cyclic point group C_n of rotations about the origin
// by multiples of 2*pi/n. It panics for n < 1.
func CyclicGroup(n int) PointGroup {
	if n < 1 {
		panic("tilings: CyclicGroup requires n >= 1")
	}
	els := make([]Isometry, n)
	for k := 0; k < n; k++ {
		els[k] = Rotation(2 * math.Pi * float64(k) / float64(n))
	}
	return PointGroup{Name: cyclicName(n), Elements: els}
}

// DihedralGroup returns the dihedral point group D_n of order 2n: the n
// rotations of C_n together with n reflections whose axes pass through the
// origin, the first reflection axis being the x-axis. It panics for n < 1.
func DihedralGroup(n int) PointGroup {
	if n < 1 {
		panic("tilings: DihedralGroup requires n >= 1")
	}
	els := make([]Isometry, 0, 2*n)
	for k := 0; k < n; k++ {
		els = append(els, Rotation(2*math.Pi*float64(k)/float64(n)))
	}
	for k := 0; k < n; k++ {
		els = append(els, Reflection(math.Pi*float64(k)/float64(n)))
	}
	return PointGroup{Name: dihedralName(n), Elements: els}
}

func cyclicName(n int) string   { return "C" + itoa(n) }
func dihedralName(n int) string { return "D" + itoa(n) }

// itoa is a tiny non-allocating-ish integer formatter for small positive n.
func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}

// RotationOrders returns the sorted set of rotation orders present in the point
// group, where a rotation by 2*pi/k contributes k. The identity is excluded.
func (g PointGroup) RotationOrders() []int {
	seen := map[int]bool{}
	for _, e := range g.Elements {
		if !e.IsDirect() {
			continue
		}
		a := e.RotationAngle()
		if approxEqualScalar(a, 0, 1e-9) {
			continue
		}
		k := int(math.Round(2 * math.Pi / math.Abs(a)))
		if k >= 2 {
			seen[k] = true
		}
	}
	out := make([]int, 0, len(seen))
	for k := range seen {
		out = append(out, k)
	}
	sort.Ints(out)
	return out
}

// MaxRotationOrder returns the highest rotation order in the point group, or 1
// if the group contains no nontrivial rotation.
func (g PointGroup) MaxRotationOrder() int {
	max := 1
	for _, k := range g.RotationOrders() {
		if k > max {
			max = k
		}
	}
	return max
}

// NumReflections returns the number of reflections (opposite isometries) in the
// point group.
func (g PointGroup) NumReflections() int {
	n := 0
	for _, e := range g.Elements {
		if !e.IsDirect() {
			n++
		}
	}
	return n
}

// HasReflection reports whether the point group contains any reflection.
func (g PointGroup) HasReflection() bool { return g.NumReflections() > 0 }

// Contains reports whether the point group contains an isometry equal to a
// within eps.
func (g PointGroup) Contains(a Isometry, eps float64) bool {
	for _, e := range g.Elements {
		if e.ApproxEqual(a, eps) {
			return true
		}
	}
	return false
}

// IsClosed reports whether the point group is closed under composition to
// within eps: the product of any two elements is again an element.
func (g PointGroup) IsClosed(eps float64) bool {
	for _, x := range g.Elements {
		for _, y := range g.Elements {
			if !g.Contains(x.Compose(y), eps) {
				return false
			}
		}
	}
	return true
}

// Orbit returns the distinct images of the point p under every element of the
// point group, deduplicated to within eps.
func (g PointGroup) Orbit(p Point, eps float64) []Point {
	var out []Point
	for _, e := range g.Elements {
		q := e.Apply(p)
		dup := false
		for _, r := range out {
			if r.ApproxEqual(q, eps) {
				dup = true
				break
			}
		}
		if !dup {
			out = append(out, q)
		}
	}
	return out
}

// StabilizerOrder returns the number of group elements fixing the point p to
// within eps.
func (g PointGroup) StabilizerOrder(p Point, eps float64) int {
	n := 0
	for _, e := range g.Elements {
		if e.Apply(p).ApproxEqual(p, eps) {
			n++
		}
	}
	return n
}

// CrystallographicPointGroups returns the ten two-dimensional crystallographic
// point groups (C1, C2, C3, C4, C6 and D1, D2, D3, D4, D6) in that order.
func CrystallographicPointGroups() []PointGroup {
	orders := []int{1, 2, 3, 4, 6}
	out := make([]PointGroup, 0, 10)
	for _, n := range orders {
		out = append(out, CyclicGroup(n))
	}
	for _, n := range orders {
		out = append(out, DihedralGroup(n))
	}
	return out
}

// IsCrystallographicRestriction reports whether n is an allowed rotation order
// for a two-dimensional crystallographic (lattice-preserving) symmetry, i.e.
// n is 1, 2, 3, 4 or 6.
func IsCrystallographicRestriction(n int) bool {
	switch n {
	case 1, 2, 3, 4, 6:
		return true
	default:
		return false
	}
}
