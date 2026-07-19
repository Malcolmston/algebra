package simplicial

import (
	"sort"
	"strconv"
	"strings"
)

// Simplex is a finite set of distinct integer vertices, stored in ascending
// order. A simplex with k+1 vertices is a k-simplex and has dimension k; the
// empty simplex has dimension −1. Simplices are value types: methods never
// mutate the receiver, they return fresh simplices.
type Simplex struct {
	verts []int
}

// NewSimplex returns the simplex on the given vertices. Repeated vertices are
// collapsed and the result is sorted ascending, so NewSimplex(2,0,2,1) and
// NewSimplex(0,1,2) denote the same 2-simplex.
func NewSimplex(verts ...int) Simplex {
	if len(verts) == 0 {
		return Simplex{}
	}
	cp := append([]int(nil), verts...)
	sort.Ints(cp)
	out := cp[:0]
	var last int
	for i, v := range cp {
		if i == 0 || v != last {
			out = append(out, v)
			last = v
		}
	}
	return Simplex{verts: out}
}

// Vertex returns the 0-simplex on the single vertex v.
func Vertex(v int) Simplex { return Simplex{verts: []int{v}} }

// Edge returns the 1-simplex on vertices a and b. If a equals b the result is
// the 0-simplex {a}.
func Edge(a, b int) Simplex { return NewSimplex(a, b) }

// Triangle returns the 2-simplex on vertices a, b and c.
func Triangle(a, b, c int) Simplex { return NewSimplex(a, b, c) }

// Tetrahedron returns the 3-simplex on vertices a, b, c and d.
func Tetrahedron(a, b, c, d int) Simplex { return NewSimplex(a, b, c, d) }

// EmptySimplex returns the (−1)-dimensional simplex with no vertices.
func EmptySimplex() Simplex { return Simplex{} }

// Dim returns the dimension of the simplex, one less than its number of
// vertices. The empty simplex has dimension −1.
func (s Simplex) Dim() int { return len(s.verts) - 1 }

// Len returns the number of vertices of the simplex.
func (s Simplex) Len() int { return len(s.verts) }

// IsEmpty reports whether the simplex has no vertices.
func (s Simplex) IsEmpty() bool { return len(s.verts) == 0 }

// Vertices returns a copy of the simplex's vertex list in ascending order.
func (s Simplex) Vertices() []int { return append([]int(nil), s.verts...) }

// Copy returns an independent copy of the simplex.
func (s Simplex) Copy() Simplex { return Simplex{verts: append([]int(nil), s.verts...)} }

// Vertex returns the i-th vertex (0-indexed) of the simplex in ascending
// order. It panics if i is out of range, matching slice indexing.
func (s Simplex) Vertex(i int) int { return s.verts[i] }

// Min returns the smallest vertex, or 0 for the empty simplex.
func (s Simplex) Min() int {
	if len(s.verts) == 0 {
		return 0
	}
	return s.verts[0]
}

// Max returns the largest vertex, or 0 for the empty simplex.
func (s Simplex) Max() int {
	if len(s.verts) == 0 {
		return 0
	}
	return s.verts[len(s.verts)-1]
}

// Contains reports whether v is a vertex of the simplex.
func (s Simplex) Contains(v int) bool {
	i := sort.SearchInts(s.verts, v)
	return i < len(s.verts) && s.verts[i] == v
}

// Key returns a canonical string key identifying the simplex, suitable as a map
// key. Distinct simplices have distinct keys.
func (s Simplex) Key() string {
	if len(s.verts) == 0 {
		return "∅"
	}
	var b strings.Builder
	for i, v := range s.verts {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(strconv.Itoa(v))
	}
	return b.String()
}

// String returns a human-readable representation such as "{0,1,2}".
func (s Simplex) String() string {
	if len(s.verts) == 0 {
		return "{}"
	}
	return "{" + s.Key() + "}"
}

// Equal reports whether s and t have exactly the same vertex set.
func (s Simplex) Equal(t Simplex) bool {
	if len(s.verts) != len(t.verts) {
		return false
	}
	for i := range s.verts {
		if s.verts[i] != t.verts[i] {
			return false
		}
	}
	return true
}

// ContainsSimplex reports whether t is a (not necessarily proper) face of s,
// that is, whether every vertex of t is a vertex of s.
func (s Simplex) ContainsSimplex(t Simplex) bool {
	if len(t.verts) > len(s.verts) {
		return false
	}
	i := 0
	for _, v := range t.verts {
		for i < len(s.verts) && s.verts[i] < v {
			i++
		}
		if i >= len(s.verts) || s.verts[i] != v {
			return false
		}
		i++
	}
	return true
}

// IsFaceOf reports whether s is a face of t (equivalent to
// t.ContainsSimplex(s)).
func (s Simplex) IsFaceOf(t Simplex) bool { return t.ContainsSimplex(s) }

// IsProperFaceOf reports whether s is a face of t of strictly lower dimension.
func (s Simplex) IsProperFaceOf(t Simplex) bool {
	return len(s.verts) < len(t.verts) && t.ContainsSimplex(s)
}

// Face returns the codimension-1 face obtained by deleting the i-th vertex
// (0-indexed, ascending order) of the simplex.
func (s Simplex) Face(i int) Simplex {
	nv := make([]int, 0, len(s.verts)-1)
	nv = append(nv, s.verts[:i]...)
	nv = append(nv, s.verts[i+1:]...)
	return Simplex{verts: nv}
}

// Faces returns the codimension-1 faces of the simplex: for a k-simplex the
// k+1 sub-simplices obtained by omitting one vertex. A 0-simplex has no faces
// (an empty slice is returned).
func (s Simplex) Faces() []Simplex {
	if len(s.verts) <= 1 {
		return nil
	}
	out := make([]Simplex, len(s.verts))
	for i := range s.verts {
		out[i] = s.Face(i)
	}
	return out
}

// BoundaryTerm is one oriented summand of a simplex's boundary: a face together
// with the sign (+1 or −1) with which it appears.
type BoundaryTerm struct {
	// Sign is +1 or −1.
	Sign int
	// Face is the codimension-1 face.
	Face Simplex
}

// Boundary returns the oriented algebraic boundary ∂s of the simplex, the
// alternating sum Σ_i (−1)^i [v_0,…,v̂_i,…,v_k] of its codimension-1 faces. A
// 0-simplex has empty boundary.
func (s Simplex) Boundary() []BoundaryTerm {
	if len(s.verts) <= 1 {
		return nil
	}
	out := make([]BoundaryTerm, len(s.verts))
	for i := range s.verts {
		sign := 1
		if i%2 == 1 {
			sign = -1
		}
		out[i] = BoundaryTerm{Sign: sign, Face: s.Face(i)}
	}
	return out
}

// Subsets returns every non-empty subset of the vertices of s as a simplex,
// including s itself. A k-simplex has 2^{k+1}−1 non-empty subsets.
func (s Simplex) Subsets() []Simplex {
	n := len(s.verts)
	if n == 0 {
		return nil
	}
	out := make([]Simplex, 0, (1<<uint(n))-1)
	for mask := 1; mask < (1 << uint(n)); mask++ {
		vs := make([]int, 0, n)
		for i := 0; i < n; i++ {
			if mask&(1<<uint(i)) != 0 {
				vs = append(vs, s.verts[i])
			}
		}
		out = append(out, Simplex{verts: vs})
	}
	return out
}

// Closure returns all faces of the simplex, itself included: every non-empty
// subset of its vertices. It is an alias for [Simplex.Subsets].
func (s Simplex) Closure() []Simplex { return s.Subsets() }

// ProperFaces returns every proper face of the simplex (all non-empty subsets
// except s itself).
func (s Simplex) ProperFaces() []Simplex {
	all := s.Subsets()
	out := all[:0:0]
	for _, f := range all {
		if !f.Equal(s) {
			out = append(out, f)
		}
	}
	return out
}

// Join returns the simplex whose vertex set is the union of those of s and t.
func (s Simplex) Join(t Simplex) Simplex {
	vs := append(append([]int(nil), s.verts...), t.verts...)
	return NewSimplex(vs...)
}

// Intersect returns the simplex on the common vertices of s and t. The result
// may be the empty simplex.
func (s Simplex) Intersect(t Simplex) Simplex {
	var out []int
	i, j := 0, 0
	for i < len(s.verts) && j < len(t.verts) {
		switch {
		case s.verts[i] < t.verts[j]:
			i++
		case s.verts[i] > t.verts[j]:
			j++
		default:
			out = append(out, s.verts[i])
			i++
			j++
		}
	}
	return Simplex{verts: out}
}

// CompareSimplices provides a total order on simplices: first by dimension,
// then lexicographically by vertex list. It returns −1, 0 or +1.
func CompareSimplices(a, b Simplex) int {
	if da, db := a.Dim(), b.Dim(); da != db {
		if da < db {
			return -1
		}
		return 1
	}
	for i := range a.verts {
		if a.verts[i] != b.verts[i] {
			if a.verts[i] < b.verts[i] {
				return -1
			}
			return 1
		}
	}
	return 0
}

// SortSimplices sorts a slice of simplices in place by [CompareSimplices].
func SortSimplices(ss []Simplex) {
	sort.Slice(ss, func(i, j int) bool { return CompareSimplices(ss[i], ss[j]) < 0 })
}
