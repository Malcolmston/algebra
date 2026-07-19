package simplicial

import (
	"sort"
	"strings"
)

// Complex is an abstract simplicial complex: a finite, downward-closed family
// of [Simplex] values. Downward-closed means that whenever a simplex belongs to
// the complex, so does every one of its faces; the mutating constructors below
// preserve this invariant automatically.
type Complex struct {
	simplices map[string]Simplex
}

// NewComplex returns a new, empty simplicial complex.
func NewComplex() *Complex {
	return &Complex{simplices: make(map[string]Simplex)}
}

// FromSimplices builds the smallest simplicial complex containing all of the
// given simplices, i.e. their union together with every face.
func FromSimplices(ss ...Simplex) *Complex {
	c := NewComplex()
	for _, s := range ss {
		c.AddSimplex(s)
	}
	return c
}

// AddSimplex inserts the simplex and every one of its faces into the complex,
// preserving downward closure. Adding a simplex already present is a no-op.
func (c *Complex) AddSimplex(s Simplex) {
	if s.IsEmpty() {
		return
	}
	for _, f := range s.Closure() {
		c.simplices[f.Key()] = f
	}
}

// AddSimplices inserts several simplices and their faces.
func (c *Complex) AddSimplices(ss ...Simplex) {
	for _, s := range ss {
		c.AddSimplex(s)
	}
}

// AddVertex inserts the 0-simplex {v}.
func (c *Complex) AddVertex(v int) { c.AddSimplex(Vertex(v)) }

// AddEdge inserts the 1-simplex {a,b} and its two endpoints.
func (c *Complex) AddEdge(a, b int) { c.AddSimplex(Edge(a, b)) }

// AddTriangle inserts the 2-simplex {a,b,c} and all of its faces.
func (c *Complex) AddTriangle(a, b, c2 int) { c.AddSimplex(Triangle(a, b, c2)) }

// Has reports whether the given simplex belongs to the complex.
func (c *Complex) Has(s Simplex) bool {
	_, ok := c.simplices[s.Key()]
	return ok
}

// HasVertex reports whether the 0-simplex {v} belongs to the complex.
func (c *Complex) HasVertex(v int) bool { return c.Has(Vertex(v)) }

// NumSimplices returns the total number of simplices of every dimension.
func (c *Complex) NumSimplices() int { return len(c.simplices) }

// Simplices returns all simplices of the complex sorted by [CompareSimplices].
func (c *Complex) Simplices() []Simplex {
	out := make([]Simplex, 0, len(c.simplices))
	for _, s := range c.simplices {
		out = append(out, s)
	}
	SortSimplices(out)
	return out
}

// Dimension returns the largest dimension of any simplex in the complex, or −1
// for the empty complex.
func (c *Complex) Dimension() int {
	d := -1
	for _, s := range c.simplices {
		if s.Dim() > d {
			d = s.Dim()
		}
	}
	return d
}

// SimplicesOfDim returns all k-simplices of the complex sorted by
// [CompareSimplices].
func (c *Complex) SimplicesOfDim(k int) []Simplex {
	var out []Simplex
	for _, s := range c.simplices {
		if s.Dim() == k {
			out = append(out, s)
		}
	}
	SortSimplices(out)
	return out
}

// NumSimplicesOfDim returns the number of k-simplices in the complex.
func (c *Complex) NumSimplicesOfDim(k int) int {
	n := 0
	for _, s := range c.simplices {
		if s.Dim() == k {
			n++
		}
	}
	return n
}

// FVector returns the f-vector of the complex: entry k is the number of
// k-simplices, for k from 0 up to the dimension of the complex. The empty
// complex yields an empty slice.
func (c *Complex) FVector() []int {
	d := c.Dimension()
	if d < 0 {
		return nil
	}
	f := make([]int, d+1)
	for _, s := range c.simplices {
		f[s.Dim()]++
	}
	return f
}

// Vertices returns the sorted list of vertices appearing in the complex.
func (c *Complex) Vertices() []int {
	var vs []int
	for _, s := range c.simplices {
		if s.Dim() == 0 {
			vs = append(vs, s.verts[0])
		}
	}
	sort.Ints(vs)
	return vs
}

// NumVertices returns the number of 0-simplices (vertices) of the complex.
func (c *Complex) NumVertices() int { return c.NumSimplicesOfDim(0) }

// EulerCharacteristic returns the alternating sum Σ_k (−1)^k f_k of the
// f-vector, the Euler characteristic of the complex.
func (c *Complex) EulerCharacteristic() int {
	chi := 0
	sign := 1
	for _, fk := range c.FVector() {
		chi += sign * fk
		sign = -sign
	}
	return chi
}

// Skeleton returns the k-skeleton of the complex: a new complex containing
// exactly the simplices of dimension at most k. A negative k yields the empty
// complex.
func (c *Complex) Skeleton(k int) *Complex {
	out := NewComplex()
	for key, s := range c.simplices {
		if s.Dim() <= k {
			out.simplices[key] = s
		}
	}
	return out
}

// Copy returns an independent copy of the complex.
func (c *Complex) Copy() *Complex {
	out := NewComplex()
	for k, v := range c.simplices {
		out.simplices[k] = v.Copy()
	}
	return out
}

// Equal reports whether c and d contain exactly the same simplices.
func (c *Complex) Equal(d *Complex) bool {
	if len(c.simplices) != len(d.simplices) {
		return false
	}
	for k := range c.simplices {
		if _, ok := d.simplices[k]; !ok {
			return false
		}
	}
	return true
}

// IsValid reports whether the complex is genuinely downward-closed. The
// mutating methods maintain this, so it should always be true; it exists for
// testing and for complexes assembled by hand.
func (c *Complex) IsValid() bool {
	for _, s := range c.simplices {
		for _, f := range s.Faces() {
			if _, ok := c.simplices[f.Key()]; !ok {
				return false
			}
		}
	}
	return true
}

// Cofaces returns every simplex of the complex that has s as a face, s itself
// included when present.
func (c *Complex) Cofaces(s Simplex) []Simplex {
	var out []Simplex
	for _, t := range c.simplices {
		if t.ContainsSimplex(s) {
			out = append(out, t)
		}
	}
	SortSimplices(out)
	return out
}

// Star returns the open star of s: all simplices of the complex having s as a
// face. It is an alias for [Complex.Cofaces].
func (c *Complex) Star(s Simplex) []Simplex { return c.Cofaces(s) }

// ClosedStar returns the closure of the star of s as a subcomplex: the smallest
// subcomplex containing every coface of s.
func (c *Complex) ClosedStar(s Simplex) *Complex {
	out := NewComplex()
	for _, t := range c.Cofaces(s) {
		out.AddSimplex(t)
	}
	return out
}

// Link returns the link of s: the subcomplex of faces τ that are disjoint from
// s yet whose join with s lies in the complex. The link of a vertex in a
// surface triangulation is a cycle.
func (c *Complex) Link(s Simplex) *Complex {
	out := NewComplex()
	for _, t := range c.simplices {
		if !t.ContainsSimplex(s) {
			continue
		}
		// t contains s, so τ = t∖s is disjoint from s and τ∪s = t ∈ K.
		rest := removeVertices(t, s)
		if rest.IsEmpty() {
			continue
		}
		out.AddSimplex(rest)
	}
	return out
}

// removeVertices returns the simplex t with the vertices of s deleted.
func removeVertices(t, s Simplex) Simplex {
	var vs []int
	for _, v := range t.verts {
		if !s.Contains(v) {
			vs = append(vs, v)
		}
	}
	return Simplex{verts: vs}
}

// RemoveSimplex deletes s and every coface of s from the complex, keeping it
// downward-closed.
func (c *Complex) RemoveSimplex(s Simplex) {
	for _, t := range c.Cofaces(s) {
		delete(c.simplices, t.Key())
	}
}

// MaximalSimplices returns the maximal simplices of the complex — those that
// are not a proper face of any other simplex. A complex is determined by its
// maximal simplices.
func (c *Complex) MaximalSimplices() []Simplex {
	var out []Simplex
	for _, s := range c.simplices {
		maximal := true
		for _, t := range c.simplices {
			if !t.Equal(s) && t.ContainsSimplex(s) {
				maximal = false
				break
			}
		}
		if maximal {
			out = append(out, s)
		}
	}
	SortSimplices(out)
	return out
}

// Union returns the simplicial complex whose simplices are those of c or d.
func (c *Complex) Union(d *Complex) *Complex {
	out := c.Copy()
	for k, v := range d.simplices {
		out.simplices[k] = v
	}
	return out
}

// String returns a compact multi-line description listing the simplices by
// dimension.
func (c *Complex) String() string {
	var b strings.Builder
	b.WriteString("Complex{")
	ss := c.Simplices()
	for i, s := range ss {
		if i > 0 {
			b.WriteByte(' ')
		}
		b.WriteString(s.Key())
	}
	b.WriteByte('}')
	return b.String()
}
