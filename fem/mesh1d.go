package fem

import (
	"errors"
	"math"
	"sort"
)

// Mesh1D is a one-dimensional mesh of an interval, described by sorted node
// coordinates. Elements connect consecutive nodes.
type Mesh1D struct {
	Nodes []float64
}

// NewMesh1D builds a mesh from the given node coordinates. The coordinates are
// copied and sorted; duplicate coordinates return an error.
func NewMesh1D(nodes []float64) (*Mesh1D, error) {
	if len(nodes) < 2 {
		return nil, errors.New("fem: a 1D mesh needs at least two nodes")
	}
	c := make([]float64, len(nodes))
	copy(c, nodes)
	sort.Float64s(c)
	for i := 1; i < len(c); i++ {
		if c[i] == c[i-1] {
			return nil, errors.New("fem: duplicate node in 1D mesh")
		}
	}
	return &Mesh1D{Nodes: c}, nil
}

// NewUniformMesh1D returns a uniform mesh of the interval [a,b] with n elements
// (n+1 equally spaced nodes). It panics if n < 1 or b <= a.
func NewUniformMesh1D(a, b float64, n int) *Mesh1D {
	if n < 1 {
		panic("fem: NewUniformMesh1D requires n >= 1")
	}
	if b <= a {
		panic("fem: NewUniformMesh1D requires b > a")
	}
	nodes := make([]float64, n+1)
	h := (b - a) / float64(n)
	for i := 0; i <= n; i++ {
		nodes[i] = a + float64(i)*h
	}
	nodes[n] = b
	return &Mesh1D{Nodes: nodes}
}

// NumNodes returns the number of nodes.
func (m *Mesh1D) NumNodes() int { return len(m.Nodes) }

// NumElements returns the number of elements.
func (m *Mesh1D) NumElements() int { return len(m.Nodes) - 1 }

// ElementNodes returns the global node indices (left, right) of element e.
func (m *Mesh1D) ElementNodes(e int) (int, int) { return e, e + 1 }

// ElementLength returns the length of element e.
func (m *Mesh1D) ElementLength(e int) float64 { return m.Nodes[e+1] - m.Nodes[e] }

// ElementMidpoint returns the midpoint coordinate of element e.
func (m *Mesh1D) ElementMidpoint(e int) float64 { return 0.5 * (m.Nodes[e] + m.Nodes[e+1]) }

// Domain returns the endpoints of the meshed interval.
func (m *Mesh1D) Domain() (a, b float64) { return m.Nodes[0], m.Nodes[len(m.Nodes)-1] }

// MinElementLength returns the length of the shortest element.
func (m *Mesh1D) MinElementLength() float64 {
	min := m.ElementLength(0)
	for e := 1; e < m.NumElements(); e++ {
		if h := m.ElementLength(e); h < min {
			min = h
		}
	}
	return min
}

// MaxElementLength returns the length of the longest element (the mesh size h).
func (m *Mesh1D) MaxElementLength() float64 {
	max := m.ElementLength(0)
	for e := 1; e < m.NumElements(); e++ {
		if h := m.ElementLength(e); h > max {
			max = h
		}
	}
	return max
}

// BoundaryNodes returns the indices of the two boundary nodes (the endpoints).
func (m *Mesh1D) BoundaryNodes() []int {
	return []int{0, len(m.Nodes) - 1}
}

// InteriorNodes returns the indices of the interior nodes.
func (m *Mesh1D) InteriorNodes() []int {
	out := make([]int, 0, len(m.Nodes)-2)
	for i := 1; i < len(m.Nodes)-1; i++ {
		out = append(out, i)
	}
	return out
}

// Refine returns a new mesh obtained by bisecting every element of m.
func (m *Mesh1D) Refine() *Mesh1D {
	nodes := make([]float64, 0, 2*len(m.Nodes)-1)
	for e := 0; e < m.NumElements(); e++ {
		nodes = append(nodes, m.Nodes[e], m.ElementMidpoint(e))
	}
	nodes = append(nodes, m.Nodes[len(m.Nodes)-1])
	return &Mesh1D{Nodes: nodes}
}

// RefineN applies Refine k times.
func (m *Mesh1D) RefineN(k int) *Mesh1D {
	out := m
	for i := 0; i < k; i++ {
		out = out.Refine()
	}
	return out
}

// P2Nodes returns the coordinates of all P2 degrees of freedom: the mesh nodes
// followed by the element midpoints. The returned slice has length
// NumNodes + NumElements and P2Connectivity indexes into it.
func (m *Mesh1D) P2Nodes() []float64 {
	nn := m.NumNodes()
	out := make([]float64, nn+m.NumElements())
	copy(out, m.Nodes)
	for e := 0; e < m.NumElements(); e++ {
		out[nn+e] = m.ElementMidpoint(e)
	}
	return out
}

// P2Connectivity returns, for each element, the three global P2 degree-of-
// freedom indices (left node, right node, midpoint) indexing into P2Nodes.
func (m *Mesh1D) P2Connectivity() [][3]int {
	nn := m.NumNodes()
	conn := make([][3]int, m.NumElements())
	for e := 0; e < m.NumElements(); e++ {
		conn[e] = [3]int{e, e + 1, nn + e}
	}
	return conn
}

// GradedMesh1D returns a mesh of [a,b] with n elements whose node spacing is
// graded towards a by the power exponent (exponent 1 is uniform, exponents > 1
// cluster nodes near a). It panics if n < 1, b <= a or exponent <= 0.
func GradedMesh1D(a, b float64, n int, exponent float64) *Mesh1D {
	if n < 1 {
		panic("fem: GradedMesh1D requires n >= 1")
	}
	if b <= a {
		panic("fem: GradedMesh1D requires b > a")
	}
	if exponent <= 0 {
		panic("fem: GradedMesh1D requires exponent > 0")
	}
	nodes := make([]float64, n+1)
	for i := 0; i <= n; i++ {
		t := float64(i) / float64(n)
		nodes[i] = a + (b-a)*math.Pow(t, exponent)
	}
	nodes[n] = b
	return &Mesh1D{Nodes: nodes}
}
