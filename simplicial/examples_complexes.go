package simplicial

// PointComplex returns the complex consisting of a single vertex {0}: a
// contractible space with the homology of a point.
func PointComplex() *Complex { return FromSimplices(Vertex(0)) }

// StandardSimplex returns the full n-simplex on vertices 0..n together with all
// of its faces. It is contractible for every n ≥ 0.
func StandardSimplex(n int) *Complex {
	vs := make([]int, n+1)
	for i := range vs {
		vs[i] = i
	}
	return FromSimplices(NewSimplex(vs...))
}

// BoundarySphere returns the boundary of the standard (n+1)-simplex: the
// n-skeleton of an (n+1)-simplex, a triangulation of the n-sphere S^n on n+2
// vertices. For n ≥ 1 it has the homology of S^n.
func BoundarySphere(n int) *Complex {
	return StandardSimplex(n + 1).Skeleton(n)
}

// SphereComplex returns a triangulation of the n-sphere S^n as the boundary of
// the standard (n+1)-simplex. SphereComplex(0) is two points; SphereComplex(1)
// is a triangle boundary (a circle); SphereComplex(2) is a hollow tetrahedron.
func SphereComplex(n int) *Complex { return BoundarySphere(n) }

// CycleComplex returns the cycle graph C_n: n vertices 0..n−1 joined in a ring
// by n edges. For n ≥ 3 it is a triangulation of the circle S^1.
func CycleComplex(n int) *Complex {
	c := NewComplex()
	if n <= 0 {
		return c
	}
	if n == 1 {
		c.AddVertex(0)
		return c
	}
	for i := 0; i < n; i++ {
		c.AddEdge(i, (i+1)%n)
	}
	return c
}

// PathComplex returns the path graph on n vertices 0..n−1 with n−1 edges. It is
// contractible for n ≥ 1.
func PathComplex(n int) *Complex {
	c := NewComplex()
	if n <= 0 {
		return c
	}
	c.AddVertex(0)
	for i := 0; i+1 < n; i++ {
		c.AddEdge(i, i+1)
	}
	return c
}

// CompleteGraphComplex returns the 1-skeleton of the full simplex on n vertices:
// the complete graph K_n with all n(n−1)/2 edges but no higher faces.
func CompleteGraphComplex(n int) *Complex {
	c := NewComplex()
	for i := 0; i < n; i++ {
		c.AddVertex(i)
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			c.AddEdge(i, j)
		}
	}
	return c
}

// DiscreteComplex returns n isolated vertices 0..n−1 with no edges; its degree-0
// homology has rank n and all higher homology vanishes.
func DiscreteComplex(n int) *Complex {
	c := NewComplex()
	for i := 0; i < n; i++ {
		c.AddVertex(i)
	}
	return c
}

// TorusComplex returns Möbius's minimal 7-vertex triangulation of the torus
// T^2, generated on vertices 0..6 (indices taken mod 7) by the triangles
// {i,i+1,i+3} and {i,i+2,i+3}. It has 7 vertices, 21 edges and 14 triangles
// (Euler characteristic 0) and its homology is H_0 = Z, H_1 = Z^2, H_2 = Z.
func TorusComplex() *Complex {
	c := NewComplex()
	for i := 0; i < 7; i++ {
		c.AddTriangle(i%7, (i+1)%7, (i+3)%7)
		c.AddTriangle(i%7, (i+2)%7, (i+3)%7)
	}
	return c
}

// ProjectivePlaneComplex returns the minimal 6-vertex triangulation of the real
// projective plane RP^2 (the hemi-icosahedron) on vertices 0..5, with 15 edges
// and 10 triangles (Euler characteristic 1). Over Z its homology is H_0 = Z,
// H_1 = Z/2, H_2 = 0; over GF(2) the Betti numbers are (1,1,1). It is the
// smallest complex in the package with non-trivial torsion.
func ProjectivePlaneComplex() *Complex {
	tris := [][3]int{
		{0, 1, 2}, {0, 1, 3}, {0, 2, 4}, {0, 3, 5}, {0, 4, 5},
		{1, 2, 5}, {1, 3, 4}, {1, 4, 5}, {2, 3, 4}, {2, 3, 5},
	}
	c := NewComplex()
	for _, t := range tris {
		c.AddTriangle(t[0], t[1], t[2])
	}
	return c
}

// DisjointUnion returns the disjoint union of c and d. The vertices of d are
// relabelled by adding the given offset so that the two vertex sets do not
// collide; choose an offset larger than any vertex of c.
func DisjointUnion(c, d *Complex, offset int) *Complex {
	out := c.Copy()
	for _, s := range d.Simplices() {
		vs := make([]int, s.Len())
		for i, v := range s.Vertices() {
			vs[i] = v + offset
		}
		out.AddSimplex(NewSimplex(vs...))
	}
	return out
}

// Cone returns the cone on c with a new apex vertex: every simplex σ of c is
// joined to apex, together with the apex itself. The cone is always
// contractible. The apex must differ from every vertex of c.
func Cone(c *Complex, apex int) *Complex {
	out := c.Copy()
	out.AddVertex(apex)
	for _, s := range c.Simplices() {
		out.AddSimplex(s.Join(Vertex(apex)))
	}
	return out
}

// Suspension returns the suspension of c: two cones on c glued along c, using
// the two new apex vertices north and south. The suspension of S^n is S^{n+1}.
// The apex vertices must differ from each other and from every vertex of c.
func Suspension(c *Complex, north, south int) *Complex {
	out := c.Copy()
	out.AddVertex(north)
	out.AddVertex(south)
	for _, s := range c.Simplices() {
		out.AddSimplex(s.Join(Vertex(north)))
		out.AddSimplex(s.Join(Vertex(south)))
	}
	return out
}
