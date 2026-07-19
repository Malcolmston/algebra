package fem

import "math"

// AssembleStiffness1D assembles the global P1 stiffness matrix for the mesh.
func AssembleStiffness1D(m *Mesh1D) *SparseMatrix {
	n := m.NumNodes()
	A := NewSparseMatrix(n)
	for e := 0; e < m.NumElements(); e++ {
		i, j := m.ElementNodes(e)
		ke := ElementStiffnessP1Interval(m.ElementLength(e))
		idx := [2]int{i, j}
		for a := 0; a < 2; a++ {
			for b := 0; b < 2; b++ {
				A.AddEntry(idx[a], idx[b], ke[a][b])
			}
		}
	}
	return A
}

// AssembleMass1D assembles the global P1 mass matrix for the mesh.
func AssembleMass1D(m *Mesh1D) *SparseMatrix {
	n := m.NumNodes()
	M := NewSparseMatrix(n)
	for e := 0; e < m.NumElements(); e++ {
		i, j := m.ElementNodes(e)
		me := ElementMassP1Interval(m.ElementLength(e))
		idx := [2]int{i, j}
		for a := 0; a < 2; a++ {
			for b := 0; b < 2; b++ {
				M.AddEntry(idx[a], idx[b], me[a][b])
			}
		}
	}
	return M
}

// AssembleLoad1D assembles the global P1 load vector for source f.
func AssembleLoad1D(m *Mesh1D, f func(float64) float64) Vector {
	n := m.NumNodes()
	b := make(Vector, n)
	for e := 0; e < m.NumElements(); e++ {
		i, j := m.ElementNodes(e)
		fe := ElementLoadP1Interval(f, m.Nodes[i], m.Nodes[j], 3)
		b[i] += fe[0]
		b[j] += fe[1]
	}
	return b
}

// AssembleStiffness1DP2 assembles the global P2 stiffness matrix. The returned
// dimension equals NumNodes + NumElements; use Mesh1D.P2Nodes for coordinates.
func AssembleStiffness1DP2(m *Mesh1D) *SparseMatrix {
	dof := m.NumNodes() + m.NumElements()
	A := NewSparseMatrix(dof)
	conn := m.P2Connectivity()
	for e := 0; e < m.NumElements(); e++ {
		ke := ElementStiffnessP2Interval(m.ElementLength(e))
		idx := conn[e]
		for a := 0; a < 3; a++ {
			for b := 0; b < 3; b++ {
				A.AddEntry(idx[a], idx[b], ke[a][b])
			}
		}
	}
	return A
}

// AssembleMass1DP2 assembles the global P2 mass matrix.
func AssembleMass1DP2(m *Mesh1D) *SparseMatrix {
	dof := m.NumNodes() + m.NumElements()
	M := NewSparseMatrix(dof)
	conn := m.P2Connectivity()
	for e := 0; e < m.NumElements(); e++ {
		me := ElementMassP2Interval(m.ElementLength(e))
		idx := conn[e]
		for a := 0; a < 3; a++ {
			for b := 0; b < 3; b++ {
				M.AddEntry(idx[a], idx[b], me[a][b])
			}
		}
	}
	return M
}

// AssembleLoad1DP2 assembles the global P2 load vector for source f.
func AssembleLoad1DP2(m *Mesh1D, f func(float64) float64) Vector {
	dof := m.NumNodes() + m.NumElements()
	b := make(Vector, dof)
	conn := m.P2Connectivity()
	for e := 0; e < m.NumElements(); e++ {
		i, j := m.ElementNodes(e)
		fe := ElementLoadP2Interval(f, m.Nodes[i], m.Nodes[j], 4)
		idx := conn[e]
		for a := 0; a < 3; a++ {
			b[idx[a]] += fe[a]
		}
	}
	return b
}

// AssembleStiffness2D assembles the global P1 stiffness matrix for the mesh.
func AssembleStiffness2D(m *Mesh2D) *SparseMatrix {
	n := m.NumNodes()
	A := NewSparseMatrix(n)
	for t := 0; t < m.NumTriangles(); t++ {
		v1, v2, v3 := m.TriangleVertices(t)
		ke := ElementStiffnessP1Triangle(v1, v2, v3)
		tri := m.Triangles[t]
		for a := 0; a < 3; a++ {
			for b := 0; b < 3; b++ {
				A.AddEntry(tri[a], tri[b], ke[a][b])
			}
		}
	}
	return A
}

// AssembleMass2D assembles the global P1 mass matrix for the mesh.
func AssembleMass2D(m *Mesh2D) *SparseMatrix {
	n := m.NumNodes()
	M := NewSparseMatrix(n)
	for t := 0; t < m.NumTriangles(); t++ {
		v1, v2, v3 := m.TriangleVertices(t)
		me := ElementMassP1Triangle(v1, v2, v3)
		tri := m.Triangles[t]
		for a := 0; a < 3; a++ {
			for b := 0; b < 3; b++ {
				M.AddEntry(tri[a], tri[b], me[a][b])
			}
		}
	}
	return M
}

// AssembleLoad2D assembles the global P1 load vector for source f.
func AssembleLoad2D(m *Mesh2D, f func(x, y float64) float64) Vector {
	b := make(Vector, m.NumNodes())
	for t := 0; t < m.NumTriangles(); t++ {
		v1, v2, v3 := m.TriangleVertices(t)
		fe := ElementLoadP1Triangle(f, v1, v2, v3, 3)
		tri := m.Triangles[t]
		for a := 0; a < 3; a++ {
			b[tri[a]] += fe[a]
		}
	}
	return b
}

// AssembleReaction2D assembles kappa*M, the reaction (zeroth-order) contribution
// with constant coefficient kappa, into a new sparse matrix.
func AssembleReaction2D(m *Mesh2D, kappa float64) *SparseMatrix {
	M := AssembleMass2D(m)
	for i := 0; i < M.Dim(); i++ {
		M.ScaleRow(i, kappa)
	}
	return M
}

// AddSparse returns a+b for two sparse matrices of equal dimension.
func AddSparse(a, b *SparseMatrix) *SparseMatrix {
	if a.Dim() != b.Dim() {
		panic("fem: dimension mismatch in AddSparse")
	}
	out := NewSparseMatrix(a.Dim())
	for i := 0; i < a.Dim(); i++ {
		for j, v := range a.rows[i] {
			out.AddEntry(i, j, v)
		}
		for j, v := range b.rows[i] {
			out.AddEntry(i, j, v)
		}
	}
	return out
}

// ScaleSparse returns s*a for scalar s.
func ScaleSparse(a *SparseMatrix, s float64) *SparseMatrix {
	out := NewSparseMatrix(a.Dim())
	for i := 0; i < a.Dim(); i++ {
		for j, v := range a.rows[i] {
			out.AddEntry(i, j, s*v)
		}
	}
	return out
}

// ApplyDirichlet imposes the Dirichlet condition u[nodes[k]] = values[k] on the
// linear system (A, b) in place using the symmetric elimination technique that
// preserves symmetry of A. nodes and values must have equal length.
func ApplyDirichlet(A *SparseMatrix, b Vector, nodes []int, values []float64) {
	if len(nodes) != len(values) {
		panic("fem: nodes and values length mismatch in ApplyDirichlet")
	}
	fixed := make(map[int]float64, len(nodes))
	for k, node := range nodes {
		fixed[node] = values[k]
	}
	// Eliminate columns of fixed dofs from the right-hand side.
	for i := 0; i < A.Dim(); i++ {
		if _, isFixed := fixed[i]; isFixed {
			continue
		}
		for j, g := range fixed {
			if v, ok := A.rows[i][j]; ok {
				b[i] -= v * g
				delete(A.rows[i], j)
			}
		}
	}
	// Replace each fixed row with the identity row and set the value.
	for node, g := range fixed {
		A.ClearRow(node)
		A.rows[node][node] = 1
		b[node] = g
	}
}

// ApplyNeumann1D adds the natural (Neumann) boundary flux value at the given
// global node to the load vector b in place. For -u”=f the outward flux g is
// added directly.
func ApplyNeumann1D(b Vector, node int, flux float64) {
	b[node] += flux
}

// ApplyRobin1D imposes a Robin condition alpha*u = g (contributing to the
// diagonal) at the given global node by adding alpha to A[node,node] and g to
// b[node] in place.
func ApplyRobin1D(A *SparseMatrix, b Vector, node int, alpha, g float64) {
	A.AddEntry(node, node, alpha)
	b[node] += g
}

// ApplyNeumann2D adds the Neumann flux contribution of g over the given
// boundary edges to the load vector b (P1). The integral of g*N along each edge
// is distributed to its two endpoints.
func ApplyNeumann2D(m *Mesh2D, b Vector, edges []Edge, g func(x, y float64) float64) {
	q := GaussLegendre1D(3)
	for _, e := range edges {
		pa, pb := m.Nodes[e.A], m.Nodes[e.B]
		dx := pb[0] - pa[0]
		dy := pb[1] - pa[1]
		L := math.Sqrt(dx*dx + dy*dy)
		nodes, weights := q.MapToInterval(0, 1)
		for k := range nodes {
			s := nodes[k]
			x := pa[0] + s*dx
			y := pa[1] + s*dy
			gv := g(x, y) * L
			b[e.A] += weights[k] * gv * (1 - s)
			b[e.B] += weights[k] * gv * s
		}
	}
}

// ApplyRobin2D adds Robin contributions alpha*u = g over the given boundary
// edges, augmenting both the stiffness matrix A (edge mass scaled by alpha) and
// the load vector b (edge integral of g), in place (P1).
func ApplyRobin2D(m *Mesh2D, A *SparseMatrix, b Vector, edges []Edge, alpha float64, g func(x, y float64) float64) {
	for _, e := range edges {
		pa, pb := m.Nodes[e.A], m.Nodes[e.B]
		dx := pb[0] - pa[0]
		dy := pb[1] - pa[1]
		L := math.Sqrt(dx*dx + dy*dy)
		// Edge mass matrix alpha*(L/6)[[2,1],[1,2]].
		c := alpha * L / 6
		A.AddEntry(e.A, e.A, 2*c)
		A.AddEntry(e.A, e.B, c)
		A.AddEntry(e.B, e.A, c)
		A.AddEntry(e.B, e.B, 2*c)
		// Load integral of g*N.
		fe := edgeLoadP1(pa, pb, g)
		b[e.A] += fe[0]
		b[e.B] += fe[1]
	}
}

func edgeLoadP1(pa, pb [2]float64, g func(x, y float64) float64) [2]float64 {
	q := GaussLegendre1D(3)
	dx := pb[0] - pa[0]
	dy := pb[1] - pa[1]
	L := math.Sqrt(dx*dx + dy*dy)
	nodes, weights := q.MapToInterval(0, 1)
	var out [2]float64
	for k := range nodes {
		s := nodes[k]
		x := pa[0] + s*dx
		y := pa[1] + s*dy
		gv := g(x, y) * L
		out[0] += weights[k] * gv * (1 - s)
		out[1] += weights[k] * gv * s
	}
	return out
}
