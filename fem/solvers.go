package fem

// DirichletData describes Dirichlet boundary values by global node index.
type DirichletData struct {
	Nodes  []int
	Values []float64
}

// NewDirichletZero returns Dirichlet data setting the given nodes to zero.
func NewDirichletZero(nodes []int) DirichletData {
	return DirichletData{Nodes: nodes, Values: make([]float64, len(nodes))}
}

// DirichletFromFunc evaluates g at the coordinates of the given 1D nodes to
// build Dirichlet data.
func DirichletFromFunc1D(m *Mesh1D, nodes []int, g func(float64) float64) DirichletData {
	vals := make([]float64, len(nodes))
	for k, n := range nodes {
		vals[k] = g(m.Nodes[n])
	}
	return DirichletData{Nodes: nodes, Values: vals}
}

// DirichletFromFunc2D evaluates g at the coordinates of the given 2D nodes to
// build Dirichlet data.
func DirichletFromFunc2D(m *Mesh2D, nodes []int, g func(x, y float64) float64) DirichletData {
	vals := make([]float64, len(nodes))
	for k, n := range nodes {
		vals[k] = g(m.Nodes[n][0], m.Nodes[n][1])
	}
	return DirichletData{Nodes: nodes, Values: vals}
}

// SolvePoisson1D solves -u” = f on the meshed interval with Dirichlet data,
// using P1 elements. It returns the nodal solution values.
func SolvePoisson1D(m *Mesh1D, f func(float64) float64, bc DirichletData) (Vector, error) {
	A := AssembleStiffness1D(m)
	b := AssembleLoad1D(m, f)
	ApplyDirichlet(A, b, bc.Nodes, bc.Values)
	return SolveSPD(A, b)
}

// SolvePoisson1DP2 solves -u” = f on the meshed interval with Dirichlet data
// using P2 elements. The returned vector has length NumNodes+NumElements with
// midpoint dofs after the vertex dofs (see Mesh1D.P2Nodes).
func SolvePoisson1DP2(m *Mesh1D, f func(float64) float64, bc DirichletData) (Vector, error) {
	A := AssembleStiffness1DP2(m)
	b := AssembleLoad1DP2(m, f)
	ApplyDirichlet(A, b, bc.Nodes, bc.Values)
	return SolveSPD(A, b)
}

// SolveReactionDiffusion1D solves -d*u” + c*u = f on the meshed interval with
// Dirichlet data using P1 elements, for constant coefficients d>0 and c>=0.
func SolveReactionDiffusion1D(m *Mesh1D, d, c float64, f func(float64) float64, bc DirichletData) (Vector, error) {
	A := ScaleSparse(AssembleStiffness1D(m), d)
	if c != 0 {
		A = AddSparse(A, ScaleSparse(AssembleMass1D(m), c))
	}
	b := AssembleLoad1D(m, f)
	ApplyDirichlet(A, b, bc.Nodes, bc.Values)
	return SolveSPD(A, b)
}

// SolvePoisson2D solves -Δu = f on the meshed domain with Dirichlet data using
// P1 elements. It returns the nodal solution values.
func SolvePoisson2D(m *Mesh2D, f func(x, y float64) float64, bc DirichletData) (Vector, error) {
	A := AssembleStiffness2D(m)
	b := AssembleLoad2D(m, f)
	ApplyDirichlet(A, b, bc.Nodes, bc.Values)
	return SolveSPD(A, b)
}

// SolveReactionDiffusion2D solves -d*Δu + c*u = f on the meshed domain with
// Dirichlet data using P1 elements, for constant coefficients d>0 and c>=0.
func SolveReactionDiffusion2D(m *Mesh2D, d, c float64, f func(x, y float64) float64, bc DirichletData) (Vector, error) {
	A := ScaleSparse(AssembleStiffness2D(m), d)
	if c != 0 {
		A = AddSparse(A, ScaleSparse(AssembleMass2D(m), c))
	}
	b := AssembleLoad2D(m, f)
	ApplyDirichlet(A, b, bc.Nodes, bc.Values)
	return SolveSPD(A, b)
}

// SolvePoisson2DNeumann solves -Δu + u = f with the natural Neumann condition
// du/dn = g on the boundary using P1 elements. The reaction term u makes the
// problem well posed without Dirichlet data. It returns the nodal solution.
func SolvePoisson2DNeumann(m *Mesh2D, f func(x, y float64) float64, g func(x, y float64) float64) (Vector, error) {
	A := AddSparse(AssembleStiffness2D(m), AssembleMass2D(m))
	b := AssembleLoad2D(m, f)
	if g != nil {
		ApplyNeumann2D(m, b, m.BoundaryEdges(), g)
	}
	return SolveSPD(A, b)
}

// SolveHelmholtzReal solves the real, sign-definite Helmholtz-type problem
// -Δu + k2*u = f with homogeneous Dirichlet data on the boundary (k2 >= 0).
func SolveHelmholtzReal(m *Mesh2D, k2 float64, f func(x, y float64) float64) (Vector, error) {
	return SolveReactionDiffusion2D(m, 1, k2, f, NewDirichletZero(m.BoundaryNodes()))
}

// ElasticityParams holds the material parameters for isotropic linear
// elasticity: Young's modulus E and Poisson ratio nu, along with a flag for the
// plane-stress (true) or plane-strain (false) assumption.
type ElasticityParams struct {
	E, Nu       float64
	PlaneStress bool
}

// LameParameters returns the Lamé parameters (lambda, mu) for the material,
// respecting the plane-stress / plane-strain convention.
func (p ElasticityParams) LameParameters() (lambda, mu float64) {
	mu = p.E / (2 * (1 + p.Nu))
	if p.PlaneStress {
		lambda = p.E * p.Nu / (1 - p.Nu*p.Nu)
	} else {
		lambda = p.E * p.Nu / ((1 + p.Nu) * (1 - 2*p.Nu))
	}
	return lambda, mu
}

// ConstitutiveMatrix returns the 3×3 plane elasticity constitutive matrix D
// relating stress (sigma_xx, sigma_yy, sigma_xy) to strain
// (eps_xx, eps_yy, gamma_xy).
func (p ElasticityParams) ConstitutiveMatrix() [3][3]float64 {
	if p.PlaneStress {
		c := p.E / (1 - p.Nu*p.Nu)
		return [3][3]float64{
			{c, c * p.Nu, 0},
			{c * p.Nu, c, 0},
			{0, 0, c * (1 - p.Nu) / 2},
		}
	}
	c := p.E / ((1 + p.Nu) * (1 - 2*p.Nu))
	return [3][3]float64{
		{c * (1 - p.Nu), c * p.Nu, 0},
		{c * p.Nu, c * (1 - p.Nu), 0},
		{0, 0, c * (1 - 2*p.Nu) / 2},
	}
}

// ElementStiffnessElasticity returns the 6×6 element stiffness matrix for a
// linear (constant-strain) triangle with the given vertices under plane
// elasticity. Degrees of freedom are ordered (u1x,u1y,u2x,u2y,u3x,u3y).
func ElementStiffnessElasticity(v1, v2, v3 [2]float64, p ElasticityParams) [6][6]float64 {
	grads, area := TriangleGradients(v1, v2, v3)
	a := area
	if a < 0 {
		a = -a
	}
	// Strain-displacement matrix B is 3×6.
	var B [3][6]float64
	for i := 0; i < 3; i++ {
		bx := grads[i][0]
		by := grads[i][1]
		B[0][2*i] = bx
		B[1][2*i+1] = by
		B[2][2*i] = by
		B[2][2*i+1] = bx
	}
	D := p.ConstitutiveMatrix()
	// K = area * B^T D B.
	var DB [3][6]float64
	for r := 0; r < 3; r++ {
		for c := 0; c < 6; c++ {
			var s float64
			for k := 0; k < 3; k++ {
				s += D[r][k] * B[k][c]
			}
			DB[r][c] = s
		}
	}
	var K [6][6]float64
	for i := 0; i < 6; i++ {
		for j := 0; j < 6; j++ {
			var s float64
			for k := 0; k < 3; k++ {
				s += B[k][i] * DB[k][j]
			}
			K[i][j] = a * s
		}
	}
	return K
}

// AssembleElasticity assembles the global 2*N × 2*N stiffness matrix for plane
// linear elasticity on the mesh, with dof 2*n and 2*n+1 being the x and y
// displacement of node n.
func AssembleElasticity(m *Mesh2D, p ElasticityParams) *SparseMatrix {
	dof := 2 * m.NumNodes()
	A := NewSparseMatrix(dof)
	for t := 0; t < m.NumTriangles(); t++ {
		v1, v2, v3 := m.TriangleVertices(t)
		ke := ElementStiffnessElasticity(v1, v2, v3, p)
		tri := m.Triangles[t]
		var gdof [6]int
		for a := 0; a < 3; a++ {
			gdof[2*a] = 2 * tri[a]
			gdof[2*a+1] = 2*tri[a] + 1
		}
		for a := 0; a < 6; a++ {
			for b := 0; b < 6; b++ {
				A.AddEntry(gdof[a], gdof[b], ke[a][b])
			}
		}
	}
	return A
}

// AssembleBodyForce assembles the load vector for a body force (bx,by) per unit
// area on the mesh, returning a vector of length 2*NumNodes.
func AssembleBodyForce(m *Mesh2D, force func(x, y float64) (bx, by float64)) Vector {
	b := make(Vector, 2*m.NumNodes())
	for t := 0; t < m.NumTriangles(); t++ {
		v1, v2, v3 := m.TriangleVertices(t)
		fx := ElementLoadP1Triangle(func(x, y float64) float64 { bx, _ := force(x, y); return bx }, v1, v2, v3, 3)
		fy := ElementLoadP1Triangle(func(x, y float64) float64 { _, by := force(x, y); return by }, v1, v2, v3, 3)
		tri := m.Triangles[t]
		for a := 0; a < 3; a++ {
			b[2*tri[a]] += fx[a]
			b[2*tri[a]+1] += fy[a]
		}
	}
	return b
}

// SolveElasticity2D solves a plane linear elasticity problem on the mesh with a
// body force and Dirichlet displacement constraints (dofs indexed as 2*node+dir
// with dir 0 for x and 1 for y). It returns the nodal displacement vector of
// length 2*NumNodes.
func SolveElasticity2D(m *Mesh2D, p ElasticityParams, force func(x, y float64) (bx, by float64), constraints DirichletData) (Vector, error) {
	A := AssembleElasticity(m, p)
	var b Vector
	if force != nil {
		b = AssembleBodyForce(m, force)
	} else {
		b = make(Vector, 2*m.NumNodes())
	}
	ApplyDirichlet(A, b, constraints.Nodes, constraints.Values)
	return SolveSPD(A, b)
}
