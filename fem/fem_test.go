package fem

import (
	"fmt"
	"math"
	"testing"
)

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func TestGaussLegendreExactness(t *testing.T) {
	// n-point Gauss is exact for polynomials up to degree 2n-1.
	tests := []struct {
		name string
		f    func(float64) float64
		a, b float64
		n    int
		want float64
	}{
		{"x^2 on [0,1]", func(x float64) float64 { return x * x }, 0, 1, 2, 1.0 / 3},
		{"x^3 on [0,1]", func(x float64) float64 { return x * x * x }, 0, 1, 2, 1.0 / 4},
		{"x^5 on [0,2]", func(x float64) float64 { return math.Pow(x, 5) }, 0, 2, 3, 64.0 / 6},
		{"exp on [0,1]", func(x float64) float64 { return math.Exp(x) }, 0, 1, 8, math.E - 1},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := IntegrateInterval(tc.f, tc.a, tc.b, tc.n)
			if !approx(got, tc.want, 1e-10) {
				t.Fatalf("got %v want %v", got, tc.want)
			}
		})
	}
}

func TestTriangleQuadrature(t *testing.T) {
	v1 := [2]float64{0, 0}
	v2 := [2]float64{1, 0}
	v3 := [2]float64{0, 1}
	tests := []struct {
		name   string
		f      func(x, y float64) float64
		degree int
		want   float64
	}{
		{"constant", func(x, y float64) float64 { return 1 }, 1, 0.5},
		{"x", func(x, y float64) float64 { return x }, 2, 1.0 / 6},
		{"xy", func(x, y float64) float64 { return x * y }, 3, 1.0 / 24},
		{"x^2*y^2", func(x, y float64) float64 { return x * x * y * y }, 5, 1.0 / 180},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := IntegrateTriangle(tc.f, v1, v2, v3, tc.degree)
			if !approx(got, tc.want, 1e-12) {
				t.Fatalf("got %v want %v", got, tc.want)
			}
		})
	}
}

func TestElementMatrices1D(t *testing.T) {
	h := 2.0
	k := ElementStiffnessP1Interval(h)
	if !approx(k[0][0], 0.5, 1e-12) || !approx(k[0][1], -0.5, 1e-12) {
		t.Fatalf("P1 stiffness wrong: %v", k)
	}
	m := ElementMassP1Interval(h)
	if !approx(m[0][0], 2*h/6, 1e-12) || !approx(m[0][1], h/6, 1e-12) {
		t.Fatalf("P1 mass wrong: %v", m)
	}
	k2 := ElementStiffnessP2Interval(h)
	if !approx(k2[0][0], 7.0/(3*h), 1e-12) || !approx(k2[2][2], 16.0/(3*h), 1e-12) {
		t.Fatalf("P2 stiffness wrong: %v", k2)
	}
	// Stiffness rows sum to zero (constant is in the kernel).
	for i := 0; i < 3; i++ {
		s := k2[i][0] + k2[i][1] + k2[i][2]
		if !approx(s, 0, 1e-12) {
			t.Fatalf("P2 stiffness row %d does not sum to zero: %v", i, s)
		}
	}
}

func TestElementMassP1TriangleSum(t *testing.T) {
	v1 := [2]float64{0, 0}
	v2 := [2]float64{2, 0}
	v3 := [2]float64{0, 3}
	area := TriangleArea(v1, v2, v3)
	m := ElementMassP1Triangle(v1, v2, v3)
	var sum float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			sum += m[i][j]
		}
	}
	if !approx(sum, area, 1e-12) {
		t.Fatalf("mass sum %v want area %v", sum, area)
	}
}

func TestStiffnessKernel2D(t *testing.T) {
	// Row sums of the element stiffness vanish for both P1 and P2.
	v1 := [2]float64{0.2, 0.1}
	v2 := [2]float64{1.3, 0.4}
	v3 := [2]float64{0.5, 1.2}
	k := ElementStiffnessP1Triangle(v1, v2, v3)
	for i := 0; i < 3; i++ {
		if !approx(k[i][0]+k[i][1]+k[i][2], 0, 1e-12) {
			t.Fatalf("P1 row %d sum nonzero", i)
		}
	}
	k2 := ElementStiffnessP2Triangle(v1, v2, v3)
	for i := 0; i < 6; i++ {
		var s float64
		for j := 0; j < 6; j++ {
			s += k2[i][j]
		}
		if !approx(s, 0, 1e-10) {
			t.Fatalf("P2 row %d sum nonzero: %v", i, s)
		}
	}
	// P2 mass matrix total sum equals the area (partition of unity).
	m2 := ElementMassP2Triangle(v1, v2, v3)
	var msum float64
	for i := 0; i < 6; i++ {
		for j := 0; j < 6; j++ {
			msum += m2[i][j]
		}
	}
	if !approx(msum, TriangleArea(v1, v2, v3), 1e-12) {
		t.Fatalf("P2 mass sum %v want %v", msum, TriangleArea(v1, v2, v3))
	}
}

func TestLUAndCG(t *testing.T) {
	a := MatrixFromRows([][]float64{
		{4, 1, 0},
		{1, 3, 1},
		{0, 1, 2},
	})
	b := VectorOf(1, 2, 3)
	x, err := SolveDense(a, b)
	if err != nil {
		t.Fatal(err)
	}
	r := a.MulVec(x).Sub(b)
	if r.NormInf() > 1e-12 {
		t.Fatalf("LU residual too large: %v", r.NormInf())
	}
	// Same system through the sparse SPD path.
	s := NewSparseMatrix(3)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			s.AddEntry(i, j, a.At(i, j))
		}
	}
	xc, iters, err := ConjugateGradient(s, b, 1e-12, 100)
	if err != nil {
		t.Fatalf("CG failed after %d iters: %v", iters, err)
	}
	if xc.Sub(x).NormInf() > 1e-9 {
		t.Fatalf("CG disagrees with LU: %v vs %v", xc, x)
	}
}

func TestPoisson1DConvergence(t *testing.T) {
	f := func(x float64) float64 { return math.Pi * math.Pi * math.Sin(math.Pi*x) }
	exact := func(x float64) float64 { return math.Sin(math.Pi * x) }
	errs := make([]float64, 0, 2)
	for _, n := range []int{20, 40} {
		mesh := NewUniformMesh1D(0, 1, n)
		bc := NewDirichletZero(mesh.BoundaryNodes())
		u, err := SolvePoisson1D(mesh, f, bc)
		if err != nil {
			t.Fatal(err)
		}
		e := L2Error1D(mesh, u, exact, 5)
		errs = append(errs, e)
	}
	if errs[0] > 1e-2 {
		t.Fatalf("coarse error too big: %v", errs[0])
	}
	rate := errs[0] / errs[1]
	if rate < 3.5 || rate > 4.5 {
		t.Fatalf("L2 convergence rate not ~4: %v", rate)
	}
}

func TestPoisson1DP2(t *testing.T) {
	f := func(x float64) float64 { return math.Pi * math.Pi * math.Sin(math.Pi*x) }
	exact := func(x float64) float64 { return math.Sin(math.Pi * x) }
	mesh := NewUniformMesh1D(0, 1, 10)
	bc := NewDirichletZero(mesh.BoundaryNodes())
	u, err := SolvePoisson1DP2(mesh, f, bc)
	if err != nil {
		t.Fatal(err)
	}
	// Evaluate the P2 solution at the domain midpoint.
	got := EvalP2_1D(mesh, u, 0.5)
	if !approx(got, exact(0.5), 5e-4) {
		t.Fatalf("P2 midpoint value %v want ~%v", got, exact(0.5))
	}
}

func TestReactionDiffusion1D(t *testing.T) {
	// -u'' + u = (1+pi^2) sin(pi x), u(0)=u(1)=0 -> u = sin(pi x).
	d, c := 1.0, 1.0
	f := func(x float64) float64 { return (1 + math.Pi*math.Pi) * math.Sin(math.Pi*x) }
	exact := func(x float64) float64 { return math.Sin(math.Pi * x) }
	mesh := NewUniformMesh1D(0, 1, 40)
	bc := NewDirichletZero(mesh.BoundaryNodes())
	u, err := SolveReactionDiffusion1D(mesh, d, c, f, bc)
	if err != nil {
		t.Fatal(err)
	}
	if e := L2Error1D(mesh, u, exact, 5); e > 5e-3 {
		t.Fatalf("reaction-diffusion L2 error too big: %v", e)
	}
}

func TestPoisson2D(t *testing.T) {
	// -Δu = 2π² sin(πx) sin(πy) with zero BC -> u = sin(πx) sin(πy).
	f := func(x, y float64) float64 {
		return 2 * math.Pi * math.Pi * math.Sin(math.Pi*x) * math.Sin(math.Pi*y)
	}
	exact := func(x, y float64) float64 { return math.Sin(math.Pi*x) * math.Sin(math.Pi*y) }
	mesh := UnitSquareMesh(16)
	bc := NewDirichletZero(mesh.BoundaryNodes())
	u, err := SolvePoisson2D(mesh, f, bc)
	if err != nil {
		t.Fatal(err)
	}
	e := L2Error2D(mesh, u, exact, 4)
	if e > 8e-3 {
		t.Fatalf("Poisson2D L2 error too big: %v", e)
	}
}

func TestPoisson2DConvergence(t *testing.T) {
	f := func(x, y float64) float64 {
		return 2 * math.Pi * math.Pi * math.Sin(math.Pi*x) * math.Sin(math.Pi*y)
	}
	exact := func(x, y float64) float64 { return math.Sin(math.Pi*x) * math.Sin(math.Pi*y) }
	var prev float64
	for i, n := range []int{8, 16} {
		mesh := UnitSquareMesh(n)
		bc := NewDirichletZero(mesh.BoundaryNodes())
		u, err := SolvePoisson2D(mesh, f, bc)
		if err != nil {
			t.Fatal(err)
		}
		e := L2Error2D(mesh, u, exact, 4)
		if i == 1 {
			rate := prev / e
			if rate < 3.0 || rate > 5.0 {
				t.Fatalf("2D L2 rate not ~4: %v", rate)
			}
		}
		prev = e
	}
}

func TestMassTotalArea2D(t *testing.T) {
	mesh := UnitSquareMesh(8)
	M := AssembleMass2D(mesh)
	var sum float64
	for i := 0; i < M.Dim(); i++ {
		for j := 0; j < M.Dim(); j++ {
			sum += M.At(i, j)
		}
	}
	if !approx(sum, 1.0, 1e-10) {
		t.Fatalf("total mass %v want 1", sum)
	}
}

func TestMeshRefine(t *testing.T) {
	m1 := NewUniformMesh1D(0, 1, 4)
	m2 := m1.Refine()
	if m2.NumElements() != 8 {
		t.Fatalf("1D refine elements %d want 8", m2.NumElements())
	}
	tm := UnitSquareMesh(2)
	nt := tm.NumTriangles()
	rt := tm.Refine()
	if rt.NumTriangles() != 4*nt {
		t.Fatalf("2D refine triangles %d want %d", rt.NumTriangles(), 4*nt)
	}
	// Refinement preserves total area.
	if !approx(rt.TotalArea(), tm.TotalArea(), 1e-12) {
		t.Fatalf("refined area changed: %v vs %v", rt.TotalArea(), tm.TotalArea())
	}
}

func TestBoundaryNodes2D(t *testing.T) {
	mesh := UnitSquareMesh(3)
	bn := mesh.BoundaryNodes()
	// A 4x4 node grid has 12 boundary nodes.
	if len(bn) != 12 {
		t.Fatalf("boundary nodes %d want 12", len(bn))
	}
	if len(mesh.InteriorNodes()) != mesh.NumNodes()-12 {
		t.Fatalf("interior node count mismatch")
	}
}

func TestElasticityRigidBody(t *testing.T) {
	mesh := UnitSquareMesh(3)
	p := ElasticityParams{E: 210e9, Nu: 0.3, PlaneStress: true}
	K := AssembleElasticity(mesh, p)
	// Rigid translation in x must produce zero internal force.
	d := make(Vector, 2*mesh.NumNodes())
	for n := 0; n < mesh.NumNodes(); n++ {
		d[2*n] = 1
	}
	force := K.MulVec(d)
	if force.NormInf() > 1e-3*p.E {
		t.Fatalf("rigid translation produced force: %v", force.NormInf())
	}
}

func TestConstitutiveMatrix(t *testing.T) {
	p := ElasticityParams{E: 1, Nu: 0.3, PlaneStress: true}
	D := p.ConstitutiveMatrix()
	c := 1.0 / (1 - 0.09)
	if !approx(D[0][0], c, 1e-12) || !approx(D[0][1], c*0.3, 1e-12) || !approx(D[2][2], c*0.35, 1e-12) {
		t.Fatalf("plane-stress D wrong: %v", D)
	}
	lam, mu := p.LameParameters()
	if !approx(mu, 1.0/(2*1.3), 1e-12) {
		t.Fatalf("mu wrong: %v", mu)
	}
	if lam <= 0 {
		t.Fatalf("lambda should be positive: %v", lam)
	}
}

func TestElasticityBeam(t *testing.T) {
	// Clamp the left edge of a rectangle and pull it in x with a body force;
	// every displacement should be finite and the clamped nodes stay fixed.
	mesh := RectangleMesh(0, 0, 4, 1, 8, 2)
	p := ElasticityParams{E: 1000, Nu: 0.3, PlaneStress: true}
	var clampNodes []int
	var clampVals []float64
	for n, xy := range mesh.Nodes {
		if approx(xy[0], 0, 1e-12) {
			clampNodes = append(clampNodes, 2*n, 2*n+1)
			clampVals = append(clampVals, 0, 0)
		}
	}
	bc := DirichletData{Nodes: clampNodes, Values: clampVals}
	u, err := SolveElasticity2D(mesh, p, func(x, y float64) (float64, float64) { return 1, 0 }, bc)
	if err != nil {
		t.Fatal(err)
	}
	for _, n := range clampNodes {
		if math.Abs(u[n]) > 1e-9 {
			t.Fatalf("clamped dof %d moved: %v", n, u[n])
		}
	}
	// The free (right) end should have moved in the positive x direction.
	var maxUx float64
	for n, xy := range mesh.Nodes {
		if approx(xy[0], 4, 1e-12) && u[2*n] > maxUx {
			maxUx = u[2*n]
		}
	}
	if maxUx <= 0 {
		t.Fatalf("beam did not stretch: %v", maxUx)
	}
}

func TestBarycentricAndShape(t *testing.T) {
	v1 := [2]float64{0, 0}
	v2 := [2]float64{1, 0}
	v3 := [2]float64{0, 1}
	l := Barycentric(v1, v2, v3, [2]float64{0.25, 0.25})
	if !approx(l[0], 0.5, 1e-12) || !approx(l[1], 0.25, 1e-12) || !approx(l[2], 0.25, 1e-12) {
		t.Fatalf("barycentric wrong: %v", l)
	}
	// P2 shape functions form a partition of unity.
	n := ShapeP2Triangle(l)
	var sum float64
	for _, v := range n {
		sum += v
	}
	if !approx(sum, 1, 1e-12) {
		t.Fatalf("P2 partition of unity failed: %v", sum)
	}
}

func TestNeumann2D(t *testing.T) {
	// -Δu + u = f with exact u = x, so f = x and du/dn = n_x on the boundary.
	mesh := UnitSquareMesh(12)
	exact := func(x, y float64) float64 { return x }
	f := func(x, y float64) float64 { return x }
	g := func(x, y float64) float64 {
		// outward normal flux du/dn for u=x on the unit square boundary
		switch {
		case approx(x, 1, 1e-9):
			return 1
		case approx(x, 0, 1e-9):
			return -1
		default:
			return 0
		}
	}
	u, err := SolvePoisson2DNeumann(mesh, f, g)
	if err != nil {
		t.Fatal(err)
	}
	if e := L2Error2D(mesh, u, exact, 4); e > 1e-2 {
		t.Fatalf("Neumann problem error too big: %v", e)
	}
}

func TestP2MeshConnectivity(t *testing.T) {
	mesh := UnitSquareMesh(2)
	nodes, conn := mesh.P2Mesh()
	if len(conn) != mesh.NumTriangles() {
		t.Fatalf("conn length %d want %d", len(conn), mesh.NumTriangles())
	}
	// Every midpoint node coordinate must be the average of its edge endpoints.
	for t0, c := range conn {
		v1, v2, v3 := mesh.TriangleVertices(t0)
		verts := [3][2]float64{v1, v2, v3}
		// midpoint opposite v1 is average of v2,v3 (index 3)
		mid := nodes[c[3]]
		want := [2]float64{0.5 * (verts[1][0] + verts[2][0]), 0.5 * (verts[1][1] + verts[2][1])}
		if !approx(mid[0], want[0], 1e-12) || !approx(mid[1], want[1], 1e-12) {
			t.Fatalf("midpoint mismatch tri %d: %v want %v", t0, mid, want)
		}
	}
}

func ExampleSolvePoisson1D() {
	// Solve -u'' = π² sin(πx) on (0,1) with u(0)=u(1)=0; exact solution sin(πx).
	f := func(x float64) float64 { return math.Pi * math.Pi * math.Sin(math.Pi*x) }
	exact := func(x float64) float64 { return math.Sin(math.Pi * x) }
	mesh := NewUniformMesh1D(0, 1, 50)
	bc := NewDirichletZero(mesh.BoundaryNodes())
	u, _ := SolvePoisson1D(mesh, f, bc)
	err := L2Error1D(mesh, u, exact, 5)
	fmt.Printf("L2 error below 1e-3: %v\n", err < 1e-3)
	// Output: L2 error below 1e-3: true
}

func ExampleIntegrateTriangle() {
	// Integrate f(x,y)=x over the reference triangle; the exact value is 1/6.
	v1 := [2]float64{0, 0}
	v2 := [2]float64{1, 0}
	v3 := [2]float64{0, 1}
	got := IntegrateTriangle(func(x, y float64) float64 { return x }, v1, v2, v3, 2)
	fmt.Printf("%.6f\n", got)
	// Output: 0.166667
}
