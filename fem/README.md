# fem — implements the finite element method (FEM) for elliptic boundary value problems in one

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/fem.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/fem)

Package fem implements the finite element method (FEM) for elliptic boundary
value problems in one and two spatial dimensions.

- Dense and sparse linear algebra (LU factorisation, conjugate gradients). -
Numerical quadrature on intervals (Gauss–Legendre) and on triangles
(symmetric Dunavant/Strang rules up to degree five). - Lagrange P1 and P2
basis functions on the reference interval and the reference triangle, together
with their gradients. - Interval and triangular mesh generation, boundary
extraction and uniform mesh refinement. - Element stiffness, mass and load
matrices and their assembly into global sparse operators. - Application of
Dirichlet, Neumann and Robin boundary conditions. - Ready-made solvers for the
Poisson equation, reaction–diffusion problems and two–dimensional linear
elasticity. - L2 and H1 (semi)norms and the corresponding error estimators.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/fem"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func ApplyDirichlet(A *SparseMatrix, b Vector, nodes []int, values []float64)` | ApplyDirichlet imposes the Dirichlet condition u[nodes`k`] = values`k` on the linear system (A, b) in place using the symmetric elimination technique… |
| `func ApplyNeumann1D(b Vector, node int, flux float64)` | ApplyNeumann1D adds the natural (Neumann) boundary flux value at the given global node to the load vector b in place. |
| `func ApplyNeumann2D(m *Mesh2D, b Vector, edges []Edge, g func(x, y float64) float64)` | ApplyNeumann2D adds the Neumann flux contribution of g over the given boundary edges to the load vector b (P1). |
| `func ApplyRobin1D(A *SparseMatrix, b Vector, node int, alpha, g float64)` | ApplyRobin1D imposes a Robin condition alpha*u = g (contributing to the diagonal) at the given global node by adding alpha to A[node,node] and g to… |
| `func ApplyRobin2D(m *Mesh2D, A *SparseMatrix, b Vector, edges []Edge, alpha float64, g func(x, y float64) float64)` | ApplyRobin2D adds Robin contributions alpha*u = g over the given boundary edges, augmenting both the stiffness matrix A (edge mass scaled by alpha)… |
| `func Barycentric(v1, v2, v3, p [2]float64) [3]float64` | Barycentric returns the barycentric coordinates (L1,L2,L3) of the point p with respect to the triangle with vertices v1,v2,v3. |
| `func Determinant(a *Matrix) float64` | Determinant returns the determinant of a square matrix, or 0 if it is singular. |
| `func ElementLoadP1Interval(f func(float64) float64, x0, x1 float64, gauss int) [2]float64` | ElementLoadP1Interval returns the 2-vector load contribution of the source f over the interval element [x0,x1], integrated with a Gauss rule of the… |
| `func ElementLoadP1Triangle(f func(x, y float64) float64, v1, v2, v3 [2]float64, degree int) [3]float64` | ElementLoadP1Triangle returns the 3-vector load contribution of f over the triangle, integrated with a quadrature rule of the given degree. |
| `func ElementLoadP2Interval(f func(float64) float64, x0, x1 float64, gauss int) [3]float64` | ElementLoadP2Interval returns the 3-vector load contribution of f over the interval element [x0,x1] in the node ordering (left, right, midpoint). |
| `func ElementLoadP2Triangle(f func(x, y float64) float64, v1, v2, v3 [2]float64) [6]float64` | ElementLoadP2Triangle returns the 6-vector load contribution of f over the triangle, integrated with a degree-5 quadrature rule. |
| `func ElementLumpedMassP1Triangle(v1, v2, v3 [2]float64) [3]float64` | ElementLumpedMassP1Triangle returns the diagonal (row-sum lumped) P1 mass matrix for the triangle, each diagonal entry being area/3. |
| `func ElementMassP1Interval(h float64) [2][2]float64` | ElementMassP1Interval returns the 2×2 P1 mass matrix for an interval element of length h, the integral of products of basis functions. |
| `func ElementMassP1Triangle(v1, v2, v3 [2]float64) [3][3]float64` | ElementMassP1Triangle returns the 3×3 P1 (consistent) mass matrix for the triangle with the given vertices. |
| `func ElementMassP2Interval(h float64) [3][3]float64` | ElementMassP2Interval returns the 3×3 P2 mass matrix for an interval element of length h in the node ordering (left, right, midpoint). |
| `func ElementMassP2Triangle(v1, v2, v3 [2]float64) [6][6]float64` | ElementMassP2Triangle returns the 6×6 P2 mass matrix for the triangle with the given vertices, computed by quadrature. |
| `func ElementStiffnessElasticity(v1, v2, v3 [2]float64, p ElasticityParams) [6][6]float64` | ElementStiffnessElasticity returns the 6×6 element stiffness matrix for a linear (constant-strain) triangle with the given vertices under plane… |
| `func ElementStiffnessP1Interval(h float64) [2][2]float64` | ElementStiffnessP1Interval returns the 2×2 P1 stiffness matrix for an interval element of length h, the integral of the product of basis-function… |
| `func ElementStiffnessP1Triangle(v1, v2, v3 [2]float64) [3][3]float64` | ElementStiffnessP1Triangle returns the 3×3 P1 stiffness matrix for the triangle with the given vertices. |
| `func ElementStiffnessP2Interval(h float64) [3][3]float64` | ElementStiffnessP2Interval returns the 3×3 P2 stiffness matrix for an interval element of length h in the node ordering (left, right, midpoint). |
| `func ElementStiffnessP2Triangle(v1, v2, v3 [2]float64) [6][6]float64` | ElementStiffnessP2Triangle returns the 6×6 P2 stiffness matrix for the triangle with the given vertices, computed by quadrature. |
| `func EvalP1_1D(m *Mesh1D, u Vector, x float64) float64` | EvalP1_1D evaluates the P1 finite element function with the given nodal values at coordinate x. |
| `func EvalP2_1D(m *Mesh1D, u Vector, x float64) float64` | EvalP2_1D evaluates the P2 finite element function with the given dof values (as returned by SolvePoisson1DP2) at coordinate x. |
| `func GradShapeP1Interval() [2]float64` | GradShapeP1Interval returns the derivatives d/dxi of the linear basis functions on the reference interval. |
| `func GradShapeP1TriangleRef() [3][2]float64` | GradShapeP1TriangleRef returns the reference gradients (d/dxi, d/deta) of the linear triangle basis functions on the reference triangle. |
| `func GradShapeP2Interval(xi float64) [3]float64` | GradShapeP2Interval returns the derivatives d/dxi of the quadratic basis functions on the reference interval at xi. |
| `func GradShapeP2TriangleRef(xi, eta float64) [6][2]float64` | GradShapeP2TriangleRef returns the reference gradients (d/dxi, d/deta) of the six quadratic basis functions, evaluated at reference coordinates… |
| `func H1Error1D(m *Mesh1D, u Vector, exact, exactDeriv func(float64) float64, gauss int) float64` | H1Error1D returns the full H1 norm of the error, sqrt(L2^2 + seminorm^2). |
| `func H1Error2D(m *Mesh2D, u Vector, exact func(x, y float64) float64, exactGrad func(x, y float64) (gx, gy float64), degree int) float64` | H1Error2D returns the full H1 norm of the error, sqrt(L2^2 + seminorm^2). |
| `func H1SeminormError1D(m *Mesh1D, u Vector, exactDeriv func(float64) float64, gauss int) float64` | H1SeminormError1D returns the H1 seminorm of the error (the L2 norm of the derivative error) for the P1 solution u. |
| `func H1SeminormError2D(m *Mesh2D, u Vector, exactGrad func(x, y float64) (gx, gy float64), degree int) float64` | H1SeminormError2D returns the H1 seminorm of the error (the L2 norm of the gradient error) for the P1 solution u, given the exact gradient. |
| `func IntegrateInterval(f func(float64) float64, a, b float64, n int) float64` | IntegrateInterval approximates the integral of f over [a, b] using the n-point Gauss–Legendre rule. |
| `func IntegrateTriangle(f func(x, y float64) float64, v1, v2, v3 [2]float64, degree int) float64` | IntegrateTriangle approximates the integral of f over the triangle with the given vertices using a rule of the requested degree. |
| `func L2Error1D(m *Mesh1D, u Vector, exact func(float64) float64, gauss int) float64` | L2Error1D returns the L2 norm of the difference between the P1 solution u and the exact function exact, integrated element by element. |
| `func L2Error2D(m *Mesh2D, u Vector, exact func(x, y float64) float64, degree int) float64` | L2Error2D returns the L2 norm of the difference between the P1 solution u and the exact function exact on the triangular mesh. |
| `func L2Norm1D(m *Mesh1D, f func(float64) float64, gauss int) float64` | L2Norm1D returns the L2 norm of a function over the meshed interval. |
| `func L2Norm2D(m *Mesh2D, f func(x, y float64) float64, degree int) float64` | L2Norm2D returns the L2 norm of a function over the triangular mesh. |
| `func PhysicalGradP2Triangle(v1, v2, v3 [2]float64, xi, eta float64) [6][2]float64` | PhysicalGradP2Triangle maps the reference gradients of the P2 basis functions at reference coordinates (xi,eta) to physical gradients on the triangle… |
| `func ShapeP1Interval(xi float64) [2]float64` | ShapeP1Interval returns the two linear Lagrange basis functions evaluated at the reference coordinate xi in [0,1]. |
| `func ShapeP1Triangle(l [3]float64) [3]float64` | ShapeP1Triangle returns the three linear basis functions given the barycentric coordinates l = (L1,L2,L3). |
| `func ShapeP2Interval(xi float64) [3]float64` | ShapeP2Interval returns the three quadratic Lagrange basis functions on the reference interval [0,1]. |
| `func ShapeP2Triangle(l [3]float64) [6]float64` | ShapeP2Triangle returns the six quadratic basis functions for barycentric coordinates l = (L1,L2,L3). |
| `func TriangleArea(v1, v2, v3 [2]float64) float64` | TriangleArea returns the (unsigned) area of the triangle with the given vertices. |
| `func TriangleGradients(v1, v2, v3 [2]float64) (grads [3][2]float64, area float64)` | TriangleGradients returns the constant physical gradients of the three barycentric coordinates (equivalently the P1 basis functions) on the triangle… |
| `func TriangleJacobian(v1, v2, v3 [2]float64) (jac [2][2]float64, det float64, invT [2][2]float64)` | TriangleJacobian returns the 2×2 Jacobian of the affine map from the reference triangle to the physical triangle, its determinant and its inverse… |
| `func TriangleSignedArea(v1, v2, v3 [2]float64) float64` | TriangleSignedArea returns the signed area of the triangle, positive when the vertices are counter-clockwise. |

### Types

| Type | What it is |
| --- | --- |
| `DirichletData` | DirichletData describes Dirichlet boundary values by global node index. |
| `Edge` | Edge is an undirected mesh edge identified by its two node indices with the smaller index stored first. |
| `ElasticityParams` | ElasticityParams holds the material parameters for isotropic linear elasticity: Young's modulus E and Poisson ratio nu, along with a flag for the… |
| `LU` | LU holds an in-place LU factorisation with partial pivoting. |
| `Matrix` | Matrix is a dense row-major real matrix. |
| `Mesh1D` | Mesh1D is a one-dimensional mesh of an interval, described by sorted node coordinates. |
| `Mesh2D` | Mesh2D is a two-dimensional triangular mesh: a list of node coordinates and a list of triangles given by triples of node indices. |
| `QuadRule1D` | QuadRule1D is a quadrature rule on the reference interval [-1, 1]. |
| `SparseMatrix` | SparseMatrix is a square sparse matrix in dictionary-of-keys form. |
| `TriQuadRule` | TriQuadRule is a symmetric quadrature rule on a triangle. |
| `Vector` | Vector is a dense real vector stored as a slice of float64 values. |

<details>
<summary><code>DirichletData</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DirichletFromFunc1D(m *Mesh1D, nodes []int, g func(float64) float64) DirichletData` | DirichletFromFunc1D evaluates g at the coordinates of the given 1D nodes to build Dirichlet data. |
| `func DirichletFromFunc2D(m *Mesh2D, nodes []int, g func(x, y float64) float64) DirichletData` | DirichletFromFunc2D evaluates g at the coordinates of the given 2D nodes to build Dirichlet data. |
| `func NewDirichletZero(nodes []int) DirichletData` | NewDirichletZero returns Dirichlet data setting the given nodes to zero. |

</details>

<details>
<summary><code>Edge</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func MakeEdge(a, b int) Edge` | MakeEdge returns the canonical Edge for nodes a and b (smaller index first). |

</details>

<details>
<summary><code>ElasticityParams</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (p ElasticityParams) ConstitutiveMatrix() [3][3]float64` | ConstitutiveMatrix returns the 3×3 plane elasticity constitutive matrix D relating stress (sigma_xx, sigma_yy, sigma_xy) to strain (eps_xx, eps_yy,… |
| `func (p ElasticityParams) LameParameters() (lambda, mu float64)` | LameParameters returns the Lamé parameters (lambda, mu) for the material, respecting the plane-stress / plane-strain convention. |

</details>

<details>
<summary><code>LU</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LUDecompose(a *Matrix) (*LU, error)` | LUDecompose computes the LU factorisation of a square matrix using partial pivoting. |
| `func (f *LU) Det() float64` | Det returns the determinant of the factorised matrix. |
| `func (f *LU) Solve(b Vector) Vector` | Solve solves A x = b using the stored factorisation. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Identity(n int) *Matrix` | Identity returns the n×n identity matrix. |
| `func Inverse(a *Matrix) (*Matrix, error)` | Inverse returns the inverse of a square matrix, or an error if it is singular. |
| `func MatrixFromRows(rows [][]float64) *Matrix` | MatrixFromRows builds a matrix from a slice of equal-length rows. |
| `func NewMatrix(r, c int) *Matrix` | NewMatrix returns a zero r×c matrix. |
| `func Zeros(r, c int) *Matrix` | Zeros is an alias for NewMatrix returning an r×c zero matrix. |
| `func (m *Matrix) Add(i, j int, v float64)` | Add adds v to the entry at row i, column j (scatter-add used in assembly). |
| `func (m *Matrix) AddMatrix(b *Matrix) *Matrix` | AddMatrix returns the sum m+b. |
| `func (m *Matrix) At(i, j int) float64` | At returns the entry at row i, column j. |
| `func (m *Matrix) Clone() *Matrix` | Clone returns an independent copy of the matrix. |
| `func (m *Matrix) Col(j int) Vector` | Col returns a copy of column j as a Vector. |
| `func (m *Matrix) Cols() int` | Cols returns the number of columns. |
| `func (m *Matrix) IsSymmetric(tol float64) bool` | IsSymmetric reports whether the matrix is square and symmetric within tol. |
| `func (m *Matrix) Mul(b *Matrix) *Matrix` | Mul returns the matrix product m*b. |
| `func (m *Matrix) MulVec(x Vector) Vector` | MulVec returns the matrix–vector product m*x. |
| `func (m *Matrix) Row(i int) Vector` | Row returns a copy of row i as a Vector. |
| `func (m *Matrix) Rows() int` | Rows returns the number of rows. |
| `func (m *Matrix) Scale(a float64) *Matrix` | Scale returns a*m as a new matrix. |
| `func (m *Matrix) Set(i, j int, v float64)` | Set assigns v to the entry at row i, column j. |
| `func (m *Matrix) Transpose() *Matrix` | Transpose returns the transpose of the matrix. |

</details>

<details>
<summary><code>Mesh1D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GradedMesh1D(a, b float64, n int, exponent float64) *Mesh1D` | GradedMesh1D returns a mesh of [a,b] with n elements whose node spacing is graded towards a by the power exponent (exponent 1 is uniform, exponents >… |
| `func NewMesh1D(nodes []float64) (*Mesh1D, error)` | NewMesh1D builds a mesh from the given node coordinates. |
| `func NewUniformMesh1D(a, b float64, n int) *Mesh1D` | NewUniformMesh1D returns a uniform mesh of the interval [a,b] with n elements (n+1 equally spaced nodes). |
| `func (m *Mesh1D) BoundaryNodes() []int` | BoundaryNodes returns the indices of the two boundary nodes (the endpoints). |
| `func (m *Mesh1D) Domain() (a, b float64)` | Domain returns the endpoints of the meshed interval. |
| `func (m *Mesh1D) ElementLength(e int) float64` | ElementLength returns the length of element e. |
| `func (m *Mesh1D) ElementMidpoint(e int) float64` | ElementMidpoint returns the midpoint coordinate of element e. |
| `func (m *Mesh1D) ElementNodes(e int) (int, int)` | ElementNodes returns the global node indices (left, right) of element e. |
| `func (m *Mesh1D) InteriorNodes() []int` | InteriorNodes returns the indices of the interior nodes. |
| `func (m *Mesh1D) MaxElementLength() float64` | MaxElementLength returns the length of the longest element (the mesh size h). |
| `func (m *Mesh1D) MinElementLength() float64` | MinElementLength returns the length of the shortest element. |
| `func (m *Mesh1D) NumElements() int` | NumElements returns the number of elements. |
| `func (m *Mesh1D) NumNodes() int` | NumNodes returns the number of nodes. |
| `func (m *Mesh1D) P2Connectivity() [][3]int` | P2Connectivity returns, for each element, the three global P2 degree-of- freedom indices (left node, right node, midpoint) indexing into P2Nodes. |
| `func (m *Mesh1D) P2Nodes() []float64` | P2Nodes returns the coordinates of all P2 degrees of freedom: the mesh nodes followed by the element midpoints. |
| `func (m *Mesh1D) Refine() *Mesh1D` | Refine returns a new mesh obtained by bisecting every element of m. |
| `func (m *Mesh1D) RefineN(k int) *Mesh1D` | RefineN applies Refine k times. |

</details>

<details>
<summary><code>Mesh2D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewMesh2D(nodes [][2]float64, triangles [][3]int) (*Mesh2D, error)` | NewMesh2D builds and validates a triangular mesh from nodes and triangles. |
| `func RectangleMesh(x0, y0, x1, y1 float64, nx, ny int) *Mesh2D` | RectangleMesh returns a structured triangular mesh of the rectangle [x0,x1]×[y0,y1] using nx by ny cells, each split into two triangles. |
| `func UnitSquareMesh(n int) *Mesh2D` | UnitSquareMesh returns a structured mesh of the unit square with n cells per side. |
| `func (m *Mesh2D) AllEdges() []Edge` | AllEdges returns every distinct edge of the mesh, sorted. |
| `func (m *Mesh2D) BoundaryEdges() []Edge` | BoundaryEdges returns the edges that belong to exactly one triangle, sorted for reproducibility. |
| `func (m *Mesh2D) BoundaryNodes() []int` | BoundaryNodes returns the sorted indices of nodes lying on the boundary. |
| `func (m *Mesh2D) Centroid(t int) [2]float64` | Centroid returns the centroid coordinate of triangle t. |
| `func (m *Mesh2D) EnsureCounterClockwise() *Mesh2D` | EnsureCounterClockwise returns a copy of the mesh with every triangle re-ordered to have counter-clockwise (positive) orientation. |
| `func (m *Mesh2D) InteriorNodes() []int` | InteriorNodes returns the sorted indices of nodes not on the boundary. |
| `func (m *Mesh2D) MaxDiameter() float64` | MaxDiameter returns the largest triangle diameter (longest edge length), a common mesh-size measure. |
| `func (m *Mesh2D) NumNodes() int` | NumNodes returns the number of nodes. |
| `func (m *Mesh2D) NumTriangles() int` | NumTriangles returns the number of triangles. |
| `func (m *Mesh2D) P2Mesh() (nodes [][2]float64, conn [][6]int)` | P2Mesh converts the P1 triangular mesh into P2 data: coordinates of all P2 degrees of freedom (vertices followed by edge midpoints) and, for each… |
| `func (m *Mesh2D) Refine() *Mesh2D` | Refine returns the uniformly (red) refined mesh: each triangle is split into four by connecting its edge midpoints. |
| `func (m *Mesh2D) RefineN(k int) *Mesh2D` | RefineN applies Refine k times. |
| `func (m *Mesh2D) TotalArea() float64` | TotalArea returns the sum of all triangle areas (the meshed domain area). |
| `func (m *Mesh2D) TriangleArea(t int) float64` | TriangleArea returns the unsigned area of triangle t. |
| `func (m *Mesh2D) TriangleVertices(t int) (v1, v2, v3 [2]float64)` | TriangleVertices returns the three vertex coordinates of triangle t. |

</details>

<details>
<summary><code>QuadRule1D</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GaussLegendre1D(n int) QuadRule1D` | GaussLegendre1D returns the n-point Gauss–Legendre rule on [-1, 1], exact for polynomials of degree up to 2n-1. |
| `func (q QuadRule1D) MapToInterval(a, b float64) (nodes, weights []float64)` | MapToInterval returns nodes and weights of the rule mapped from [-1,1] to the interval [a, b]. |

</details>

<details>
<summary><code>SparseMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AddSparse(a, b *SparseMatrix) *SparseMatrix` | AddSparse returns a+b for two sparse matrices of equal dimension. |
| `func AssembleElasticity(m *Mesh2D, p ElasticityParams) *SparseMatrix` | AssembleElasticity assembles the global 2*N × 2*N stiffness matrix for plane linear elasticity on the mesh, with dof 2*n and 2*n+1 being the x and y… |
| `func AssembleMass1D(m *Mesh1D) *SparseMatrix` | AssembleMass1D assembles the global P1 mass matrix for the mesh. |
| `func AssembleMass1DP2(m *Mesh1D) *SparseMatrix` | AssembleMass1DP2 assembles the global P2 mass matrix. |
| `func AssembleMass2D(m *Mesh2D) *SparseMatrix` | AssembleMass2D assembles the global P1 mass matrix for the mesh. |
| `func AssembleReaction2D(m *Mesh2D, kappa float64) *SparseMatrix` | AssembleReaction2D assembles kappa*M, the reaction (zeroth-order) contribution with constant coefficient kappa, into a new sparse matrix. |
| `func AssembleStiffness1D(m *Mesh1D) *SparseMatrix` | AssembleStiffness1D assembles the global P1 stiffness matrix for the mesh. |
| `func AssembleStiffness1DP2(m *Mesh1D) *SparseMatrix` | AssembleStiffness1DP2 assembles the global P2 stiffness matrix. |
| `func AssembleStiffness2D(m *Mesh2D) *SparseMatrix` | AssembleStiffness2D assembles the global P1 stiffness matrix for the mesh. |
| `func NewSparseMatrix(n int) *SparseMatrix` | NewSparseMatrix returns an n×n zero sparse matrix. |
| `func ScaleSparse(a *SparseMatrix, s float64) *SparseMatrix` | ScaleSparse returns s*a for scalar s. |
| `func (s *SparseMatrix) AddEntry(i, j int, v float64)` | AddEntry scatters v into entry (i,j), the fundamental assembly operation. |
| `func (s *SparseMatrix) At(i, j int) float64` | At returns the entry at row i, column j (zero if absent). |
| `func (s *SparseMatrix) ClearRow(i int)` | ClearRow removes every entry of row i. |
| `func (s *SparseMatrix) Diagonal() Vector` | Diagonal returns the main diagonal of the sparse matrix as a Vector. |
| `func (s *SparseMatrix) Dim() int` | Dim returns the dimension n of the n×n matrix. |
| `func (s *SparseMatrix) IsSymmetric(tol float64) bool` | IsSymmetric reports whether the sparse matrix is symmetric within tol. |
| `func (s *SparseMatrix) MulVec(x Vector) Vector` | MulVec returns the sparse matrix–vector product s*x. |
| `func (s *SparseMatrix) NNZ() int` | NNZ returns the number of stored (structurally non-zero) entries. |
| `func (s *SparseMatrix) ScaleRow(i int, a float64)` | ScaleRow multiplies every stored entry of row i by a. |
| `func (s *SparseMatrix) Set(i, j int, v float64)` | Set assigns v to entry (i,j), removing the entry if v is exactly zero. |
| `func (s *SparseMatrix) ToDense() *Matrix` | ToDense returns the dense representation of the sparse matrix. |

</details>

<details>
<summary><code>TriQuadRule</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func TriangleQuadrature(degree int) TriQuadRule` | TriangleQuadrature returns a symmetric quadrature rule on the reference triangle exact for polynomials up to the requested degree (1, 2, 3 or 5; a… |

</details>

<details>
<summary><code>Vector</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AXPY(a float64, x, y Vector) Vector` | AXPY returns a*x+y, the classic BLAS operation. |
| `func AssembleBodyForce(m *Mesh2D, force func(x, y float64) (bx, by float64)) Vector` | AssembleBodyForce assembles the load vector for a body force (bx,by) per unit area on the mesh, returning a vector of length 2*NumNodes. |
| `func AssembleLoad1D(m *Mesh1D, f func(float64) float64) Vector` | AssembleLoad1D assembles the global P1 load vector for source f. |
| `func AssembleLoad1DP2(m *Mesh1D, f func(float64) float64) Vector` | AssembleLoad1DP2 assembles the global P2 load vector for source f. |
| `func AssembleLoad2D(m *Mesh2D, f func(x, y float64) float64) Vector` | AssembleLoad2D assembles the global P1 load vector for source f. |
| `func ConjugateGradient(a *SparseMatrix, b Vector, tol float64, maxIter int) (Vector, int, error)` | ConjugateGradient solves the symmetric positive-definite system A x = b using the (unpreconditioned) conjugate gradient method. |
| `func InterpolateNodal1D(m *Mesh1D, f func(float64) float64) Vector` | InterpolateNodal1D returns the nodal interpolant of f on the mesh (its value at each node). |
| `func InterpolateNodal2D(m *Mesh2D, f func(x, y float64) float64) Vector` | InterpolateNodal2D returns the nodal interpolant of f on the mesh. |
| `func NewVector(n int) Vector` | NewVector returns a zero Vector of length n. |
| `func PCGJacobi(a *SparseMatrix, b Vector, tol float64, maxIter int) (Vector, int, error)` | PCGJacobi solves A x = b with Jacobi (diagonal) preconditioned conjugate gradients. |
| `func SolveDense(a *Matrix, b Vector) (Vector, error)` | SolveDense solves the dense linear system A x = b by LU factorisation. |
| `func SolveElasticity2D(m *Mesh2D, p ElasticityParams, force func(x, y float64) (bx, by float64), constraints DirichletData) (Vector, error)` | SolveElasticity2D solves a plane linear elasticity problem on the mesh with a body force and Dirichlet displacement constraints (dofs indexed as… |
| `func SolveHelmholtzReal(m *Mesh2D, k2 float64, f func(x, y float64) float64) (Vector, error)` | SolveHelmholtzReal solves the real, sign-definite Helmholtz-type problem -Δu + k2*u = f with homogeneous Dirichlet data on the boundary (k2 >= 0). |
| `func SolvePoisson1D(m *Mesh1D, f func(float64) float64, bc DirichletData) (Vector, error)` | SolvePoisson1D solves -u” = f on the meshed interval with Dirichlet data, using P1 elements. |
| `func SolvePoisson1DP2(m *Mesh1D, f func(float64) float64, bc DirichletData) (Vector, error)` | SolvePoisson1DP2 solves -u” = f on the meshed interval with Dirichlet data using P2 elements. |
| `func SolvePoisson2D(m *Mesh2D, f func(x, y float64) float64, bc DirichletData) (Vector, error)` | SolvePoisson2D solves -Δu = f on the meshed domain with Dirichlet data using P1 elements. |
| `func SolvePoisson2DNeumann(m *Mesh2D, f func(x, y float64) float64, g func(x, y float64) float64) (Vector, error)` | SolvePoisson2DNeumann solves -Δu + u = f with the natural Neumann condition du/dn = g on the boundary using P1 elements. |
| `func SolveReactionDiffusion1D(m *Mesh1D, d, c float64, f func(float64) float64, bc DirichletData) (Vector, error)` | SolveReactionDiffusion1D solves -d*u” + c*u = f on the meshed interval with Dirichlet data using P1 elements, for constant coefficients d>0 and… |
| `func SolveReactionDiffusion2D(m *Mesh2D, d, c float64, f func(x, y float64) float64, bc DirichletData) (Vector, error)` | SolveReactionDiffusion2D solves -d*Δu + c*u = f on the meshed domain with Dirichlet data using P1 elements, for constant coefficients d>0 and c>=0. |
| `func SolveSPD(a *SparseMatrix, b Vector) (Vector, error)` | SolveSPD solves a symmetric positive-definite sparse system, preferring the dense LU path for small problems and Jacobi-preconditioned CG otherwise. |
| `func VectorOf(vals ...float64) Vector` | VectorOf returns a Vector containing the supplied values. |
| `func (v Vector) Add(w Vector) Vector` | Add returns the elementwise sum v+w. |
| `func (v Vector) Clone() Vector` | Clone returns an independent copy of the vector. |
| `func (v Vector) Dot(w Vector) float64` | Dot returns the Euclidean inner product of v and w. |
| `func (v Vector) Len() int` | Len returns the number of entries in the vector. |
| `func (v Vector) Norm2() float64` | Norm2 returns the Euclidean (L2) norm of the vector. |
| `func (v Vector) NormInf() float64` | NormInf returns the maximum absolute value of the vector entries. |
| `func (v Vector) Scale(a float64) Vector` | Scale returns a new vector equal to a*v. |
| `func (v Vector) Sub(w Vector) Vector` | Sub returns the elementwise difference v-w. |
| `func (v Vector) Sum() float64` | Sum returns the sum of the vector entries. |

</details>

### Variables

`ErrNotConverged`, `ErrSingular`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/fem).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
