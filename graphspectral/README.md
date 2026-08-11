# graphspectral — implements spectral graph theory using only the Go standard library

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/graphspectral.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/graphspectral)

Package graphspectral implements spectral graph theory using only the Go
standard library.

The package centres on two data types. `Matrix` is a small dense real matrix
with the linear-algebra operations required for spectral analysis: addition,
multiplication, transposition, determinants, linear solves, inverses,
symmetric eigendecomposition (via the cyclic Jacobi method) and the leading
eigenpair by power iteration. `Graph` is a finite undirected, optionally
weighted, graph stored as a symmetric adjacency matrix; it produces the
adjacency, degree, (combinatorial) Laplacian, normalized Laplacian,
random-walk Laplacian, signless Laplacian, transition and incidence matrices,
and answers structural questions such as connectivity, bipartiteness,
regularity and density.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/graphspectral"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AdjacencySpectrum(g *Graph) ([]float64, error)` | AdjacencySpectrum returns the eigenvalues of the adjacency matrix in ascending order. |
| `func AlgebraicConnectivity(g *Graph) (float64, error)` | AlgebraicConnectivity returns the second-smallest eigenvalue of the Laplacian (the Fiedler value). |
| `func ArgMax(v []float64) int` | ArgMax returns the index of the (first) maximum entry of v, or -1 if v is empty. |
| `func ArgMin(v []float64) int` | ArgMin returns the index of the (first) minimum entry of v, or -1 if v is empty. |
| `func BFSDistances(g *Graph, source int) []float64` | BFSDistances returns the number of edges on a shortest (unweighted) path from source to every vertex, ignoring edge weights. |
| `func CheegerConstant(g *Graph) (float64, error)` | CheegerConstant returns the exact Cheeger (conductance) constant, the minimum conductance over all nonempty proper vertex subsets, computed by brute… |
| `func CheegerLowerBound(g *Graph) (float64, error)` | CheegerLowerBound returns λ/2, the lower Cheeger bound on the conductance constant, where λ is the second-smallest eigenvalue of the normalized… |
| `func CheegerUpperBound(g *Graph) (float64, error)` | CheegerUpperBound returns sqrt(2λ), the upper Cheeger bound on the conductance constant, where λ is the second-smallest eigenvalue of the… |
| `func ClosenessCentrality(g *Graph) []float64` | ClosenessCentrality returns the closeness centrality of every vertex using weighted shortest-path distances. |
| `func Conductance(g *Graph, s map[int]bool) float64` | Conductance returns the conductance of the cut (S, V\S): \|∂S\| / min(vol(S), vol(V\S)). |
| `func CosineSimilarity(a, b []float64) float64` | CosineSimilarity returns the cosine of the angle between a and b. |
| `func CutSize(g *Graph, part []int) float64` | CutSize returns the total weight of edges crossing between the two blocks of a bipartition given by part (entries 0 or 1). |
| `func DegreeCentrality(g *Graph) []float64` | DegreeCentrality returns the degree centrality of every vertex: its weighted degree divided by n-1, so that in an unweighted graph a vertex adjacent… |
| `func Determinant(m *Matrix) (float64, error)` | Determinant returns the determinant of a square matrix, computed by Gaussian elimination with partial pivoting. |
| `func Diameter(g *Graph) float64` | Diameter returns the greatest weighted shortest-path distance between any two vertices of a connected graph. |
| `func DijkstraDistances(g *Graph, source int) []float64` | DijkstraDistances returns the weighted shortest-path distance from source to every vertex, treating each edge weight as a positive length. |
| `func Dot(a, b []float64) float64` | Dot returns the Euclidean inner product of a and b. |
| `func EdgeBoundary(g *Graph, s map[int]bool) float64` | EdgeBoundary returns the total weight of edges with exactly one endpoint in the vertex set S (given as a membership map). |
| `func EffectiveResistance(g *Graph, s, t int) (float64, error)` | EffectiveResistance returns the effective (electrical) resistance between vertices s and t, treating each edge as a unit conductance times its weight. |
| `func Eigenvalues(m *Matrix) ([]float64, error)` | Eigenvalues returns the eigenvalues of a symmetric matrix in ascending order. |
| `func EigenvectorCentrality(g *Graph, maxIter int, tol float64) ([]float64, error)` | EigenvectorCentrality returns the eigenvector centrality: the (non-negative) dominant eigenvector of the adjacency matrix, scaled to unit Euclidean… |
| `func EstradaIndex(g *Graph) (float64, error)` | EstradaIndex returns the Estrada index Σ exp(λ_i) over the adjacency eigenvalues. |
| `func FiedlerValue(g *Graph) (float64, error)` | FiedlerValue is a synonym for `AlgebraicConnectivity`. |
| `func FiedlerVector(g *Graph) ([]float64, error)` | FiedlerVector returns the eigenvector of the Laplacian associated with the second-smallest eigenvalue. |
| `func GraphEnergy(g *Graph) (float64, error)` | GraphEnergy returns the graph energy: the sum of the absolute values of the adjacency eigenvalues. |
| `func HarmonicCentrality(g *Graph) []float64` | HarmonicCentrality returns the harmonic centrality of every vertex, Σ_{u≠v} 1/d(v,u) divided by n-1, using weighted shortest-path distances. |
| `func KMeans(points [][]float64, k, maxIter int) []int` | KMeans runs Lloyd's algorithm on the given points (each a fixed-length coordinate slice) with k clusters and at most maxIter iterations, returning a… |
| `func KatzCentrality(g *Graph, alpha, beta float64) ([]float64, error)` | KatzCentrality returns the Katz centrality x = (I - α·A)^{-1}·β·1, solving the corresponding linear system. |
| `func KirchhoffIndex(g *Graph) (float64, error)` | KirchhoffIndex returns the Kirchhoff index: the sum of effective resistances over all unordered vertex pairs, equal to n·Σ 1/μ_i over the nonzero… |
| `func LaplacianEnergy(g *Graph) (float64, error)` | LaplacianEnergy returns the Laplacian energy Σ \|μ_i - 2m/n\|, where the μ_i are the Laplacian eigenvalues and 2m/n is the average (weighted)… |
| `func LaplacianSpectralGap(g *Graph) (float64, error)` | LaplacianSpectralGap returns the difference between the second-smallest and smallest Laplacian eigenvalues. |
| `func LaplacianSpectrum(g *Graph) ([]float64, error)` | LaplacianSpectrum returns the eigenvalues of the combinatorial Laplacian in ascending order. |
| `func LargestEigenpair(a *Matrix, maxIter int, tol float64) (float64, []float64, error)` | LargestEigenpair returns the algebraically largest eigenvalue of a real symmetric matrix together with a unit eigenvector, using the power method on… |
| `func Norm1(v []float64) float64` | Norm1 returns the L1 (taxicab) norm of v, the sum of absolute values. |
| `func Norm2(v []float64) float64` | Norm2 returns the Euclidean (L2) norm of v. |
| `func NormInf(v []float64) float64` | NormInf returns the L-infinity (maximum absolute value) norm of v. |
| `func Normalize(v []float64) []float64` | Normalize returns a new unit-length copy of v under the Euclidean norm. |
| `func NormalizeL1(v []float64) []float64` | NormalizeL1 returns a new copy of v scaled so that its L1 norm is one. |
| `func NormalizedAlgebraicConnectivity(g *Graph) (float64, error)` | NormalizedAlgebraicConnectivity returns the second-smallest eigenvalue of the normalized Laplacian, the spectral quantity appearing in the Cheeger… |
| `func NormalizedCut(g *Graph, part []int) float64` | NormalizedCut returns the normalized-cut value of a bipartition (part entries 0 or 1): cut(A,B)/vol(A) + cut(A,B)/vol(B). |
| `func NormalizedLaplacianSpectrum(g *Graph) ([]float64, error)` | NormalizedLaplacianSpectrum returns the eigenvalues of the symmetric normalized Laplacian in ascending order. |
| `func NumberOfComponentsSpectral(g *Graph, tol float64) (int, error)` | NumberOfComponentsSpectral returns the number of connected components inferred from the spectrum: the multiplicity of the zero eigenvalue of the… |
| `func NumberOfSpanningTrees(g *Graph) (int, error)` | NumberOfSpanningTrees returns `SpanningTreeCount` rounded to the nearest integer, appropriate for unweighted graphs. |
| `func Ones(n int) []float64` | Ones returns a new vector of length n whose entries are all one. |
| `func PageRank(g *Graph, damping float64, maxIter int, tol float64) ([]float64, error)` | PageRank returns the PageRank vector of the graph by power iteration on the Google matrix with the given damping factor (typically 0.85). |
| `func PersonalizedPageRank(g *Graph, damping float64, restart []float64, maxIter int, tol float64) ([]float64, error)` | PersonalizedPageRank returns PageRank with a custom teleportation distribution: at each step the random surfer restarts according to restart (which… |
| `func PowerIteration(a *Matrix, maxIter int, tol float64) (float64, []float64, error)` | PowerIteration estimates the dominant eigenpair (largest magnitude eigenvalue) of a real symmetric matrix. |
| `func Rank(m *Matrix, tol float64) int` | Rank returns the numerical rank of a matrix: the number of pivots found during Gaussian elimination whose magnitude exceeds tol. |
| `func RayleighQuotient(a *Matrix, x []float64) float64` | RayleighQuotient returns xᵀAx / xᵀx, an estimate of an eigenvalue associated with the direction x. |
| `func SignlessLaplacianSpectrum(g *Graph) ([]float64, error)` | SignlessLaplacianSpectrum returns the eigenvalues of the signless Laplacian Q = D + A in ascending order. |
| `func SmallestEigenpair(a *Matrix, maxIter int, tol float64) (float64, []float64, error)` | SmallestEigenpair returns the algebraically smallest eigenvalue of a real symmetric matrix together with a unit eigenvector. |
| `func SolveLinear(a *Matrix, b []float64) ([]float64, error)` | SolveLinear solves the system A·x = b for a square matrix A, using Gaussian elimination with partial pivoting. |
| `func SortedFloats(v []float64) []float64` | SortedFloats returns a sorted (ascending) copy of v. |
| `func SpanningTreeCount(g *Graph) (float64, error)` | SpanningTreeCount returns the (weighted) number of spanning trees via the Matrix-Tree theorem: the determinant of the Laplacian with its last row and… |
| `func SpectralBisection(g *Graph) ([]int, error)` | SpectralBisection partitions the vertices into two blocks using the sign of the Fiedler vector: vertices with a negative Fiedler coordinate are… |
| `func SpectralBisectionMedian(g *Graph) ([]int, error)` | SpectralBisectionMedian partitions the vertices into two balanced blocks by thresholding the Fiedler vector at its median: the vertices whose Fiedler… |
| `func SpectralClustering(g *Graph, k int) ([]int, error)` | SpectralClustering partitions the vertices into k clusters by embedding each vertex into the space spanned by the k eigenvectors of the normalized… |
| `func SpectralGap(g *Graph) (float64, error)` | SpectralGap returns the adjacency spectral gap: the difference between the largest and second-largest adjacency eigenvalues. |
| `func SpectralOrdering(g *Graph) ([]int, error)` | SpectralOrdering returns a permutation of the vertices sorted by their Fiedler vector value in ascending order. |
| `func SpectralRadius(m *Matrix) (float64, error)` | SpectralRadius returns the largest absolute eigenvalue of a symmetric matrix. |
| `func SpectralRadiusAdjacency(g *Graph) (float64, error)` | SpectralRadiusAdjacency returns the largest absolute adjacency eigenvalue, also known as the spectral radius of the graph. |
| `func VecAXPY(a float64, x, y []float64) []float64` | VecAXPY returns a*x + y (a scalar times x, plus y). |
| `func VecAdd(a, b []float64) []float64` | VecAdd returns the element-wise sum a+b. |
| `func VecApproxEqual(a, b []float64, tol float64) bool` | VecApproxEqual reports whether a and b agree to within absolute tolerance tol in every coordinate. |
| `func VecClone(v []float64) []float64` | VecClone returns a copy of v. |
| `func VecDistance(a, b []float64) float64` | VecDistance returns the Euclidean distance between a and b, or +Inf on a length mismatch. |
| `func VecEqual(a, b []float64) bool` | VecEqual reports whether a and b are element-wise identical. |
| `func VecMax(v []float64) float64` | VecMax returns the maximum entry of v, or 0 for an empty slice. |
| `func VecMaxInt(v []int) int` | VecMaxInt returns the maximum entry of an int slice, or -1 for an empty slice. |
| `func VecMean(v []float64) float64` | VecMean returns the arithmetic mean of the entries of v, or 0 for an empty slice. |
| `func VecMin(v []float64) float64` | VecMin returns the minimum entry of v, or 0 for an empty slice. |
| `func VecScale(v []float64, s float64) []float64` | VecScale returns a new vector equal to s*v. |
| `func VecSub(a, b []float64) []float64` | VecSub returns the element-wise difference a-b. |
| `func VecSum(v []float64) float64` | VecSum returns the sum of the entries of v. |
| `func Volume(g *Graph, s map[int]bool) float64` | Volume returns the sum of the weighted degrees of the vertices in S. |
| `func Zeros(n int) []float64` | Zeros returns a new zero vector of length n. |

### Types

| Type | What it is |
| --- | --- |
| `EigenResult` | EigenResult holds the result of a symmetric eigendecomposition: the eigenvalues and their corresponding eigenvectors, stored as the columns of a… |
| `Graph` | Graph is a finite, undirected, optionally weighted graph on the vertex set {0, 1, …, n-1}. |
| `Matrix` | Matrix is a dense matrix of float64 values stored in row-major order. |

<details>
<summary><code>EigenResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func EigenSymmetric(m *Matrix) (*EigenResult, error)` | EigenSymmetric computes all eigenvalues and eigenvectors of a real symmetric matrix using the cyclic Jacobi method. |
| `func (e *EigenResult) Largest() (float64, []float64)` | Largest returns the largest eigenvalue and its eigenvector. |
| `func (e *EigenResult) Len() int` | Len returns the number of eigenpairs. |
| `func (e *EigenResult) Pair(i int) (float64, []float64)` | Pair returns the i-th eigenvalue together with its eigenvector. |
| `func (e *EigenResult) Smallest() (float64, []float64)` | Smallest returns the smallest eigenvalue and its eigenvector. |
| `func (e *EigenResult) SortAscending()` | SortAscending reorders the eigenpairs so the eigenvalues increase. |
| `func (e *EigenResult) SortDescending()` | SortDescending reorders the eigenpairs so the eigenvalues decrease. |
| `func (e *EigenResult) Vector(i int) []float64` | Vector returns the i-th eigenvector (column i of Vectors). |

</details>

<details>
<summary><code>Graph</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CompleteBipartiteGraph(m, n int) *Graph` | CompleteBipartiteGraph returns K_{m,n}: vertices 0..m-1 form one side, vertices m..m+n-1 form the other, and every cross pair is joined. |
| `func CompleteGraph(n int) *Graph` | CompleteGraph returns the complete graph K_n in which every pair of distinct vertices is joined by an edge. |
| `func CycleGraph(n int) *Graph` | CycleGraph returns the cycle C_n on n vertices (n >= 3 for a genuine cycle; for n < 3 it degenerates to a path or single edge). |
| `func EmptyGraph(n int) *Graph` | EmptyGraph returns the edgeless graph on n vertices. |
| `func GridGraph(r, c int) *Graph` | GridGraph returns the r-by-c grid (lattice) graph, the Cartesian product of P_r and P_c, with vertex (i, j) indexed as i*c + j. |
| `func NewGraph(n int) *Graph` | NewGraph returns an edgeless graph on n vertices. |
| `func NewGraphFromAdjacency(a *Matrix) (*Graph, error)` | NewGraphFromAdjacency builds a graph from a square symmetric adjacency matrix. |
| `func PathGraph(n int) *Graph` | PathGraph returns the path P_n: vertices 0..n-1 joined in a line 0-1-2-…-(n-1). |
| `func StarGraph(n int) *Graph` | StarGraph returns the star K_{1,n-1}: a central vertex 0 joined to each of the remaining n-1 vertices. |
| `func WheelGraph(n int) *Graph` | WheelGraph returns the wheel W_n: a cycle on vertices 1..n-1 plus a hub vertex 0 joined to every cycle vertex. |
| `func (g *Graph) AddEdge(u, v int)` | AddEdge adds (or, for a weighted graph, sets to weight 1) an undirected edge between u and v. |
| `func (g *Graph) Adjacency() *Matrix` | Adjacency returns a copy of the adjacency matrix A, where A`i``j` is the weight of the edge between i and j. |
| `func (g *Graph) Complement() *Graph` | Complement returns the complement of an unweighted simple graph: vertices are adjacent in the result exactly when they are non-adjacent (and… |
| `func (g *Graph) ConnectedComponents() []int` | ConnectedComponents returns a labelling of the vertices into connected components, using a breadth-first flood fill over the (unweighted) edge set. |
| `func (g *Graph) Degree(v int) int` | Degree returns the (integer) number of neighbours of v in an unweighted graph, counting a self-loop once. |
| `func (g *Graph) DegreeMatrix() *Matrix` | DegreeMatrix returns the diagonal degree matrix D whose (v, v) entry is the weighted degree of v. |
| `func (g *Graph) Degrees() []float64` | Degrees returns the weighted degree of every vertex. |
| `func (g *Graph) Density() float64` | Density returns the edge density of an unweighted view: the number of edges (ignoring self-loops) divided by the number of possible edges n(n-1)/2. |
| `func (g *Graph) HasEdge(u, v int) bool` | HasEdge reports whether there is a (nonzero-weight) edge between u and v. |
| `func (g *Graph) IncidenceMatrix() *Matrix` | IncidenceMatrix returns the (oriented) vertex-edge incidence matrix B of an unweighted view of the graph. |
| `func (g *Graph) IsBipartite() (bool, []int)` | IsBipartite reports whether the graph is bipartite, and if so returns a two-colouring (side 0 or 1) of its vertices. |
| `func (g *Graph) IsConnected() bool` | IsConnected reports whether the graph has exactly one connected component. |
| `func (g *Graph) IsRegular() (bool, float64)` | IsRegular reports whether every vertex has the same weighted degree (to within 1e-9), and returns that common degree. |
| `func (g *Graph) Laplacian() *Matrix` | Laplacian returns the combinatorial Laplacian L = D - A. |
| `func (g *Graph) MaxDegree() float64` | MaxDegree returns the largest weighted degree, or 0 for an empty graph. |
| `func (g *Graph) MinDegree() float64` | MinDegree returns the smallest weighted degree, or 0 for an empty graph. |
| `func (g *Graph) Neighbors(v int) []int` | Neighbors returns the vertices adjacent to v, in ascending order. |
| `func (g *Graph) NormalizedLaplacian() *Matrix` | NormalizedLaplacian returns the symmetric normalized Laplacian Lsym = I - D^{-1/2} A D^{-1/2}. |
| `func (g *Graph) NumConnectedComponents() int` | NumConnectedComponents returns the number of connected components. |
| `func (g *Graph) Order() int` | Order returns the number of vertices. |
| `func (g *Graph) RandomWalkLaplacian() *Matrix` | RandomWalkLaplacian returns the random-walk Laplacian Lrw = I - D^{-1} A. Isolated vertices contribute a unit diagonal entry. |
| `func (g *Graph) RemoveEdge(u, v int)` | RemoveEdge deletes the edge between u and v. |
| `func (g *Graph) SetWeight(u, v int, w float64)` | SetWeight sets the weight of the undirected edge between u and v (and, by symmetry, between v and u). |
| `func (g *Graph) SignlessLaplacian() *Matrix` | SignlessLaplacian returns the signless Laplacian Q = D + A. |
| `func (g *Graph) Size() int` | Size returns the number of edges. |
| `func (g *Graph) TotalWeight() float64` | TotalWeight returns the sum of all edge weights, counting each undirected edge once. |
| `func (g *Graph) TransitionMatrix() *Matrix` | TransitionMatrix returns the random-walk transition matrix P = D^{-1} A, whose row i is the probability distribution of a single step from vertex i. |
| `func (g *Graph) Weight(u, v int) float64` | Weight returns the weight of the edge between u and v (0 if there is none). |
| `func (g *Graph) WeightedDegree(v int) float64` | WeightedDegree returns the weighted degree of vertex v: the sum of the weights of its incident edges. |

</details>

<details>
<summary><code>Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AllPairsShortestPaths(g *Graph) *Matrix` | AllPairsShortestPaths returns the matrix of weighted shortest-path distances between every ordered pair of vertices, computed by running Dijkstra… |
| `func DiagMatrix(d []float64) *Matrix` | DiagMatrix returns a square diagonal matrix whose diagonal is d. |
| `func IdentityMatrix(n int) *Matrix` | IdentityMatrix returns the n-by-n identity matrix. |
| `func Inverse(a *Matrix) (*Matrix, error)` | Inverse returns the inverse of a square matrix using Gauss-Jordan elimination with partial pivoting. |
| `func LaplacianPseudoinverse(g *Graph, tol float64) (*Matrix, error)` | LaplacianPseudoinverse returns the Moore-Penrose pseudoinverse L⁺ of the combinatorial Laplacian, computed from its eigendecomposition by inverting… |
| `func NewMatrix(rows, cols int) *Matrix` | NewMatrix returns a new rows-by-cols matrix with every entry zero. |
| `func NewMatrixFromRows(rows [][]float64) (*Matrix, error)` | NewMatrixFromRows builds a matrix from a slice of equal-length rows. |
| `func OuterProduct(a, b []float64) *Matrix` | OuterProduct returns the rank-one matrix a·bᵀ. |
| `func (m *Matrix) Add(i, j int, v float64)` | Add increments the entry in row i, column j by v. |
| `func (m *Matrix) Apply(f func(float64) float64) *Matrix` | Apply returns a new matrix whose entries are f applied element-wise. |
| `func (m *Matrix) ApproxEqual(o *Matrix, tol float64) bool` | ApproxEqual reports whether m and o have the same shape and agree entry-wise to within absolute tolerance tol. |
| `func (m *Matrix) At(i, j int) float64` | At returns the entry in row i, column j. |
| `func (m *Matrix) Clone() *Matrix` | Clone returns an independent copy of the matrix. |
| `func (m *Matrix) Col(j int) []float64` | Col returns a copy of column j. |
| `func (m *Matrix) ColSums() []float64` | ColSums returns the sum of each column. |
| `func (m *Matrix) Cols() int` | Cols returns the number of columns. |
| `func (m *Matrix) Diagonal() []float64` | Diagonal returns a copy of the main diagonal. |
| `func (m *Matrix) Equals(o *Matrix) bool` | Equals reports whether m and o have the same shape and identical entries. |
| `func (m *Matrix) FrobeniusNorm() float64` | FrobeniusNorm returns the square root of the sum of squared entries. |
| `func (m *Matrix) IsSquare() bool` | IsSquare reports whether the matrix has equal numbers of rows and columns. |
| `func (m *Matrix) IsSymmetric(tol float64) bool` | IsSymmetric reports whether the matrix is square and symmetric to within absolute tolerance tol. |
| `func (m *Matrix) MaxAbs() float64` | MaxAbs returns the largest absolute entry. |
| `func (m *Matrix) Minus(o *Matrix) (*Matrix, error)` | Minus returns the difference m-o. |
| `func (m *Matrix) Mul(o *Matrix) (*Matrix, error)` | Mul returns the matrix product m*o. |
| `func (m *Matrix) MulVec(v []float64) ([]float64, error)` | MulVec returns the matrix-vector product m*v. |
| `func (m *Matrix) Plus(o *Matrix) (*Matrix, error)` | Plus returns the sum m+o. |
| `func (m *Matrix) Pow(p int) (*Matrix, error)` | Pow returns m raised to the non-negative integer power p by repeated squaring. |
| `func (m *Matrix) Row(i int) []float64` | Row returns a copy of row i. |
| `func (m *Matrix) RowSums() []float64` | RowSums returns the sum of each row. |
| `func (m *Matrix) Rows() int` | Rows returns the number of rows. |
| `func (m *Matrix) Scale(s float64) *Matrix` | Scale returns a new matrix equal to s*m. |
| `func (m *Matrix) Set(i, j int, v float64)` | Set assigns v to the entry in row i, column j. |
| `func (m *Matrix) SetRow(i int, r []float64) error` | SetRow overwrites row i with the values in r. |
| `func (m *Matrix) String() string` | String renders the matrix with each entry formatted using %g. |
| `func (m *Matrix) SubMatrix(dropRows, dropCols map[int]bool) *Matrix` | SubMatrix returns the matrix obtained by deleting the given set of row indices and column indices. |
| `func (m *Matrix) Symmetrize() *Matrix` | Symmetrize returns (m + mᵀ)/2, the symmetric part of a square matrix. |
| `func (m *Matrix) ToRows() [][]float64` | ToRows returns the matrix as a fresh slice of row slices. |
| `func (m *Matrix) Trace() float64` | Trace returns the sum of the diagonal entries. |
| `func (m *Matrix) Transpose() *Matrix` | Transpose returns the transpose of the matrix. |

</details>

### Variables

`ErrDimensionMismatch`, `ErrNotSquare`, `ErrNotSymmetric`, `ErrSingular`, `ErrEmpty`, `ErrOutOfRange`, `ErrInvalidArgument`, `ErrNotConnected`, `ErrNoConvergence`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/graphspectral).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
