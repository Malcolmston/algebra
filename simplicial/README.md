# simplicial — implements abstract simplicial complexes and the core machinery of computational

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/simplicial.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/simplicial)

Package simplicial implements abstract simplicial complexes and the core
machinery of computational topology over the Go standard library alone.

The fundamental combinatorial object is a `Simplex`: a finite, non-empty set
of integer vertices, stored in ascending order. A k-simplex has k+1 vertices
and dimension k. Simplices are glued together into an abstract `Complex`, a
downward-closed family of simplices — whenever a simplex belongs to the
complex so does every one of its faces. From a complex the package builds the
chain complex of boundary maps

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/simplicial"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func ChebyshevDistance(a, b []float64) float64` | ChebyshevDistance returns the L∞ (maximum coordinate) distance. |
| `func Circumsphere(points [][]float64) (center []float64, radius float64, ok bool)` | Circumsphere returns the centre and radius of the unique sphere passing through all of the given (affinely independent) points, with the centre lying… |
| `func CompareSimplices(a, b Simplex) int` | CompareSimplices provides a total order on simplices: first by dimension, then lexicographically by vertex list. |
| `func EnclosingRadius(points [][]float64) float64` | EnclosingRadius returns the radius of the minimal enclosing ball of the points, a convenience wrapper around `MinimalEnclosingBall`. |
| `func EuclideanDistance(a, b []float64) float64` | EuclideanDistance returns the L2 distance between two equal-length points. |
| `func ManhattanDistance(a, b []float64) float64` | ManhattanDistance returns the L1 (taxicab) distance. |
| `func MinimalEnclosingBall(points [][]float64) (center []float64, radius float64)` | MinimalEnclosingBall returns the centre and radius of the smallest ball containing every one of the given points, computed exactly with Welzl's… |
| `func MinkowskiDistance(a, b []float64, p float64) float64` | MinkowskiDistance returns the L_p distance for p ≥ 1. |
| `func SortSimplices(ss []Simplex)` | SortSimplices sorts a slice of simplices in place by `CompareSimplices`. |
| `func SquaredEuclideanDistance(a, b []float64) float64` | SquaredEuclideanDistance returns the squared L2 distance, avoiding the square root. |

### Types

| Type | What it is |
| --- | --- |
| `Bar` | Bar is a half-open interval [Birth, Death) in a persistence barcode. |
| `BoundaryTerm` | BoundaryTerm is one oriented summand of a simplex's boundary: a face together with the sign (+1 or −1) with which it appears. |
| `Complex` | Complex is an abstract simplicial complex: a finite, downward-closed family of `Simplex` values. |
| `FilteredSimplex` | FilteredSimplex is a simplex together with the scalar value at which it enters a filtration. |
| `Filtration` | Filtration is a simplicial complex whose simplices are equipped with non-decreasing entrance values: every face of a simplex must enter no later than… |
| `GF2Matrix` | GF2Matrix is a dense matrix over the two-element field GF(2) = {0,1}, with arithmetic performed modulo 2. |
| `HomologyGroup` | HomologyGroup describes a finitely generated abelian group H_k = Z^r ⊕ (Z/t_1 ⊕ … ⊕ Z/t_m) by its free rank r and its torsion coefficients… |
| `IntMatrix` | IntMatrix is a dense matrix over the integers Z, with entries stored exactly as *big.Int. |
| `Metric` | Metric is a distance function between two points of equal length. |
| `Persistence` | Persistence holds the result of a persistent-homology computation: the list of birth/death pairs across all dimensions. |
| `PersistencePair` | PersistencePair records the birth and death of a single homology class of a given dimension during a filtration. |
| `PointCloud` | PointCloud is a finite set of points in Euclidean space, all of the same ambient dimension. |
| `RatMatrix` | RatMatrix is a dense matrix over the field of rationals Q, with entries stored exactly as *big.Rat. |
| `Simplex` | Simplex is a finite set of distinct integer vertices, stored in ascending order. |
| `SmithResult` | SmithResult holds the Smith normal form decomposition U·A·V = D of an integer matrix A, where U and V are unimodular (invertible over Z) and D is… |

<details>
<summary><code>Bar</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (b Bar) IsInfinite() bool` | IsInfinite reports whether the bar extends to +Inf. |
| `func (b Bar) Length() float64` | Length returns the length of the bar (+Inf for an essential bar). |

</details>

<details>
<summary><code>Complex</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoundarySphere(n int) *Complex` | BoundarySphere returns the boundary of the standard (n+1)-simplex: the n-skeleton of an (n+1)-simplex, a triangulation of the n-sphere S^n on n+2… |
| `func Cech(pc *PointCloud, r float64, maxDim int) *Complex` | Cech builds the Čech complex of a point cloud at radius r up to dimension maxDim. |
| `func CompleteGraphComplex(n int) *Complex` | CompleteGraphComplex returns the 1-skeleton of the full simplex on n vertices: the complete graph K_n with all n(n−1)/2 edges but no higher faces. |
| `func Cone(c *Complex, apex int) *Complex` | Cone returns the cone on c with a new apex vertex: every simplex σ of c is joined to apex, together with the apex itself. |
| `func CycleComplex(n int) *Complex` | CycleComplex returns the cycle graph C_n: n vertices 0..n−1 joined in a ring by n edges. |
| `func DiscreteComplex(n int) *Complex` | DiscreteComplex returns n isolated vertices 0..n−1 with no edges; its degree-0 homology has rank n and all higher homology vanishes. |
| `func DisjointUnion(c, d *Complex, offset int) *Complex` | DisjointUnion returns the disjoint union of c and d. |
| `func FromSimplices(ss ...Simplex) *Complex` | FromSimplices builds the smallest simplicial complex containing all of the given simplices, i.e. |
| `func NeighborhoodGraph(pc *PointCloud, epsilon float64) *Complex` | NeighborhoodGraph returns the 1-skeleton (graph) of the Vietoris–Rips complex of a point cloud at scale epsilon: vertices for every point and an… |
| `func NewComplex() *Complex` | NewComplex returns a new, empty simplicial complex. |
| `func PathComplex(n int) *Complex` | PathComplex returns the path graph on n vertices 0..n−1 with n−1 edges. |
| `func PointComplex() *Complex` | PointComplex returns the complex consisting of a single vertex {0}: a contractible space with the homology of a point. |
| `func ProjectivePlaneComplex() *Complex` | ProjectivePlaneComplex returns the minimal 6-vertex triangulation of the real projective plane RP^2 (the hemi-icosahedron) on vertices 0..5, with 15… |
| `func SphereComplex(n int) *Complex` | SphereComplex returns a triangulation of the n-sphere S^n as the boundary of the standard (n+1)-simplex. |
| `func StandardSimplex(n int) *Complex` | StandardSimplex returns the full n-simplex on vertices 0..n together with all of its faces. |
| `func Suspension(c *Complex, north, south int) *Complex` | Suspension returns the suspension of c: two cones on c glued along c, using the two new apex vertices north and south. |
| `func TorusComplex() *Complex` | TorusComplex returns Möbius's minimal 7-vertex triangulation of the torus T^2, generated on vertices 0..6 (indices taken mod 7) by the triangles… |
| `func VietorisRips(pc *PointCloud, epsilon float64, maxDim int) *Complex` | VietorisRips builds the Vietoris–Rips complex of a point cloud at scale epsilon under the Euclidean metric, up to dimension maxDim. |
| `func VietorisRipsFromDistances(dist [][]float64, epsilon float64, maxDim int) *Complex` | VietorisRipsFromDistances builds the Vietoris–Rips complex at scale epsilon directly from a symmetric distance matrix, up to dimension maxDim. |
| `func VietorisRipsMetric(pc *PointCloud, epsilon float64, maxDim int, metric Metric) *Complex` | VietorisRipsMetric builds the Vietoris–Rips complex of a point cloud at scale epsilon under an arbitrary metric, up to dimension maxDim. |
| `func (c *Complex) AddEdge(a, b int)` | AddEdge inserts the 1-simplex {a,b} and its two endpoints. |
| `func (c *Complex) AddSimplex(s Simplex)` | AddSimplex inserts the simplex and every one of its faces into the complex, preserving downward closure. |
| `func (c *Complex) AddSimplices(ss ...Simplex)` | AddSimplices inserts several simplices and their faces. |
| `func (c *Complex) AddTriangle(a, b, c2 int)` | AddTriangle inserts the 2-simplex {a,b,c} and all of its faces. |
| `func (c *Complex) AddVertex(v int)` | AddVertex inserts the 0-simplex {v}. |
| `func (c *Complex) BettiNumber(k int) int` | BettiNumber returns the k-th Betti number b_k of the complex, the rank of the k-th homology group with rational coefficients. |
| `func (c *Complex) BettiNumberGF2(k int) int` | BettiNumberGF2 returns the k-th Betti number computed with GF(2) coefficients. |
| `func (c *Complex) BettiNumbers() []int` | BettiNumbers returns the rational Betti numbers b_0,…,b_d of the complex, where d is its dimension. |
| `func (c *Complex) BettiNumbersGF2() []int` | BettiNumbersGF2 returns the GF(2) Betti numbers b_0,…,b_d of the complex. |
| `func (c *Complex) BoundaryMatrixGF2(k int) *GF2Matrix` | BoundaryMatrixGF2 returns the matrix of the boundary map ∂_k : C_k → C_{k-1} over GF(2). |
| `func (c *Complex) BoundaryMatrixQ(k int) *RatMatrix` | BoundaryMatrixQ returns the boundary map ∂_k over the rationals, identical in content to `Complex.BoundaryMatrixZ` but with rational entries. |
| `func (c *Complex) BoundaryMatrixZ(k int) *IntMatrix` | BoundaryMatrixZ returns the matrix of the boundary map ∂_k : C_k → C_{k-1} over the integers, with the standard orientation signs. |
| `func (c *Complex) BoundaryRankGF2(k int) int` | BoundaryRankGF2 returns the rank over GF(2) of the boundary map ∂_k. |
| `func (c *Complex) BoundaryRankQ(k int) int` | BoundaryRankQ returns the rank over Q of the boundary map ∂_k. |
| `func (c *Complex) ChainRank(k int) int` | ChainRank returns the dimension of the chain group C_k, i.e. |
| `func (c *Complex) ClosedStar(s Simplex) *Complex` | ClosedStar returns the closure of the star of s as a subcomplex: the smallest subcomplex containing every coface of s. |
| `func (c *Complex) Cofaces(s Simplex) []Simplex` | Cofaces returns every simplex of the complex that has s as a face, s itself included when present. |
| `func (c *Complex) ComponentOf(v int) []int` | ComponentOf returns the sorted vertex list of the connected component containing vertex v, or nil if v is not a vertex of the complex. |
| `func (c *Complex) ConnectedComponents() [][]int` | ConnectedComponents partitions the vertices into connected components and returns, for each component, its sorted vertex list. |
| `func (c *Complex) Copy() *Complex` | Copy returns an independent copy of the complex. |
| `func (c *Complex) Dimension() int` | Dimension returns the largest dimension of any simplex in the complex, or −1 for the empty complex. |
| `func (c *Complex) Equal(d *Complex) bool` | Equal reports whether c and d contain exactly the same simplices. |
| `func (c *Complex) EulerCharacteristic() int` | EulerCharacteristic returns the alternating sum Σ_k (−1)^k f_k of the f-vector, the Euler characteristic of the complex. |
| `func (c *Complex) EulerCharacteristicFromBetti() int` | EulerCharacteristicFromBetti returns the Euler characteristic computed as the alternating sum of the Betti numbers, Σ_k (−1)^k b_k. |
| `func (c *Complex) FVector() []int` | FVector returns the f-vector of the complex: entry k is the number of k-simplices, for k from 0 up to the dimension of the complex. |
| `func (c *Complex) Has(s Simplex) bool` | Has reports whether the given simplex belongs to the complex. |
| `func (c *Complex) HasVertex(v int) bool` | HasVertex reports whether the 0-simplex {v} belongs to the complex. |
| `func (c *Complex) Homology() []HomologyGroup` | Homology returns the integral homology groups H_0,…,H_d of the complex. |
| `func (c *Complex) HomologyRank(k int) int` | HomologyRank is an alias for `Complex.BettiNumber`: the free rank of H_k. |
| `func (c *Complex) HomologyZ(k int) HomologyGroup` | HomologyZ returns the k-th integral homology group H_k(c; Z) as a `HomologyGroup`, combining the Betti number (free rank) with the torsion… |
| `func (c *Complex) IsConnected() bool` | IsConnected reports whether the complex has at most one connected component. |
| `func (c *Complex) IsCycleGF2(k int, chain []int) bool` | IsCycleGF2 reports whether the given GF(2) chain — a 0/1 coefficient vector indexed like `Complex.SimplicesOfDim`(k) — is a cycle, i.e. |
| `func (c *Complex) IsValid() bool` | IsValid reports whether the complex is genuinely downward-closed. |
| `func (c *Complex) Link(s Simplex) *Complex` | Link returns the link of s: the subcomplex of faces τ that are disjoint from s yet whose join with s lies in the complex. |
| `func (c *Complex) MaximalSimplices() []Simplex` | MaximalSimplices returns the maximal simplices of the complex — those that are not a proper face of any other simplex. |
| `func (c *Complex) NumConnectedComponents() int` | NumConnectedComponents returns the number of connected components of the 1-skeleton of the complex — equivalently the rank b_0 of degree-0 homology. |
| `func (c *Complex) NumSimplices() int` | NumSimplices returns the total number of simplices of every dimension. |
| `func (c *Complex) NumSimplicesOfDim(k int) int` | NumSimplicesOfDim returns the number of k-simplices in the complex. |
| `func (c *Complex) NumVertices() int` | NumVertices returns the number of 0-simplices (vertices) of the complex. |
| `func (c *Complex) RemoveSimplex(s Simplex)` | RemoveSimplex deletes s and every coface of s from the complex, keeping it downward-closed. |
| `func (c *Complex) Simplices() []Simplex` | Simplices returns all simplices of the complex sorted by `CompareSimplices`. |
| `func (c *Complex) SimplicesOfDim(k int) []Simplex` | SimplicesOfDim returns all k-simplices of the complex sorted by `CompareSimplices`. |
| `func (c *Complex) Skeleton(k int) *Complex` | Skeleton returns the k-skeleton of the complex: a new complex containing exactly the simplices of dimension at most k. |
| `func (c *Complex) Star(s Simplex) []Simplex` | Star returns the open star of s: all simplices of the complex having s as a face. |
| `func (c *Complex) String() string` | String returns a compact multi-line description listing the simplices by dimension. |
| `func (c *Complex) TorsionCoefficients(k int) []*big.Int` | TorsionCoefficients returns the torsion coefficients of the k-th integral homology group H_k(c; Z): the elementary divisors greater than one of the… |
| `func (c *Complex) Union(d *Complex) *Complex` | Union returns the simplicial complex whose simplices are those of c or d. |
| `func (c *Complex) Vertices() []int` | Vertices returns the sorted list of vertices appearing in the complex. |

</details>

<details>
<summary><code>Filtration</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewFiltration() *Filtration` | NewFiltration returns an empty filtration. |
| `func RipsFiltration(pc *PointCloud, maxEps float64, maxDim int) *Filtration` | RipsFiltration builds the Vietoris–Rips filtration of a point cloud up to the scale maxEps and dimension maxDim. |
| `func (f *Filtration) Add(s Simplex, value float64)` | Add inserts a single simplex at the given value. |
| `func (f *Filtration) AddClosure(s Simplex, value float64)` | AddClosure inserts s and every one of its faces at the given value (faces take the smaller of their existing and the new value), preserving validity. |
| `func (f *Filtration) ComplexAt(t float64) *Complex` | ComplexAt returns the subcomplex of all simplices with value at most t. |
| `func (f *Filtration) Len() int` | Len returns the number of simplices in the filtration. |
| `func (f *Filtration) Simplices() []FilteredSimplex` | Simplices returns the filtered simplices in filtration order: ascending by value, then by dimension, then by `CompareSimplices`. |
| `func (f *Filtration) Value(s Simplex) (float64, bool)` | Value returns the filtration value of the simplex and whether it is present. |

</details>

<details>
<summary><code>GF2Matrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func GF2Identity(n int) *GF2Matrix` | GF2Identity returns the n×n identity matrix over GF(2). |
| `func NewGF2Matrix(rows, cols int) *GF2Matrix` | NewGF2Matrix returns a rows×cols zero matrix over GF(2). |
| `func (m *GF2Matrix) Add(n *GF2Matrix) (*GF2Matrix, error)` | Add returns the entrywise sum (XOR) m+n over GF(2). |
| `func (m *GF2Matrix) AddRow(dst, src int)` | AddRow adds (XORs) row src into row dst in place. |
| `func (m *GF2Matrix) At(i, j int) int` | At returns the entry in row i, column j as 0 or 1. |
| `func (m *GF2Matrix) Clone() *GF2Matrix` | Clone returns an independent copy of the matrix. |
| `func (m *GF2Matrix) Col(j int) []int` | Col returns a copy of column j as a slice of 0/1 ints. |
| `func (m *GF2Matrix) Cols() int` | Cols returns the number of columns. |
| `func (m *GF2Matrix) Equal(n *GF2Matrix) bool` | Equal reports whether m and n have the same shape and entries. |
| `func (m *GF2Matrix) IsZero() bool` | IsZero reports whether every entry is zero. |
| `func (m *GF2Matrix) KernelBasis() [][]int` | KernelBasis returns a basis of the null space {x : m·x = 0} over GF(2). |
| `func (m *GF2Matrix) Mul(n *GF2Matrix) (*GF2Matrix, error)` | Mul returns the matrix product m·n over GF(2). |
| `func (m *GF2Matrix) Nullity() int` | Nullity returns the dimension of the kernel of the matrix acting on column vectors, namely cols − rank. |
| `func (m *GF2Matrix) Rank() int` | Rank returns the rank of the matrix over GF(2). |
| `func (m *GF2Matrix) Row(i int) []int` | Row returns a copy of row i as a slice of 0/1 ints. |
| `func (m *GF2Matrix) RowEchelon() *GF2Matrix` | RowEchelon returns the reduced row-echelon form of the matrix over GF(2). |
| `func (m *GF2Matrix) Rows() int` | Rows returns the number of rows. |
| `func (m *GF2Matrix) Set(i, j, v int)` | Set stores v (reduced modulo 2) in row i, column j. |
| `func (m *GF2Matrix) String() string` | String returns a multi-line 0/1 rendering of the matrix. |
| `func (m *GF2Matrix) SwapRows(i, j int)` | SwapRows exchanges rows i and j in place. |
| `func (m *GF2Matrix) Transpose() *GF2Matrix` | Transpose returns the transpose of the matrix. |

</details>

<details>
<summary><code>HomologyGroup</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (h HomologyGroup) Betti() int` | Betti returns the free rank of the group. |
| `func (h HomologyGroup) IsFree() bool` | IsFree reports whether the group is torsion-free. |
| `func (h HomologyGroup) IsTrivial() bool` | IsTrivial reports whether the group is the zero group. |
| `func (h HomologyGroup) String() string` | String renders the group in the usual notation, e.g. |

</details>

<details>
<summary><code>IntMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func IntIdentity(n int) *IntMatrix` | IntIdentity returns the n×n identity matrix over Z. |
| `func IntMatrixFromInt64(r, c int, vals []int64) *IntMatrix` | IntMatrixFromInt64 builds an r×c integer matrix from a row-major slice of int64 values. |
| `func NewIntMatrix(rows, cols int) *IntMatrix` | NewIntMatrix returns a rows×cols zero matrix over Z. |
| `func (m *IntMatrix) Add(n *IntMatrix) (*IntMatrix, error)` | Add returns m+n. |
| `func (m *IntMatrix) At(i, j int) *big.Int` | At returns a copy of the entry in row i, column j. |
| `func (m *IntMatrix) Clone() *IntMatrix` | Clone returns an independent copy of the matrix. |
| `func (m *IntMatrix) Cols() int` | Cols returns the number of columns. |
| `func (m *IntMatrix) Determinant() (*big.Int, error)` | Determinant returns the determinant of a square integer matrix. |
| `func (m *IntMatrix) Equal(n *IntMatrix) bool` | Equal reports whether m and n have the same shape and entries. |
| `func (m *IntMatrix) InvariantFactors() []*big.Int` | InvariantFactors returns the non-zero invariant factors of the matrix. |
| `func (m *IntMatrix) IsZero() bool` | IsZero reports whether every entry is zero. |
| `func (m *IntMatrix) Mul(n *IntMatrix) (*IntMatrix, error)` | Mul returns the product m·n. |
| `func (m *IntMatrix) Neg() *IntMatrix` | Neg returns −m. |
| `func (m *IntMatrix) Rank() int` | Rank returns the rank of the integer matrix, computed over Q (equivalently the number of non-zero invariant factors of its Smith normal form). |
| `func (m *IntMatrix) Rows() int` | Rows returns the number of rows. |
| `func (m *IntMatrix) Set(i, j int, v *big.Int)` | Set stores a copy of v in row i, column j. |
| `func (m *IntMatrix) SetInt(i, j int, v int64)` | SetInt stores the integer v in row i, column j. |
| `func (m *IntMatrix) SmithNormalForm() *SmithResult` | SmithNormalForm computes the Smith normal form of the matrix and returns U, D and V with U·A·V = D. The original matrix is not modified. |
| `func (m *IntMatrix) SmithNormalFormDiagonal() []*big.Int` | SmithNormalFormDiagonal returns just the diagonal (invariant factors) of the Smith normal form. |
| `func (m *IntMatrix) String() string` | String returns a multi-line rendering of the matrix. |
| `func (m *IntMatrix) Sub(n *IntMatrix) (*IntMatrix, error)` | Sub returns m−n. |
| `func (m *IntMatrix) ToRat() *RatMatrix` | ToRat returns the matrix reinterpreted over the rationals. |
| `func (m *IntMatrix) Transpose() *IntMatrix` | Transpose returns the transpose of the matrix. |

</details>

<details>
<summary><code>Persistence</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func PersistentHomology(f *Filtration) *Persistence` | PersistentHomology runs the standard matrix-reduction algorithm over GF(2) on the given filtration and returns the resulting persistence pairs. |
| `func (p *Persistence) Barcode(k int) []Bar` | Barcode returns the persistence barcode in dimension k: one `Bar` per feature. |
| `func (p *Persistence) BettiAt(k int, t float64) int` | BettiAt returns the k-th Betti number of the subcomplex at threshold t, counting the features of dimension k that are alive at t (born at or before t… |
| `func (p *Persistence) Diagram(k int) [][2]float64` | Diagram returns the persistence diagram in dimension k as a slice of [birth, death] points; essential classes have death +Inf. |
| `func (p *Persistence) EssentialClasses() []PersistencePair` | EssentialClasses returns the pairs that never die (infinite persistence), whose count in dimension k is the k-th Betti number of the final complex. |
| `func (p *Persistence) NumPairs() int` | NumPairs returns the total number of persistence pairs. |
| `func (p *Persistence) Pairs() []PersistencePair` | Pairs returns all persistence pairs across every dimension. |
| `func (p *Persistence) PairsOfDim(k int) []PersistencePair` | PairsOfDim returns the persistence pairs of homological dimension k. |

</details>

<details>
<summary><code>PersistencePair</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (p PersistencePair) IsEssential() bool` | IsEssential reports whether the class never dies (infinite persistence). |
| `func (p PersistencePair) Persistence() float64` | Persistence returns Death − Birth, the lifetime of the feature (+Inf for essential classes). |

</details>

<details>
<summary><code>PointCloud</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CirclePoints(n int, radius float64) *PointCloud` | CirclePoints returns n points equally spaced on the circle of the given radius centred at the origin in the plane. |
| `func GridPoints(dims ...int) *PointCloud` | GridPoints returns a point cloud of the integer lattice points of an n-dimensional box: for dims = [a,b,c] it produces a·b·c points at all integer… |
| `func NewPointCloud(points ...[]float64) *PointCloud` | NewPointCloud returns a point cloud from the given points. |
| `func (pc *PointCloud) Add(p []float64) int` | Add appends a copy of point p to the cloud and returns its index. |
| `func (pc *PointCloud) Centroid() []float64` | Centroid returns the coordinatewise mean of the points, or nil for an empty cloud. |
| `func (pc *PointCloud) Diameter() float64` | Diameter returns the largest pairwise Euclidean distance in the cloud, or 0 for a cloud with fewer than two points. |
| `func (pc *PointCloud) Dim() int` | Dim returns the ambient dimension (the length of the first point), or 0 for an empty cloud. |
| `func (pc *PointCloud) Distance(i, j int) float64` | Distance returns the Euclidean distance between points i and j. |
| `func (pc *PointCloud) DistanceMatrix() [][]float64` | DistanceMatrix returns the symmetric matrix of pairwise Euclidean distances. |
| `func (pc *PointCloud) DistanceMatrixWith(metric Metric) [][]float64` | DistanceMatrixWith returns the symmetric matrix of pairwise distances under the given metric. |
| `func (pc *PointCloud) DistanceWith(i, j int, metric Metric) float64` | DistanceWith returns the distance between points i and j under the given metric. |
| `func (pc *PointCloud) Len() int` | Len returns the number of points in the cloud. |
| `func (pc *PointCloud) Point(i int) []float64` | Point returns a copy of the i-th point. |

</details>

<details>
<summary><code>RatMatrix</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewRatMatrix(rows, cols int) *RatMatrix` | NewRatMatrix returns a rows×cols zero matrix over Q. |
| `func RatIdentity(n int) *RatMatrix` | RatIdentity returns the n×n identity matrix over Q. |
| `func (m *RatMatrix) Add(n *RatMatrix) (*RatMatrix, error)` | Add returns m+n. |
| `func (m *RatMatrix) At(i, j int) *big.Rat` | At returns a copy of the entry in row i, column j. |
| `func (m *RatMatrix) Clone() *RatMatrix` | Clone returns an independent copy of the matrix. |
| `func (m *RatMatrix) Cols() int` | Cols returns the number of columns. |
| `func (m *RatMatrix) Determinant() (*big.Rat, error)` | Determinant returns the determinant of a square matrix over Q. It returns `ErrNotSquare` for a rectangular matrix. |
| `func (m *RatMatrix) Equal(n *RatMatrix) bool` | Equal reports whether m and n have the same shape and equal entries. |
| `func (m *RatMatrix) IsZero() bool` | IsZero reports whether every entry is zero. |
| `func (m *RatMatrix) KernelBasis() [][]*big.Rat` | KernelBasis returns a basis of the null space {x : m·x = 0} over Q. Each basis vector is a length-cols slice of *big.Rat. |
| `func (m *RatMatrix) Mul(n *RatMatrix) (*RatMatrix, error)` | Mul returns the product m·n. |
| `func (m *RatMatrix) Nullity() int` | Nullity returns cols − rank, the dimension of the kernel over Q. |
| `func (m *RatMatrix) Rank() int` | Rank returns the rank of the matrix over Q. |
| `func (m *RatMatrix) RowEchelon() *RatMatrix` | RowEchelon returns the reduced row-echelon form of the matrix over Q. |
| `func (m *RatMatrix) Rows() int` | Rows returns the number of rows. |
| `func (m *RatMatrix) Scale(s *big.Rat) *RatMatrix` | Scale returns the matrix with every entry multiplied by s. |
| `func (m *RatMatrix) Set(i, j int, v *big.Rat)` | Set stores a copy of v in row i, column j. |
| `func (m *RatMatrix) SetInt(i, j int, v int64)` | SetInt stores the integer v in row i, column j. |
| `func (m *RatMatrix) String() string` | String returns a multi-line rendering of the matrix. |
| `func (m *RatMatrix) Sub(n *RatMatrix) (*RatMatrix, error)` | Sub returns m−n. |
| `func (m *RatMatrix) Transpose() *RatMatrix` | Transpose returns the transpose of the matrix. |

</details>

<details>
<summary><code>Simplex</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Edge(a, b int) Simplex` | Edge returns the 1-simplex on vertices a and b. |
| `func EmptySimplex() Simplex` | EmptySimplex returns the (−1)-dimensional simplex with no vertices. |
| `func NewSimplex(verts ...int) Simplex` | NewSimplex returns the simplex on the given vertices. |
| `func Tetrahedron(a, b, c, d int) Simplex` | Tetrahedron returns the 3-simplex on vertices a, b, c and d. |
| `func Triangle(a, b, c int) Simplex` | Triangle returns the 2-simplex on vertices a, b and c. |
| `func Vertex(v int) Simplex` | Vertex returns the 0-simplex on the single vertex v. |
| `func (s Simplex) Boundary() []BoundaryTerm` | Boundary returns the oriented algebraic boundary ∂s of the simplex, the alternating sum Σ_i (−1)^i [v_0,…,v̂_i,…,v_k] of its codimension-1… |
| `func (s Simplex) Closure() []Simplex` | Closure returns all faces of the simplex, itself included: every non-empty subset of its vertices. |
| `func (s Simplex) Contains(v int) bool` | Contains reports whether v is a vertex of the simplex. |
| `func (s Simplex) ContainsSimplex(t Simplex) bool` | ContainsSimplex reports whether t is a (not necessarily proper) face of s, that is, whether every vertex of t is a vertex of s. |
| `func (s Simplex) Copy() Simplex` | Copy returns an independent copy of the simplex. |
| `func (s Simplex) Dim() int` | Dim returns the dimension of the simplex, one less than its number of vertices. |
| `func (s Simplex) Equal(t Simplex) bool` | Equal reports whether s and t have exactly the same vertex set. |
| `func (s Simplex) Face(i int) Simplex` | Face returns the codimension-1 face obtained by deleting the i-th vertex (0-indexed, ascending order) of the simplex. |
| `func (s Simplex) Faces() []Simplex` | Faces returns the codimension-1 faces of the simplex: for a k-simplex the k+1 sub-simplices obtained by omitting one vertex. |
| `func (s Simplex) Intersect(t Simplex) Simplex` | Intersect returns the simplex on the common vertices of s and t. |
| `func (s Simplex) IsEmpty() bool` | IsEmpty reports whether the simplex has no vertices. |
| `func (s Simplex) IsFaceOf(t Simplex) bool` | IsFaceOf reports whether s is a face of t (equivalent to t.ContainsSimplex(s)). |
| `func (s Simplex) IsProperFaceOf(t Simplex) bool` | IsProperFaceOf reports whether s is a face of t of strictly lower dimension. |
| `func (s Simplex) Join(t Simplex) Simplex` | Join returns the simplex whose vertex set is the union of those of s and t. |
| `func (s Simplex) Key() string` | Key returns a canonical string key identifying the simplex, suitable as a map key. |
| `func (s Simplex) Len() int` | Len returns the number of vertices of the simplex. |
| `func (s Simplex) Max() int` | Max returns the largest vertex, or 0 for the empty simplex. |
| `func (s Simplex) Min() int` | Min returns the smallest vertex, or 0 for the empty simplex. |
| `func (s Simplex) ProperFaces() []Simplex` | ProperFaces returns every proper face of the simplex (all non-empty subsets except s itself). |
| `func (s Simplex) String() string` | String returns a human-readable representation such as "{0,1,2}". |
| `func (s Simplex) Subsets() []Simplex` | Subsets returns every non-empty subset of the vertices of s as a simplex, including s itself. |
| `func (s Simplex) Vertex(i int) int` | Vertex returns the i-th vertex (0-indexed) of the simplex in ascending order. |
| `func (s Simplex) Vertices() []int` | Vertices returns a copy of the simplex's vertex list in ascending order. |

</details>

<details>
<summary><code>SmithResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (s *SmithResult) Diagonal() []*big.Int` | Diagonal returns the diagonal entries of D, i.e. |
| `func (s *SmithResult) InvariantFactors() []*big.Int` | InvariantFactors returns the non-zero diagonal entries of D, the invariant factors of the matrix. |
| `func (s *SmithResult) Rank() int` | Rank returns the number of non-zero diagonal entries of D. |
| `func (s *SmithResult) TorsionFactors() []*big.Int` | TorsionFactors returns the invariant factors of D that are strictly greater than one; these are the elementary divisors giving rise to torsion. |

</details>

### Variables

`ErrDimMismatch`, `ErrDimension`, `ErrEmptyCloud`, `ErrEmptySimplex`, `ErrNotSquare`, `ErrShape`, `ErrSingular`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/simplicial).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
