# matroids — implements matroid theory and matroid-optimization algorithms using only the Go

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/matroids.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/matroids)

Package matroids implements matroid theory and matroid-optimization algorithms
using only the Go standard library.

A matroid is a combinatorial structure that abstracts the notion of linear
independence. Every matroid in this package has a finite ground set whose
elements are the integers 0, 1, ..., n-1, where n is reported by
`Matroid.Size`. The single defining operation is the rank function
`Matroid.Rank`, which maps a subset of the ground set to a non-negative
integer satisfying the matroid rank axioms (normalisation, monotonicity, unit
increase and submodularity). All higher-level notions — independence, bases,
circuits, closure, flats, loops, coloops and connectivity — are derived from
the rank oracle by the generic free functions in this package, so they work
uniformly across every matroid representation.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/matroids"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func Bases(m Matroid) [][]int` | Bases returns every basis of m, sorted lexicographically. |
| `func BasisCount(m Matroid) int` | BasisCount returns the number of bases of m. |
| `func Cardinality(a []int) int` | Cardinality returns the number of distinct elements in a. |
| `func CheckRankAxioms(m Matroid) error` | CheckRankAxioms verifies, by brute force over all subsets, that m satisfies the matroid rank axioms: 0 ≤ rank(X) ≤ \|X\|; rank is monotone; and… |
| `func Circuits(m Matroid) [][]int` | Circuits returns every circuit of m, sorted by size then lexicographically. |
| `func Closure(m Matroid, set []int) []int` | Closure returns the closure (span) of set: every element e of the ground set such that rank(set ∪ {e}) = rank(set). |
| `func Cobases(m Matroid) [][]int` | Cobases returns the cobases of m: the bases of its dual, i.e. |
| `func Cocircuits(m Matroid) [][]int` | Cocircuits returns the cocircuits of m: the circuits of its dual, equivalently the minimal sets meeting every basis of m. |
| `func Coloops(m Matroid) []int` | Coloops returns the sorted list of coloops of m. |
| `func Complement(n int, a []int) []int` | Complement returns the sorted elements of {0,...,n-1} that are not in a. |
| `func ComponentOf(m Matroid, e int) []int` | ComponentOf returns the connected component of m that contains element e. |
| `func Components(m Matroid) [][]int` | Components returns the connected components of m as a partition of the ground set. |
| `func ConnectivityFunction(m Matroid, set []int) int` | ConnectivityFunction returns the matroid connectivity function λ(S) = r(S) + r(E\S) - r(E), where r is the rank function of m and E the ground set. |
| `func CopySet(a []int) []int` | CopySet returns a copy of a. |
| `func Corank(m Matroid, set []int) int` | Corank returns the rank of the ground set minus the rank of set, i.e. |
| `func Dependent(m Matroid, set []int) bool` | Dependent reports whether set is dependent (not independent) in m. |
| `func DependentSets(m Matroid) [][]int` | DependentSets returns every dependent set of m. |
| `func DistinctSorted(a []int) []int` | DistinctSorted returns the distinct elements of a in ascending order. |
| `func Flats(m Matroid) [][]int` | Flats returns every flat (closed set) of m, sorted by rank then lexicographically. |
| `func FlatsOfRank(m Matroid, r int) [][]int` | FlatsOfRank returns every flat of m whose rank equals r. |
| `func FullRank(m Matroid) int` | FullRank returns the rank of the whole ground set of m, i.e. |
| `func FundamentalCircuit(m Matroid, basis []int, e int) []int` | FundamentalCircuit returns the unique circuit contained in basis ∪ {e}, where basis is an independent set and e is an element not in basis whose… |
| `func FundamentalCocircuit(m Matroid, basis []int, e int) []int` | FundamentalCocircuit returns the fundamental cocircuit of element e with respect to basis: the unique cocircuit contained in (E \ basis) ∪ {e},… |
| `func Girth(m Matroid) int` | Girth returns the size of a smallest circuit of m. |
| `func Ground(m Matroid) []int` | Ground returns the ground set {0, 1, ..., m.Size()-1}. |
| `func Hyperplanes(m Matroid) [][]int` | Hyperplanes returns the hyperplanes of m: the flats of rank FullRank(m)-1. |
| `func Independent(m Matroid, set []int) bool` | Independent reports whether set is an independent set of m, i.e. |
| `func IndependentSets(m Matroid) [][]int` | IndependentSets returns every independent set of m, sorted by size then lexicographically. |
| `func Intersection(m1, m2 Matroid) []int` | Intersection returns a maximum-cardinality set that is independent in both m1 and m2, computed by repeatedly finding a shortest augmenting path in… |
| `func IntersectionSize(m1, m2 Matroid) int` | IntersectionSize returns the size of a maximum common independent set of m1 and m2. |
| `func IsBasis(m Matroid, set []int) bool` | IsBasis reports whether set is a basis of m: a maximal independent set, equivalently an independent spanning set. |
| `func IsCircuit(m Matroid, set []int) bool` | IsCircuit reports whether set is a circuit of m: a minimal dependent set. |
| `func IsClosed(m Matroid, set []int) bool` | IsClosed reports whether set equals its own closure in m; such sets are the flats of m. |
| `func IsColoop(m Matroid, e int) bool` | IsColoop reports whether element e is a coloop (isthmus, bridge) of m, i.e. |
| `func IsConnected(m Matroid) bool` | IsConnected reports whether m is connected: it has at most one component (every pair of elements lies in a common circuit). |
| `func IsFlat(m Matroid, set []int) bool` | IsFlat is a synonym for `IsClosed`. |
| `func IsLoop(m Matroid, e int) bool` | IsLoop reports whether element e is a loop of m, i.e. |
| `func IsMatroidIndependenceSystem(n int, sets [][]int) bool` | IsMatroidIndependenceSystem reports whether the given family of subsets of {0..n-1} satisfies the matroid axioms (non-empty, downward closed,… |
| `func IsPartitionable(parts []Matroid) bool` | IsPartitionable reports whether the whole ground set of the given matroids can be partitioned into independent sets, one per matroid; equivalently… |
| `func IsSeparator(m Matroid, set []int) bool` | IsSeparator reports whether set is a separator of m, i.e. |
| `func IsSimple(m Matroid) bool` | IsSimple reports whether m is a simple matroid: it has no loops and no two distinct parallel elements. |
| `func IsSpanning(m Matroid, set []int) bool` | IsSpanning reports whether set spans m, i.e. |
| `func IsSubset(a, b []int) bool` | IsSubset reports whether every distinct element of a is also in b. |
| `func IsUniform(m Matroid) bool` | IsUniform reports whether m has the rank profile of a uniform matroid, i.e. |
| `func Loops(m Matroid) []int` | Loops returns the sorted list of loops of m. |
| `func MaxWeightIndependentSet(m Matroid, weights []float64) []int` | MaxWeightIndependentSet returns the set found by `Greedy`; it is a convenience wrapper returning only the element list. |
| `func Nullity(m Matroid, set []int) int` | Nullity returns \|set\| - rank(set), the number of elements of set that are dependent on the rest. |
| `func NumComponents(m Matroid) int` | NumComponents returns the number of connected components of m. |
| `func Parallel(m Matroid, e, f int) bool` | Parallel reports whether distinct non-loop elements e and f are parallel in m, i.e. |
| `func ParallelClasses(m Matroid) [][]int` | ParallelClasses returns the partition of the non-loop elements into maximal classes of mutually parallel elements. |
| `func RangeSet(n int) []int` | RangeSet returns the slice {0, 1, ..., n-1}. |
| `func RankProfile(m Matroid) []int` | RankProfile returns a slice p where p`k` is the number of independent sets of m of size k, for k from 0 to FullRank(m). |
| `func SetContains(a []int, e int) bool` | SetContains reports whether a contains e. |
| `func SetDifference(a, b []int) []int` | SetDifference returns the sorted set of elements in a but not in b. |
| `func SetEqual(a, b []int) bool` | SetEqual reports whether a and b contain the same distinct elements. |
| `func SetInsert(a []int, e int) []int` | SetInsert returns the sorted set a with e added. |
| `func SetIntersection(a, b []int) []int` | SetIntersection returns the sorted intersection of a and b. |
| `func SetRemove(a []int, e int) []int` | SetRemove returns the sorted set a with e removed. |
| `func SetUnion(a, b []int) []int` | SetUnion returns the sorted union of a and b as distinct elements. |
| `func SpanningSets(m Matroid) [][]int` | SpanningSets returns every spanning set of m. |
| `func Spans(m Matroid, set []int, e int) bool` | Spans reports whether adding e to set does not increase its rank, i.e. |
| `func Subsets(a []int) [][]int` | Subsets returns every subset of the distinct elements of a. |
| `func SubsetsOfSize(a []int, k int) [][]int` | SubsetsOfSize returns every size-k subset of the distinct elements of a. |
| `func UnionMaxIndependentSet(parts []Matroid) []int` | UnionMaxIndependentSet returns a maximum-size set independent in the union of the given matroids. |
| `func UnionRank(parts []Matroid) int` | UnionRank returns the rank of the union of the given matroids (the maximum size of a set partitionable into independents of the parts). |
| `func WeightedIntersection(m1, m2 Matroid, weights []float64) []int` | WeightedIntersection returns a maximum-weight set that is independent in both m1 and m2. |
| `func WeightedIntersectionWeight(m1, m2 Matroid, weights []float64) float64` | WeightedIntersectionWeight returns the total weight of the set returned by `WeightedIntersection`. |

### Types

| Type | What it is |
| --- | --- |
| `BinaryMatroid` | BinaryMatroid is the column matroid of a 0/1 matrix over the two-element field GF(2). |
| `DirectSumMatroid` | DirectSumMatroid is the direct sum (1-sum) of two or more matroids on disjoint ground sets. |
| `DualMatroid` | DualMatroid is the dual of a matroid m on the same ground set. |
| `ExplicitMatroid` | ExplicitMatroid is a matroid given directly by the collection of its independent sets. |
| `GraphicMatroid` | GraphicMatroid is the cycle matroid (polygon matroid) of an undirected multigraph. |
| `GreedyResult` | GreedyResult holds the outcome of a greedy matroid optimisation: the selected set (sorted) and its total weight. |
| `IntSet` | IntSet is a mutable set of integers implemented as a map. |
| `LinearMatroid` | LinearMatroid is the column matroid of an exact rational matrix. |
| `Matroid` | Matroid is the interface satisfied by every matroid representation in this package. |
| `MinorMatroid` | MinorMatroid is a minor of a base matroid, obtained by deleting a set D and contracting a set C of ground elements. |
| `PartitionMatroid` | PartitionMatroid is a matroid whose ground set is partitioned into blocks, each with an independence capacity: a set is independent when it contains… |
| `TransversalMatroid` | TransversalMatroid is the transversal matroid of a family of sets. |
| `UniformMatroid` | UniformMatroid is the uniform matroid U(r, n): its ground set has n elements and a subset is independent exactly when it has at most r elements. |
| `UnionMatroid` | UnionMatroid is the matroid union M1 ∨ M2 ∨ ... |

<details>
<summary><code>BinaryMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewBinaryMatroid(matrix [][]int) *BinaryMatroid` | NewBinaryMatroid builds a GF(2) linear matroid from a rows×cols 0/1 matrix in row-major order. |
| `func (m *BinaryMatroid) NumRows() int` | NumRows returns the number of rows of the matrix. |
| `func (m *BinaryMatroid) Rank(set []int) int` | Rank returns the GF(2) rank of the columns in set, computed by Gaussian elimination over GF(2). |
| `func (m *BinaryMatroid) Size() int` | Size returns the number of columns (ground-set elements). |

</details>

<details>
<summary><code>DirectSumMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func DirectSum(parts ...Matroid) *DirectSumMatroid` | DirectSum returns the direct sum of the given matroids. |
| `func (ds *DirectSumMatroid) NumParts() int` | NumParts returns the number of summands. |
| `func (ds *DirectSumMatroid) PartOf(e int) (part, local int)` | PartOf returns the part index and the local (within-part) index of global element e. |
| `func (ds *DirectSumMatroid) Rank(set []int) int` | Rank returns the sum over parts of the rank of the elements of set falling in that part (translated to local indices). |
| `func (ds *DirectSumMatroid) Size() int` | Size returns the total number of elements across all parts. |

</details>

<details>
<summary><code>DualMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Dual(m Matroid) *DualMatroid` | Dual returns the dual matroid of m. |
| `func (d *DualMatroid) Primal() Matroid` | Primal returns the underlying matroid whose dual this is. |
| `func (d *DualMatroid) Rank(set []int) int` | Rank returns the dual rank r*(S) = \|S\| + r(E\S) - r(E). |
| `func (d *DualMatroid) Size() int` | Size returns the ground-set size, the same as the primal matroid. |

</details>

<details>
<summary><code>ExplicitMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewExplicitMatroidFromBases(n int, bases [][]int) *ExplicitMatroid` | NewExplicitMatroidFromBases builds a matroid on n elements whose bases are exactly the given sets. |
| `func NewExplicitMatroidFromCircuits(n int, circuits [][]int) *ExplicitMatroid` | NewExplicitMatroidFromCircuits builds a matroid on n elements from its circuits: a set is independent exactly when it contains no circuit. |
| `func NewExplicitMatroidFromIndependentSets(n int, sets [][]int) *ExplicitMatroid` | NewExplicitMatroidFromIndependentSets builds a matroid on n elements whose independent sets are exactly the given sets, downward-closed… |
| `func (m *ExplicitMatroid) IsIndependentSet(set []int) bool` | IsIndependentSet reports whether the exact set (as a mask) is recorded as independent. |
| `func (m *ExplicitMatroid) ListIndependentSets() [][]int` | ListIndependentSets returns every recorded independent set of m, sorted by size then lexicographically. |
| `func (m *ExplicitMatroid) Rank(set []int) int` | Rank returns the size of a largest independent subset of set, found by a downward search over the recorded independent sets. |
| `func (m *ExplicitMatroid) Size() int` | Size returns the number of ground-set elements. |
| `func (m *ExplicitMatroid) Validate() error` | Validate checks that the recorded independent sets form a matroid: the empty set is independent, independence is downward closed, and the exchange… |

</details>

<details>
<summary><code>GraphicMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func CycleMatroidOfCompleteGraph(n int) *GraphicMatroid` | CycleMatroidOfCompleteGraph returns the graphic matroid of the complete graph K_n. |
| `func CycleMatroidOfCycleGraph(n int) *GraphicMatroid` | CycleMatroidOfCycleGraph returns the graphic matroid of the cycle graph C_n (n ≥ 1 vertices joined in a single ring). |
| `func NewGraphicMatroid(verts int) *GraphicMatroid` | NewGraphicMatroid returns an empty graphic matroid on the given number of vertices (labeled 0..verts-1) with no edges. |
| `func NewGraphicMatroidFromEdges(verts int, edges [][2]int) *GraphicMatroid` | NewGraphicMatroidFromEdges returns the graphic matroid on verts vertices with the given edges. |
| `func (m *GraphicMatroid) AddEdge(u, v int) int` | AddEdge appends an edge {u, v} and returns its element index. |
| `func (m *GraphicMatroid) Edge(e int) (int, int)` | Edge returns the endpoints of edge e. |
| `func (m *GraphicMatroid) IsForest(set []int) bool` | IsForest reports whether the given edge set is independent in m (contains no cycle). |
| `func (m *GraphicMatroid) NumVertices() int` | NumVertices returns the number of vertices of the underlying graph. |
| `func (m *GraphicMatroid) Rank(set []int) int` | Rank returns verts - (number of connected components of the spanning subgraph induced by the edges in set). |
| `func (m *GraphicMatroid) Size() int` | Size returns the number of edges (ground-set elements). |

</details>

<details>
<summary><code>GreedyResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Greedy(m Matroid, weights []float64) GreedyResult` | Greedy runs the matroid greedy algorithm for a maximum-weight independent set. |
| `func GreedyMaxWeightBasis(m Matroid, weights []float64) GreedyResult` | GreedyMaxWeightBasis returns a maximum-weight basis of m using the greedy algorithm: elements are considered in order of non-increasing weight and… |
| `func GreedyMinWeightBasis(m Matroid, weights []float64) GreedyResult` | GreedyMinWeightBasis returns a minimum-weight basis of m. |

</details>

<details>
<summary><code>IntSet</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewIntSet(elems ...int) IntSet` | NewIntSet returns a new `IntSet` containing the given elements. |
| `func (s IntSet) Add(e int)` | Add inserts e into the set. |
| `func (s IntSet) AddAll(elems []int)` | AddAll inserts every element of elems into the set. |
| `func (s IntSet) Clone() IntSet` | Clone returns an independent copy of the set. |
| `func (s IntSet) Contains(e int) bool` | Contains reports whether e is a member of the set. |
| `func (s IntSet) Difference(t IntSet) IntSet` | Difference returns a new set containing the elements of s that are not in t. |
| `func (s IntSet) Equal(t IntSet) bool` | Equal reports whether s and t contain exactly the same elements. |
| `func (s IntSet) Intersect(t IntSet) IntSet` | Intersect returns a new set containing the elements common to s and t. |
| `func (s IntSet) IsSubsetOf(t IntSet) bool` | IsSubsetOf reports whether every element of s is also in t. |
| `func (s IntSet) Len() int` | Len returns the number of elements in the set. |
| `func (s IntSet) Remove(e int)` | Remove deletes e from the set if present. |
| `func (s IntSet) Slice() []int` | Slice returns the elements of the set in ascending order. |
| `func (s IntSet) Union(t IntSet) IntSet` | Union returns a new set containing every element of s or t. |

</details>

<details>
<summary><code>LinearMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewLinearMatroid(matrix [][]int64) *LinearMatroid` | NewLinearMatroid builds a linear matroid from a rows×cols matrix of int64 entries given in row-major order. |
| `func NewLinearMatroidRat(matrix [][]*big.Rat) *LinearMatroid` | NewLinearMatroidRat builds a linear matroid from a rows×cols matrix of *big.Rat entries in row-major order. |
| `func (m *LinearMatroid) Entry(i, j int) *big.Rat` | Entry returns a copy of the matrix entry in row i, column j. |
| `func (m *LinearMatroid) NumRows() int` | NumRows returns the number of rows of the matrix. |
| `func (m *LinearMatroid) Rank(set []int) int` | Rank returns the rank of the submatrix formed by the columns in set, computed by exact Gaussian elimination. |
| `func (m *LinearMatroid) Size() int` | Size returns the number of columns (ground-set elements). |

</details>

<details>
<summary><code>MinorMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Contraction(m Matroid, con []int) *MinorMatroid` | Contraction returns the minor of m obtained by contracting the elements in con. |
| `func Deletion(m Matroid, del []int) *MinorMatroid` | Deletion returns the minor of m obtained by deleting the elements in del. |
| `func Minor(m Matroid, del, con []int) *MinorMatroid` | Minor returns the minor of m obtained by deleting the elements in del and contracting the elements in con. |
| `func Restriction(m Matroid, keep []int) *MinorMatroid` | Restriction returns the restriction of m to the elements in keep, i.e. |
| `func (mm *MinorMatroid) Elements() []int` | Elements returns the original ground-set labels of the minor's elements, in new-index order (Elements()`i` is the original label of new element i). |
| `func (mm *MinorMatroid) OriginalLabel(i int) int` | OriginalLabel maps a new element index of the minor back to its original label in the base matroid. |
| `func (mm *MinorMatroid) Rank(set []int) int` | Rank returns the rank of the given set of new-indexed elements in the minor. |
| `func (mm *MinorMatroid) Size() int` | Size returns the number of elements of the minor. |

</details>

<details>
<summary><code>PartitionMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewPartitionMatroid(n int, blockOf, caps []int) *PartitionMatroid` | NewPartitionMatroid builds a partition matroid on n elements. |
| `func NewPartitionMatroidFromBlocks(n int, blocks [][]int, caps []int) *PartitionMatroid` | NewPartitionMatroidFromBlocks builds a partition matroid from an explicit list of blocks (each a slice of element indices in [0, n)) and matching… |
| `func (m *PartitionMatroid) BlockOf(e int) int` | BlockOf returns the block index of element e. |
| `func (m *PartitionMatroid) Capacity(b int) int` | Capacity returns the independence capacity of block b. |
| `func (m *PartitionMatroid) NumBlocks() int` | NumBlocks returns the number of blocks in the partition. |
| `func (m *PartitionMatroid) Rank(set []int) int` | Rank returns the sum over blocks of min(count in block, capacity). |
| `func (m *PartitionMatroid) Size() int` | Size returns the number of ground-set elements. |

</details>

<details>
<summary><code>TransversalMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewTransversalMatroid(n int, adj [][]int) *TransversalMatroid` | NewTransversalMatroid builds a transversal matroid on n ground elements. |
| `func (m *TransversalMatroid) NumTargets() int` | NumTargets returns the number of distinct target (right) vertices. |
| `func (m *TransversalMatroid) Rank(set []int) int` | Rank returns the size of a maximum matching that saturates a subset of the distinct in-range elements of set, i.e. |
| `func (m *TransversalMatroid) Size() int` | Size returns the number of ground-set elements. |
| `func (m *TransversalMatroid) SystemOfDistinctRepresentatives(set []int) (map[int]int, bool)` | SystemOfDistinctRepresentatives returns, for an independent set S (given as set), an assignment mapping each element of S to a distinct target index,… |
| `func (m *TransversalMatroid) Targets(e int) []int` | Targets returns the internal target indices associated with element e. |

</details>

<details>
<summary><code>UniformMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func FreeMatroid(n int) *UniformMatroid` | FreeMatroid is an alias constructor returning the free matroid on n elements; it is provided for symmetry with the other named constructions. |
| `func NewFreeMatroid(n int) *UniformMatroid` | NewFreeMatroid returns the free matroid on n elements, U(n, n), in which every subset is independent. |
| `func NewTrivialMatroid(n int) *UniformMatroid` | NewTrivialMatroid returns the rank-zero matroid on n elements, U(0, n), in which every element is a loop and the only independent set is empty. |
| `func NewUniformMatroid(r, n int) *UniformMatroid` | NewUniformMatroid returns the uniform matroid U(r, n) on n elements with rank r. |
| `func TrivialMatroid(n int) *UniformMatroid` | TrivialMatroid is an alias constructor returning the rank-zero matroid on n elements. |
| `func (m *UniformMatroid) Rank(set []int) int` | Rank returns min(\|distinct in-range elements of set\|, r). |
| `func (m *UniformMatroid) Size() int` | Size returns the number of ground-set elements, n. |
| `func (m *UniformMatroid) UniformRank() int` | UniformRank returns the parameter r of U(r, n), i.e. |

</details>

<details>
<summary><code>UnionMatroid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Union(parts ...Matroid) *UnionMatroid` | Union returns the matroid union of the given matroids, which must all share a ground set (equal Size). |
| `func (u *UnionMatroid) MaxIndependentSet() []int` | MaxIndependentSet returns a maximum-size subset of the whole ground set that is independent in the union (partitionable into independents of the… |
| `func (u *UnionMatroid) NumParts() int` | NumParts returns the number of matroids being combined. |
| `func (u *UnionMatroid) Partition() [][]int` | Partition returns an assignment of a maximum union-independent set into the parts: result`i` is the (sorted) list of ground elements assigned to… |
| `func (u *UnionMatroid) Rank(set []int) int` | Rank returns the union rank of set: the maximum size of a subset of set that is partitionable into independent sets of the parts. |
| `func (u *UnionMatroid) Size() int` | Size returns the size of the shared ground set. |

</details>

### Variables

`ErrDimensionMismatch`, `ErrInvalidElement`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/matroids).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
