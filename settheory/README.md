# settheory — implements finite set theory, binary relations

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/settheory.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/settheory)

Package settheory implements finite set theory, binary relations, order theory
and lattices over integer and string domains, built entirely on the Go
standard library.

The package is self-contained and deterministic: identical inputs always yield
identical outputs, iteration-order sensitivity is eliminated by returning
sorted slices from every accessor, and no third-party dependency is required.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/settheory"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func BellNumber(n int) int` | BellNumber returns the nth Bell number, the number of partitions of an n-set, computed exactly with the Bell triangle. |
| `func StirlingSecond(n, k int) int` | StirlingSecond returns the Stirling number of the second kind S(n, k), the number of ways to partition an n-set into exactly k non-empty blocks. |

### Types

| Type | What it is |
| --- | --- |
| `IntSet` | IntSet is a finite set of int values backed by a hash map. |
| `Pair` | Pair is an ordered pair of integers (From, To). |
| `Partition` | Partition is a partition of a finite set of integers into non-empty, pairwise-disjoint blocks whose union is the whole set. |
| `Poset` | Poset is a finite partially ordered set: a set of integer elements together with a reflexive, antisymmetric and transitive "less than or equal"… |
| `Relation` | Relation is a binary relation on the integers, represented as a set of ordered pairs. |
| `StringSet` | StringSet is a finite set of string values backed by a hash map. |

<details>
<summary><code>IntSet</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewIntSet(elems ...int) IntSet` | NewIntSet returns a new IntSet containing the given elements. |
| `func (s IntSet) Add(x int)` | Add inserts x into the set. |
| `func (s IntSet) CartesianProduct(other IntSet) []Pair` | CartesianProduct returns every ordered pair (a, b) with a in s and b in other, sorted lexicographically. |
| `func (s IntSet) Clone() IntSet` | Clone returns an independent copy of the set that shares no storage with the receiver. |
| `func (s IntSet) Contains(x int) bool` | Contains reports whether x is a member of the set. |
| `func (s IntSet) Difference(other IntSet) IntSet` | Difference returns a new set containing every element that is in s but not in other. |
| `func (s IntSet) Elements() []int` | Elements returns the members of the set in ascending order. |
| `func (s IntSet) Equal(other IntSet) bool` | Equal reports whether s and other contain exactly the same elements. |
| `func (s IntSet) Intersection(other IntSet) IntSet` | Intersection returns a new set containing every element that is in both s and other. |
| `func (s IntSet) IsDisjoint(other IntSet) bool` | IsDisjoint reports whether s and other share no elements. |
| `func (s IntSet) IsEmpty() bool` | IsEmpty reports whether the set has no elements. |
| `func (s IntSet) IsProperSubsetOf(other IntSet) bool` | IsProperSubsetOf reports whether s is a subset of other and the two sets are not equal. |
| `func (s IntSet) IsSubsetOf(other IntSet) bool` | IsSubsetOf reports whether every element of s is also an element of other. |
| `func (s IntSet) IsSupersetOf(other IntSet) bool` | IsSupersetOf reports whether every element of other is also an element of s. |
| `func (s IntSet) Len() int` | Len returns the cardinality of the set. |
| `func (s IntSet) Max() (int, bool)` | Max returns the largest element of the set. |
| `func (s IntSet) Min() (int, bool)` | Min returns the smallest element of the set. |
| `func (s IntSet) PowerSet() []IntSet` | PowerSet returns all 2^n subsets of the set, where n is its cardinality. |
| `func (s IntSet) Remove(x int)` | Remove deletes x from the set. |
| `func (s IntSet) String() string` | String returns a canonical textual representation of the set such as "{1, 2, 3}", with elements in ascending order. |
| `func (s IntSet) Sum() int` | Sum returns the sum of all elements of the set; the empty set sums to zero. |
| `func (s IntSet) SymmetricDifference(other IntSet) IntSet` | SymmetricDifference returns a new set containing every element that is in exactly one of s and other. |
| `func (s IntSet) Union(other IntSet) IntSet` | Union returns a new set containing every element that is in s or in other. |

</details>

<details>
<summary><code>Pair</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func (p Pair) String() string` | String renders the pair as "(from, to)". |

</details>

<details>
<summary><code>Partition</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func EquivalenceClasses(rel Relation, set IntSet) Partition` | EquivalenceClasses returns the partition of set into the equivalence classes of the relation rel, which must be an equivalence relation on set. |
| `func (p Partition) Blocks() []IntSet` | Blocks returns the blocks of the partition sorted by ascending least element, each block itself presented via IntSet.Elements when rendered. |
| `func (p Partition) IsValidPartitionOf(set IntSet) bool` | IsValidPartitionOf reports whether p is a valid partition of set: every block is non-empty, blocks are pairwise disjoint, and their union equals set. |
| `func (p Partition) Len() int` | Len returns the number of blocks in the partition. |
| `func (p Partition) Refines(q Partition) bool` | Refines reports whether p is a refinement of q, i.e. |
| `func (p Partition) Universe() IntSet` | Universe returns the union of all blocks, i.e. |

</details>

<details>
<summary><code>Poset</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewChainPoset(elems ...int) Poset` | NewChainPoset returns the totally ordered Poset on the given elements ordered by ascending integer value. |
| `func NewDivisibilityPoset(n int) (Poset, error)` | NewDivisibilityPoset returns the divisibility Poset on the positive divisors of n, where a ≤ b means a divides b. |
| `func NewPoset(set IntSet, leq Relation) (Poset, error)` | NewPoset returns a Poset over set whose order is leq. |
| `func (p Poset) AreComparable(a, b int) bool` | AreComparable reports whether a ≤ b or b ≤ a. |
| `func (p Poset) Covers(lower, upper int) bool` | Covers reports whether upper covers lower: lower < upper and there is no element c with lower < c < upper. |
| `func (p Poset) Elements() []int` | Elements returns the elements of the poset in ascending order. |
| `func (p Poset) GreatestElement() (int, bool)` | GreatestElement returns the unique element that is ≥ every element of the poset, if one exists. |
| `func (p Poset) HasseEdges() []Pair` | HasseEdges returns the covering pairs of the poset as Pairs (From = lower, To = upper) with To covering From. |
| `func (p Poset) Height() int` | Height returns the number of elements in a longest chain of the poset. |
| `func (p Poset) IsAntichain(subset []int) bool` | IsAntichain reports whether subset is an antichain: no two distinct elements of subset are comparable. |
| `func (p Poset) IsChain() bool` | IsChain reports whether the poset is totally ordered: every two elements are comparable. |
| `func (p Poset) IsLattice() bool` | IsLattice reports whether every pair of elements has both a meet and a join, i.e. |
| `func (p Poset) Join(a, b int) (int, bool)` | Join returns the least upper bound (supremum) of a and b, if it exists. |
| `func (p Poset) LeastElement() (int, bool)` | LeastElement returns the unique element that is ≤ every element of the poset, if one exists. |
| `func (p Poset) Leq(a, b int) bool` | Leq reports whether a ≤ b in the poset order. |
| `func (p Poset) Less(a, b int) bool` | Less reports whether a < b, i.e. |
| `func (p Poset) LowerBounds(subset []int) []int` | LowerBounds returns the elements that are ≤ every element of subset, in ascending order. |
| `func (p Poset) MaximalElements() []int` | MaximalElements returns the elements that have no strictly larger element, in ascending order. |
| `func (p Poset) Meet(a, b int) (int, bool)` | Meet returns the greatest lower bound (infimum) of a and b, if it exists. |
| `func (p Poset) MinimalElements() []int` | MinimalElements returns the elements that have no strictly smaller element, in ascending order. |
| `func (p Poset) Size() int` | Size returns the number of elements in the poset. |
| `func (p Poset) TopologicalOrder() []int` | TopologicalOrder returns a linear extension of the poset: a permutation of the elements in which a appears before b whenever a < b. |
| `func (p Poset) UpperBounds(subset []int) []int` | UpperBounds returns the elements that are ≥ every element of subset, in ascending order. |

</details>

<details>
<summary><code>Relation</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewRelation(pairs ...Pair) Relation` | NewRelation returns a new Relation containing the given pairs. |
| `func RelationFromPairs(pairs [][2]int) Relation` | RelationFromPairs builds a Relation from ordered pairs supplied as two-element arrays [from, to]. |
| `func RelationFromPartition(p Partition) Relation` | RelationFromPartition returns the equivalence relation induced by p: two elements are related exactly when they lie in the same block. |
| `func (r Relation) Add(a, b int)` | Add inserts the ordered pair (a, b) into the relation. |
| `func (r Relation) Clone() Relation` | Clone returns an independent copy of the relation. |
| `func (r Relation) Compose(other Relation) Relation` | Compose returns the composition r∘other defined so that (a, c) is in the result whenever (a, b) is in other and (b, c) is in r for some b. |
| `func (r Relation) Domain() IntSet` | Domain returns the set of first coordinates appearing in the relation. |
| `func (r Relation) Equal(other Relation) bool` | Equal reports whether r and other contain exactly the same pairs. |
| `func (r Relation) EquivalenceClosureOn(set IntSet) Relation` | EquivalenceClosureOn returns the smallest equivalence relation over set that contains r. |
| `func (r Relation) Field() IntSet` | Field returns the union of the domain and range, i.e. |
| `func (r Relation) Inverse() Relation` | Inverse returns the converse relation, containing (b, a) for every pair (a, b) in r. |
| `func (r Relation) IsAntisymmetric() bool` | IsAntisymmetric reports whether the presence of both (a, b) and (b, a) implies a equals b. |
| `func (r Relation) IsEquivalenceOn(set IntSet) bool` | IsEquivalenceOn reports whether the relation is reflexive on set, symmetric and transitive, i.e. |
| `func (r Relation) IsIrreflexiveOn(set IntSet) bool` | IsIrreflexiveOn reports whether (x, x) is absent for every x in set. |
| `func (r Relation) IsPartialOrderOn(set IntSet) bool` | IsPartialOrderOn reports whether the relation is reflexive on set, antisymmetric and transitive, i.e. |
| `func (r Relation) IsReflexiveOn(set IntSet) bool` | IsReflexiveOn reports whether (x, x) is in the relation for every x in set. |
| `func (r Relation) IsSymmetric() bool` | IsSymmetric reports whether (b, a) is in the relation whenever (a, b) is. |
| `func (r Relation) IsTransitive() bool` | IsTransitive reports whether (a, c) is in the relation whenever both (a, b) and (b, c) are. |
| `func (r Relation) Len() int` | Len returns the number of ordered pairs in the relation. |
| `func (r Relation) Pairs() []Pair` | Pairs returns the ordered pairs of the relation sorted lexicographically by (From, To). |
| `func (r Relation) Range() IntSet` | Range returns the set of second coordinates appearing in the relation. |
| `func (r Relation) ReflexiveClosureOn(set IntSet) Relation` | ReflexiveClosureOn returns the smallest relation containing r that is reflexive on set. |
| `func (r Relation) Related(a, b int) bool` | Related reports whether the ordered pair (a, b) is in the relation. |
| `func (r Relation) Remove(a, b int)` | Remove deletes the ordered pair (a, b) from the relation. |
| `func (r Relation) SymmetricClosure() Relation` | SymmetricClosure returns the smallest symmetric relation containing r. |
| `func (r Relation) TransitiveClosure() Relation` | TransitiveClosure returns the smallest transitive relation containing r, computed with Warshall's algorithm over the field of the relation. |

</details>

<details>
<summary><code>StringSet</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewStringSet(elems ...string) StringSet` | NewStringSet returns a new StringSet containing the given elements. |
| `func (s StringSet) Add(x string)` | Add inserts x into the set. |
| `func (s StringSet) Clone() StringSet` | Clone returns an independent copy of the set that shares no storage with the receiver. |
| `func (s StringSet) Contains(x string) bool` | Contains reports whether x is a member of the set. |
| `func (s StringSet) Difference(other StringSet) StringSet` | Difference returns a new set containing every element that is in s but not in other. |
| `func (s StringSet) Elements() []string` | Elements returns the members of the set in ascending lexicographic order. |
| `func (s StringSet) Equal(other StringSet) bool` | Equal reports whether s and other contain exactly the same elements. |
| `func (s StringSet) Intersection(other StringSet) StringSet` | Intersection returns a new set containing every element that is in both s and other. |
| `func (s StringSet) IsDisjoint(other StringSet) bool` | IsDisjoint reports whether s and other share no elements. |
| `func (s StringSet) IsEmpty() bool` | IsEmpty reports whether the set has no elements. |
| `func (s StringSet) IsProperSubsetOf(other StringSet) bool` | IsProperSubsetOf reports whether s is a subset of other and the two sets are not equal. |
| `func (s StringSet) IsSubsetOf(other StringSet) bool` | IsSubsetOf reports whether every element of s is also an element of other. |
| `func (s StringSet) IsSupersetOf(other StringSet) bool` | IsSupersetOf reports whether every element of other is also an element of s. |
| `func (s StringSet) Len() int` | Len returns the cardinality of the set. |
| `func (s StringSet) PowerSet() []StringSet` | PowerSet returns all 2^n subsets of the set, where n is its cardinality. |
| `func (s StringSet) Remove(x string)` | Remove deletes x from the set. |
| `func (s StringSet) String() string` | String returns a canonical textual representation of the set such as "{a, b, c}", with elements quoted and in ascending order. |
| `func (s StringSet) SymmetricDifference(other StringSet) StringSet` | SymmetricDifference returns a new set containing every element that is in exactly one of s and other. |
| `func (s StringSet) Union(other StringSet) StringSet` | Union returns a new set containing every element that is in s or in other. |

</details>

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/settheory).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
