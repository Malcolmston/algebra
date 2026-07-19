package simplicial

import (
	"math"
	"sort"
)

// FilteredSimplex is a simplex together with the scalar value at which it enters
// a filtration.
type FilteredSimplex struct {
	// Simplex is the simplex.
	Simplex Simplex
	// Value is the filtration value at which the simplex appears.
	Value float64
}

// Filtration is a simplicial complex whose simplices are equipped with
// non-decreasing entrance values: every face of a simplex must enter no later
// than the simplex itself. Growing the value threshold sweeps out a nested
// family of subcomplexes whose evolving homology is captured by persistent
// homology.
type Filtration struct {
	byKey map[string]FilteredSimplex
}

// NewFiltration returns an empty filtration.
func NewFiltration() *Filtration {
	return &Filtration{byKey: map[string]FilteredSimplex{}}
}

// Add inserts a single simplex at the given value. It does not automatically
// insert faces; use [Filtration.AddClosure] to keep the filtration valid.
// Re-adding a simplex keeps the smaller of the two values.
func (f *Filtration) Add(s Simplex, value float64) {
	if s.IsEmpty() {
		return
	}
	k := s.Key()
	if old, ok := f.byKey[k]; ok && old.Value <= value {
		return
	}
	f.byKey[k] = FilteredSimplex{Simplex: s, Value: value}
}

// AddClosure inserts s and every one of its faces at the given value (faces take
// the smaller of their existing and the new value), preserving validity.
func (f *Filtration) AddClosure(s Simplex, value float64) {
	for _, face := range s.Closure() {
		f.Add(face, value)
	}
}

// Len returns the number of simplices in the filtration.
func (f *Filtration) Len() int { return len(f.byKey) }

// Simplices returns the filtered simplices in filtration order: ascending by
// value, then by dimension, then by [CompareSimplices]. This order lists every
// face before the simplices it bounds.
func (f *Filtration) Simplices() []FilteredSimplex {
	out := make([]FilteredSimplex, 0, len(f.byKey))
	for _, fs := range f.byKey {
		out = append(out, fs)
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Value != out[j].Value {
			return out[i].Value < out[j].Value
		}
		return CompareSimplices(out[i].Simplex, out[j].Simplex) < 0
	})
	return out
}

// Value returns the filtration value of the simplex and whether it is present.
func (f *Filtration) Value(s Simplex) (float64, bool) {
	fs, ok := f.byKey[s.Key()]
	return fs.Value, ok
}

// ComplexAt returns the subcomplex of all simplices with value at most t.
func (f *Filtration) ComplexAt(t float64) *Complex {
	c := NewComplex()
	for _, fs := range f.byKey {
		if fs.Value <= t {
			c.simplices[fs.Simplex.Key()] = fs.Simplex
		}
	}
	return c
}

// RipsFiltration builds the Vietoris–Rips filtration of a point cloud up to the
// scale maxEps and dimension maxDim. Every simplex enters at its diameter — the
// largest pairwise distance among its vertices — with vertices entering at 0.
func RipsFiltration(pc *PointCloud, maxEps float64, maxDim int) *Filtration {
	dist := pc.DistanceMatrix()
	c := VietorisRipsFromDistances(dist, maxEps, maxDim)
	f := NewFiltration()
	for _, s := range c.Simplices() {
		vs := s.Vertices()
		var diam float64
		for a := 0; a < len(vs); a++ {
			for b := a + 1; b < len(vs); b++ {
				if d := dist[vs[a]][vs[b]]; d > diam {
					diam = d
				}
			}
		}
		f.Add(s, diam)
	}
	return f
}

// PersistencePair records the birth and death of a single homology class of a
// given dimension during a filtration. A death of +Inf marks an essential class
// that never dies.
type PersistencePair struct {
	// Dim is the homological dimension of the feature.
	Dim int
	// Birth is the value at which the class appears.
	Birth float64
	// Death is the value at which the class dies, or +Inf if it is essential.
	Death float64
	// BirthSimplex is the simplex that creates the class.
	BirthSimplex Simplex
	// DeathSimplex is the simplex that destroys the class; the zero Simplex for
	// essential classes.
	DeathSimplex Simplex
}

// Persistence returns Death − Birth, the lifetime of the feature (+Inf for
// essential classes).
func (p PersistencePair) Persistence() float64 { return p.Death - p.Birth }

// IsEssential reports whether the class never dies (infinite persistence).
func (p PersistencePair) IsEssential() bool { return math.IsInf(p.Death, 1) }

// Bar is a half-open interval [Birth, Death) in a persistence barcode.
type Bar struct {
	// Birth is the left endpoint.
	Birth float64
	// Death is the right endpoint, +Inf for an essential bar.
	Death float64
}

// Length returns the length of the bar (+Inf for an essential bar).
func (b Bar) Length() float64 { return b.Death - b.Birth }

// IsInfinite reports whether the bar extends to +Inf.
func (b Bar) IsInfinite() bool { return math.IsInf(b.Death, 1) }

// Persistence holds the result of a persistent-homology computation: the list of
// birth/death pairs across all dimensions.
type Persistence struct {
	pairs []PersistencePair
}

// symdiff returns the sorted symmetric difference of two sorted int slices,
// which is GF(2) column addition.
func symdiff(a, b []int) []int {
	out := make([]int, 0, len(a)+len(b))
	i, j := 0, 0
	for i < len(a) && j < len(b) {
		switch {
		case a[i] < b[j]:
			out = append(out, a[i])
			i++
		case a[i] > b[j]:
			out = append(out, b[j])
			j++
		default:
			i++
			j++
		}
	}
	out = append(out, a[i:]...)
	out = append(out, b[j:]...)
	return out
}

// PersistentHomology runs the standard matrix-reduction algorithm over GF(2) on
// the given filtration and returns the resulting persistence pairs.
func PersistentHomology(f *Filtration) *Persistence {
	order := f.Simplices()
	n := len(order)
	index := make(map[string]int, n)
	for i, fs := range order {
		index[fs.Simplex.Key()] = i
	}

	// boundary columns as sorted sets of row indices
	col := make([][]int, n)
	for j, fs := range order {
		var rows []int
		for _, face := range fs.Simplex.Faces() {
			if i, ok := index[face.Key()]; ok {
				rows = append(rows, i)
			}
		}
		sort.Ints(rows)
		col[j] = rows
	}

	low := func(c []int) int {
		if len(c) == 0 {
			return -1
		}
		return c[len(c)-1]
	}

	lowToCol := make(map[int]int)
	for j := 0; j < n; j++ {
		for len(col[j]) > 0 {
			l := low(col[j])
			if pj, ok := lowToCol[l]; ok && pj < j {
				col[j] = symdiff(col[j], col[pj])
			} else {
				break
			}
		}
		if len(col[j]) > 0 {
			lowToCol[low(col[j])] = j
		}
	}

	killed := make(map[int]bool) // birth indices that are paired
	var pairs []PersistencePair
	for j := 0; j < n; j++ {
		if len(col[j]) == 0 {
			continue
		}
		i := low(col[j])
		killed[i] = true
		birth := order[i].Value
		death := order[j].Value
		pairs = append(pairs, PersistencePair{
			Dim:          order[i].Simplex.Dim(),
			Birth:        birth,
			Death:        death,
			BirthSimplex: order[i].Simplex,
			DeathSimplex: order[j].Simplex,
		})
	}
	// essential classes: creators (empty reduced column) not paired
	for i := 0; i < n; i++ {
		if len(col[i]) == 0 && !killed[i] {
			pairs = append(pairs, PersistencePair{
				Dim:          order[i].Simplex.Dim(),
				Birth:        order[i].Value,
				Death:        math.Inf(1),
				BirthSimplex: order[i].Simplex,
			})
		}
	}
	sort.Slice(pairs, func(a, b int) bool {
		if pairs[a].Dim != pairs[b].Dim {
			return pairs[a].Dim < pairs[b].Dim
		}
		if pairs[a].Birth != pairs[b].Birth {
			return pairs[a].Birth < pairs[b].Birth
		}
		return pairs[a].Death < pairs[b].Death
	})
	return &Persistence{pairs: pairs}
}

// Pairs returns all persistence pairs across every dimension.
func (p *Persistence) Pairs() []PersistencePair {
	return append([]PersistencePair(nil), p.pairs...)
}

// PairsOfDim returns the persistence pairs of homological dimension k.
func (p *Persistence) PairsOfDim(k int) []PersistencePair {
	var out []PersistencePair
	for _, pr := range p.pairs {
		if pr.Dim == k {
			out = append(out, pr)
		}
	}
	return out
}

// NumPairs returns the total number of persistence pairs.
func (p *Persistence) NumPairs() int { return len(p.pairs) }

// Barcode returns the persistence barcode in dimension k: one [Bar] per feature.
func (p *Persistence) Barcode(k int) []Bar {
	var out []Bar
	for _, pr := range p.PairsOfDim(k) {
		out = append(out, Bar{Birth: pr.Birth, Death: pr.Death})
	}
	return out
}

// Diagram returns the persistence diagram in dimension k as a slice of
// [birth, death] points; essential classes have death +Inf.
func (p *Persistence) Diagram(k int) [][2]float64 {
	var out [][2]float64
	for _, pr := range p.PairsOfDim(k) {
		out = append(out, [2]float64{pr.Birth, pr.Death})
	}
	return out
}

// BettiAt returns the k-th Betti number of the subcomplex at threshold t,
// counting the features of dimension k that are alive at t (born at or before t
// and dying strictly after t).
func (p *Persistence) BettiAt(k int, t float64) int {
	count := 0
	for _, pr := range p.PairsOfDim(k) {
		if pr.Birth <= t && t < pr.Death {
			count++
		}
	}
	return count
}

// EssentialClasses returns the pairs that never die (infinite persistence),
// whose count in dimension k is the k-th Betti number of the final complex.
func (p *Persistence) EssentialClasses() []PersistencePair {
	var out []PersistencePair
	for _, pr := range p.pairs {
		if pr.IsEssential() {
			out = append(out, pr)
		}
	}
	return out
}
