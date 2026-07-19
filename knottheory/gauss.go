package knottheory

import (
	"fmt"
	"sort"
)

// GaussEntry is a single crossing encounter recorded while traversing a knot or
// link component. Crossing is the crossing label (a positive integer), Over is
// true when this pass is along the over-strand, and Sign is the sign of the
// crossing (+1 or -1) determined by the right-hand rule.
type GaussEntry struct {
	Crossing int
	Over     bool
	Sign     int
}

// GaussCode is the signed Gauss code of a single knot component: the cyclic
// sequence of crossing encounters met while traversing the knot once.
type GaussCode struct {
	Entries []GaussEntry
}

// NewGaussCode validates and returns a GaussCode. Every crossing label must
// occur exactly twice, once as an over-pass and once as an under-pass, with the
// same sign in both occurrences.
func NewGaussCode(entries []GaussEntry) (GaussCode, error) {
	gc := GaussCode{Entries: append([]GaussEntry(nil), entries...)}
	if err := gc.Validate(); err != nil {
		return GaussCode{}, err
	}
	return gc, nil
}

// Validate checks the structural consistency of the Gauss code and returns a
// descriptive error if it is malformed.
func (gc GaussCode) Validate() error {
	type rec struct {
		over, under int
		sign        int
		sawSign     bool
	}
	m := map[int]*rec{}
	for _, e := range gc.Entries {
		if e.Crossing <= 0 {
			return fmt.Errorf("knottheory: crossing labels must be positive, got %d", e.Crossing)
		}
		if e.Sign != 1 && e.Sign != -1 {
			return fmt.Errorf("knottheory: crossing %d has invalid sign %d", e.Crossing, e.Sign)
		}
		r := m[e.Crossing]
		if r == nil {
			r = &rec{}
			m[e.Crossing] = r
		}
		if r.sawSign && r.sign != e.Sign {
			return fmt.Errorf("knottheory: crossing %d has inconsistent signs", e.Crossing)
		}
		r.sign = e.Sign
		r.sawSign = true
		if e.Over {
			r.over++
		} else {
			r.under++
		}
	}
	for c, r := range m {
		if r.over != 1 || r.under != 1 {
			return fmt.Errorf("knottheory: crossing %d must appear once over and once under (got over=%d under=%d)", c, r.over, r.under)
		}
	}
	return nil
}

// Length returns the number of crossing encounters, which is twice the number
// of crossings.
func (gc GaussCode) Length() int { return len(gc.Entries) }

// CrossingLabels returns the sorted list of distinct crossing labels.
func (gc GaussCode) CrossingLabels() []int {
	seen := map[int]bool{}
	var out []int
	for _, e := range gc.Entries {
		if !seen[e.Crossing] {
			seen[e.Crossing] = true
			out = append(out, e.Crossing)
		}
	}
	sort.Ints(out)
	return out
}

// CrossingNumber returns the number of crossings of the diagram.
func (gc GaussCode) CrossingNumber() int { return len(gc.CrossingLabels()) }

// Writhe returns the writhe of the diagram, the sum of the crossing signs.
func (gc GaussCode) Writhe() int {
	seen := map[int]bool{}
	w := 0
	for _, e := range gc.Entries {
		if !seen[e.Crossing] {
			seen[e.Crossing] = true
			w += e.Sign
		}
	}
	return w
}

// Mirror returns the Gauss code of the mirror image, obtained by exchanging
// over and under at every crossing and negating every sign.
func (gc GaussCode) Mirror() GaussCode {
	out := make([]GaussEntry, len(gc.Entries))
	for i, e := range gc.Entries {
		out[i] = GaussEntry{Crossing: e.Crossing, Over: !e.Over, Sign: -e.Sign}
	}
	return GaussCode{Entries: out}
}

// ReverseOrientation returns the Gauss code obtained by traversing the knot in
// the opposite direction. The writhe is unchanged.
func (gc GaussCode) ReverseOrientation() GaussCode {
	n := len(gc.Entries)
	out := make([]GaussEntry, n)
	for i := 0; i < n; i++ {
		out[i] = gc.Entries[n-1-i]
	}
	return GaussCode{Entries: out}
}

// Diagram is a knot or link diagram given as one signed Gauss code per
// component.
type Diagram struct {
	Components []GaussCode
}

// NewDiagram validates and returns a Diagram from its component Gauss codes.
// Each crossing label must appear exactly twice across all components combined,
// once over and once under, with a consistent sign.
func NewDiagram(components ...GaussCode) (Diagram, error) {
	type rec struct {
		over, under int
		sign        int
		saw         bool
	}
	m := map[int]*rec{}
	for _, comp := range components {
		for _, e := range comp.Entries {
			if e.Crossing <= 0 {
				return Diagram{}, fmt.Errorf("knottheory: crossing labels must be positive, got %d", e.Crossing)
			}
			if e.Sign != 1 && e.Sign != -1 {
				return Diagram{}, fmt.Errorf("knottheory: crossing %d has invalid sign %d", e.Crossing, e.Sign)
			}
			r := m[e.Crossing]
			if r == nil {
				r = &rec{}
				m[e.Crossing] = r
			}
			if r.saw && r.sign != e.Sign {
				return Diagram{}, fmt.Errorf("knottheory: crossing %d has inconsistent signs", e.Crossing)
			}
			r.sign = e.Sign
			r.saw = true
			if e.Over {
				r.over++
			} else {
				r.under++
			}
		}
	}
	for c, r := range m {
		if r.over != 1 || r.under != 1 {
			return Diagram{}, fmt.Errorf("knottheory: crossing %d must appear once over and once under", c)
		}
	}
	cp := make([]GaussCode, len(components))
	copy(cp, components)
	return Diagram{Components: cp}, nil
}

// NumComponents returns the number of components of the diagram.
func (d Diagram) NumComponents() int { return len(d.Components) }

// CrossingNumber returns the number of crossings in the diagram.
func (d Diagram) CrossingNumber() int {
	seen := map[int]bool{}
	for _, comp := range d.Components {
		for _, e := range comp.Entries {
			seen[e.Crossing] = true
		}
	}
	return len(seen)
}

// Writhe returns the writhe of the whole diagram, the sum of all crossing
// signs.
func (d Diagram) Writhe() int {
	seen := map[int]bool{}
	w := 0
	for _, comp := range d.Components {
		for _, e := range comp.Entries {
			if !seen[e.Crossing] {
				seen[e.Crossing] = true
				w += e.Sign
			}
		}
	}
	return w
}

// crossingComponents returns, for each crossing label, the set of component
// indices it touches (a self-crossing touches a single component recorded once).
func (d Diagram) crossingComponents() map[int][]int {
	res := map[int][]int{}
	for ci, comp := range d.Components {
		for _, e := range comp.Entries {
			res[e.Crossing] = append(res[e.Crossing], ci)
		}
	}
	return res
}

// crossingSign returns the sign of a crossing label anywhere in the diagram.
func (d Diagram) crossingSign(label int) int {
	for _, comp := range d.Components {
		for _, e := range comp.Entries {
			if e.Crossing == label {
				return e.Sign
			}
		}
	}
	return 0
}

// LinkingNumber returns the Gauss linking number of components i and j, defined
// as half the sum of the signs of the crossings between the two components. It
// returns an error if the indices are out of range or equal.
func (d Diagram) LinkingNumber(i, j int) (int, error) {
	if i == j {
		return 0, fmt.Errorf("knottheory: linking number requires two distinct components")
	}
	if i < 0 || j < 0 || i >= len(d.Components) || j >= len(d.Components) {
		return 0, fmt.Errorf("knottheory: component index out of range")
	}
	cc := d.crossingComponents()
	sum := 0
	for label, comps := range cc {
		if len(comps) != 2 {
			continue
		}
		a, b := comps[0], comps[1]
		if (a == i && b == j) || (a == j && b == i) {
			sum += d.crossingSign(label)
		}
	}
	return sum / 2, nil
}

// TotalLinkingNumber returns the sum of the pairwise linking numbers over all
// unordered pairs of distinct components.
func (d Diagram) TotalLinkingNumber() int {
	total := 0
	for i := 0; i < len(d.Components); i++ {
		for j := i + 1; j < len(d.Components); j++ {
			lk, _ := d.LinkingNumber(i, j)
			total += lk
		}
	}
	return total
}

// SelfWrithe returns the writhe contributed only by self-crossings of the given
// component (crossings both of whose strands lie on that component).
func (d Diagram) SelfWrithe(comp int) int {
	if comp < 0 || comp >= len(d.Components) {
		return 0
	}
	cc := d.crossingComponents()
	w := 0
	for label, comps := range cc {
		if len(comps) == 2 && comps[0] == comp && comps[1] == comp {
			w += d.crossingSign(label)
		}
	}
	return w
}
