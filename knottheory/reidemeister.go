package knottheory

import "sort"

// This file implements combinatorial recognisers for the three Reidemeister
// moves acting on the signed Gauss code of a knot. The moves are the local
// diagram changes that generate ambient isotopy; a diagram is reducible by a
// move when the corresponding combinatorial pattern is present in its Gauss
// code.

// pair is an unordered pair of crossing labels used as a map key.
type pair struct{ a, b int }

func makePair(x, y int) pair {
	if x > y {
		x, y = y, x
	}
	return pair{x, y}
}

// ReidemeisterIKinks returns the crossing labels that form a Reidemeister I
// kink, that is a crossing whose two encounters are cyclically adjacent in the
// Gauss code. Such a crossing can be removed by a type I move without changing
// the knot type.
func (gc GaussCode) ReidemeisterIKinks() []int {
	n := len(gc.Entries)
	var out []int
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		if n >= 2 && gc.Entries[i].Crossing == gc.Entries[j].Crossing {
			out = append(out, gc.Entries[i].Crossing)
		}
	}
	sort.Ints(out)
	return out
}

// IsReidemeisterIReducible reports whether the diagram contains at least one
// type I kink.
func (gc GaussCode) IsReidemeisterIReducible() bool {
	return len(gc.ReidemeisterIKinks()) > 0
}

// RemoveReidemeisterIKink returns the Gauss code with the kink at the given
// crossing removed. It returns the original code and ok=false if that crossing
// is not a type I kink.
func (gc GaussCode) RemoveReidemeisterIKink(crossing int) (GaussCode, bool) {
	kinks := gc.ReidemeisterIKinks()
	found := false
	for _, k := range kinks {
		if k == crossing {
			found = true
			break
		}
	}
	if !found {
		return gc, false
	}
	out := make([]GaussEntry, 0, len(gc.Entries)-2)
	for _, e := range gc.Entries {
		if e.Crossing != crossing {
			out = append(out, e)
		}
	}
	return GaussCode{Entries: out}, true
}

// ReidemeisterIIReducible returns the pairs of crossing labels that form a
// reducible Reidemeister II bigon. Such a pair {a,b} occurs when, cyclically,
// the two crossings appear adjacent as an all-over arc (both encounters over)
// in one place and as an all-under arc (both encounters under) in another. Each
// returned slice has length two with the smaller label first.
func (gc GaussCode) ReidemeisterIIReducible() [][2]int {
	n := len(gc.Entries)
	if n < 4 {
		return nil
	}
	type flags struct{ over, under bool }
	adj := map[pair]*flags{}
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		e1, e2 := gc.Entries[i], gc.Entries[j]
		if e1.Crossing == e2.Crossing {
			continue
		}
		key := makePair(e1.Crossing, e2.Crossing)
		f := adj[key]
		if f == nil {
			f = &flags{}
			adj[key] = f
		}
		if e1.Over && e2.Over {
			f.over = true
		}
		if !e1.Over && !e2.Over {
			f.under = true
		}
	}
	var out [][2]int
	for k, f := range adj {
		if f.over && f.under {
			out = append(out, [2]int{k.a, k.b})
		}
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i][0] != out[j][0] {
			return out[i][0] < out[j][0]
		}
		return out[i][1] < out[j][1]
	})
	return out
}

// IsReidemeisterIIReducible reports whether a reducible type II bigon exists.
func (gc GaussCode) IsReidemeisterIIReducible() bool {
	return len(gc.ReidemeisterIIReducible()) > 0
}

// RemoveReidemeisterII returns the Gauss code with the two crossings of a
// reducible bigon removed. It returns ok=false when {a,b} is not a reducible
// type II pair.
func (gc GaussCode) RemoveReidemeisterII(a, b int) (GaussCode, bool) {
	ok := false
	for _, p := range gc.ReidemeisterIIReducible() {
		if (p[0] == a && p[1] == b) || (p[0] == b && p[1] == a) {
			ok = true
			break
		}
	}
	if !ok {
		return gc, false
	}
	out := make([]GaussEntry, 0, len(gc.Entries)-4)
	for _, e := range gc.Entries {
		if e.Crossing != a && e.Crossing != b {
			out = append(out, e)
		}
	}
	return GaussCode{Entries: out}, true
}

// ReidemeisterIIITriangles returns the triples of distinct crossing labels that
// are pairwise cyclically adjacent in the Gauss code. Such a triangle of chords
// is the combinatorial footprint required for a Reidemeister III move, so the
// returned triples are the candidate locations at which a type III move may be
// applied. Each triple is sorted in increasing order.
func (gc GaussCode) ReidemeisterIIITriangles() [][3]int {
	n := len(gc.Entries)
	if n < 6 {
		return nil
	}
	adjacent := map[pair]bool{}
	for i := 0; i < n; i++ {
		j := (i + 1) % n
		a, b := gc.Entries[i].Crossing, gc.Entries[j].Crossing
		if a != b {
			adjacent[makePair(a, b)] = true
		}
	}
	labels := gc.CrossingLabels()
	var out [][3]int
	for i := 0; i < len(labels); i++ {
		for j := i + 1; j < len(labels); j++ {
			for k := j + 1; k < len(labels); k++ {
				a, b, c := labels[i], labels[j], labels[k]
				if adjacent[makePair(a, b)] && adjacent[makePair(b, c)] && adjacent[makePair(a, c)] {
					out = append(out, [3]int{a, b, c})
				}
			}
		}
	}
	return out
}

// IsReidemeisterIIIApplicable reports whether the diagram contains a candidate
// type III triangle.
func (gc GaussCode) IsReidemeisterIIIApplicable() bool {
	return len(gc.ReidemeisterIIITriangles()) > 0
}

// IsReduced reports whether the diagram admits neither a type I kink nor a
// reducible type II bigon, a necessary (though not sufficient) condition for the
// diagram to be minimal.
func (gc GaussCode) IsReduced() bool {
	return !gc.IsReidemeisterIReducible() && !gc.IsReidemeisterIIReducible()
}
