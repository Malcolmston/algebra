package automata

import "sort"

// Minimize returns the unique (up to renumbering) minimal DFA equivalent to the
// receiver. It first trims unreachable states, completes the transition
// function, then merges indistinguishable states using Hopcroft's algorithm.
func (d *DFA) Minimize() *DFA {
	return Hopcroft(d)
}

// Hopcroft returns the minimal DFA equivalent to d using Hopcroft's partition
// refinement algorithm, which runs in O(n·|Σ|·log n) time. Unreachable states
// are removed before minimisation.
func Hopcroft(d *DFA) *DFA {
	work := d.RemoveUnreachable().Complete()
	n := work.NumStates
	alphabet := work.Alphabet

	// Initial partition: accepting vs non-accepting.
	accepting := NewStateSet()
	nonAccepting := NewStateSet()
	for q := 0; q < n; q++ {
		if work.Accept[q] {
			accepting[q] = true
		} else {
			nonAccepting[q] = true
		}
	}
	var partition []StateSet
	if len(accepting) > 0 {
		partition = append(partition, accepting)
	}
	if len(nonAccepting) > 0 {
		partition = append(partition, nonAccepting)
	}

	// Precompute predecessors: pred[a][q] = states with an a-edge into q.
	pred := make(map[rune]map[int][]int)
	for _, a := range alphabet {
		pred[a] = make(map[int][]int)
	}
	for q := 0; q < n; q++ {
		for _, a := range alphabet {
			if to, ok := work.Transition(q, a); ok {
				pred[a][to] = append(pred[a][to], q)
			}
		}
	}

	// Worklist of (block, symbol) pairs to process, using block copies.
	worklist := make([]StateSet, 0, len(partition))
	for _, b := range partition {
		worklist = append(worklist, b.Clone())
	}

	for len(worklist) > 0 {
		splitter := worklist[len(worklist)-1]
		worklist = worklist[:len(worklist)-1]
		for _, a := range alphabet {
			// X = set of states with an a-transition into splitter.
			X := NewStateSet()
			for q := range splitter {
				for _, p := range pred[a][q] {
					X[p] = true
				}
			}
			if X.IsEmpty() {
				continue
			}
			var newPartition []StateSet
			for _, block := range partition {
				inter := block.Intersect(X)
				if inter.IsEmpty() || inter.Len() == block.Len() {
					newPartition = append(newPartition, block)
					continue
				}
				diff := NewStateSet()
				for q := range block {
					if !X[q] {
						diff[q] = true
					}
				}
				newPartition = append(newPartition, inter, diff)
				// Update worklist: replace block if present, else add smaller.
				replaced := false
				for i, w := range worklist {
					if w.Equal(block) {
						worklist[i] = inter.Clone()
						worklist = append(worklist, diff.Clone())
						replaced = true
						break
					}
				}
				if !replaced {
					if inter.Len() <= diff.Len() {
						worklist = append(worklist, inter.Clone())
					} else {
						worklist = append(worklist, diff.Clone())
					}
				}
			}
			partition = newPartition
		}
	}

	return buildQuotient(work, partition)
}

// Moore returns the minimal DFA equivalent to d using Moore's classic O(n²·|Σ|)
// partition-refinement algorithm. It provides an independent cross-check of the
// faster Hopcroft implementation.
func Moore(d *DFA) *DFA {
	work := d.RemoveUnreachable().Complete()
	n := work.NumStates
	alphabet := work.Alphabet

	// class[q] is the current block id of state q.
	class := make([]int, n)
	for q := 0; q < n; q++ {
		if work.Accept[q] {
			class[q] = 1
		} else {
			class[q] = 0
		}
	}
	numClasses := 2

	for {
		type sig struct {
			base int
			next []int
		}
		sigKey := make(map[string]int)
		newClass := make([]int, n)
		count := 0
		for q := 0; q < n; q++ {
			s := sig{base: class[q]}
			for _, a := range alphabet {
				to, _ := work.Transition(q, a)
				s.next = append(s.next, class[to])
			}
			key := encodeSig(s.base, s.next)
			id, ok := sigKey[key]
			if !ok {
				id = count
				sigKey[key] = id
				count++
			}
			newClass[q] = id
		}
		class = newClass
		if count == numClasses {
			break
		}
		numClasses = count
	}

	// Build partition blocks from class ids.
	blocks := make(map[int]StateSet)
	for q := 0; q < n; q++ {
		if blocks[class[q]] == nil {
			blocks[class[q]] = NewStateSet()
		}
		blocks[class[q]][q] = true
	}
	ids := make([]int, 0, len(blocks))
	for id := range blocks {
		ids = append(ids, id)
	}
	sort.Ints(ids)
	partition := make([]StateSet, 0, len(blocks))
	for _, id := range ids {
		partition = append(partition, blocks[id])
	}
	return buildQuotient(work, partition)
}

// MinimizeMoore is a method form of Moore's minimisation algorithm.
func (d *DFA) MinimizeMoore() *DFA {
	return Moore(d)
}

// EquivalenceClasses returns the blocks of the Myhill–Nerode equivalence over
// the reachable, completed DFA: each block is a set of mutually
// indistinguishable states (indices refer to the completed automaton).
func (d *DFA) EquivalenceClasses() [][]int {
	m := Hopcroft(d)
	// Recover blocks by grouping original states is non-trivial after
	// renumbering; instead recompute the partition on the completed DFA.
	work := d.RemoveUnreachable().Complete()
	part := indistinguishablePartition(work)
	// Sort deterministically.
	out := make([][]int, 0, len(part))
	for _, b := range part {
		out = append(out, b.Sorted())
	}
	sort.Slice(out, func(i, j int) bool { return out[i][0] < out[j][0] })
	_ = m
	return out
}

// indistinguishablePartition computes, for a completed DFA, the partition of
// states into Myhill–Nerode equivalence classes via the table-filling method.
func indistinguishablePartition(work *DFA) []StateSet {
	n := work.NumStates
	class := make([]int, n)
	for q := 0; q < n; q++ {
		if work.Accept[q] {
			class[q] = 1
		}
	}
	numClasses := 2
	for {
		sigKey := make(map[string]int)
		newClass := make([]int, n)
		count := 0
		for q := 0; q < n; q++ {
			next := make([]int, 0, len(work.Alphabet))
			for _, a := range work.Alphabet {
				to, _ := work.Transition(q, a)
				next = append(next, class[to])
			}
			key := encodeSig(class[q], next)
			id, ok := sigKey[key]
			if !ok {
				id = count
				sigKey[key] = id
				count++
			}
			newClass[q] = id
		}
		class = newClass
		if count == numClasses {
			break
		}
		numClasses = count
	}
	blocks := make(map[int]StateSet)
	for q := 0; q < n; q++ {
		if blocks[class[q]] == nil {
			blocks[class[q]] = NewStateSet()
		}
		blocks[class[q]][q] = true
	}
	out := make([]StateSet, 0, len(blocks))
	for _, b := range blocks {
		out = append(out, b)
	}
	return out
}

// buildQuotient collapses each block of the partition into a single state,
// preserving transitions and accepting status, and renumbers so the block
// containing the start state is the new start.
func buildQuotient(work *DFA, partition []StateSet) *DFA {
	blockOf := make([]int, work.NumStates)
	for i, block := range partition {
		for q := range block {
			blockOf[q] = i
		}
	}
	startBlock := blockOf[work.Start]
	out := NewDFA(len(partition), work.Alphabet, startBlock)
	for i, block := range partition {
		rep := block.Sorted()[0]
		for _, a := range work.Alphabet {
			if to, ok := work.Transition(rep, a); ok {
				out.SetTransition(i, a, blockOf[to])
			}
		}
		if work.Accept[rep] {
			out.Accept[i] = true
		}
	}
	return out.RemoveUnreachable()
}

// encodeSig builds a stable string key from a base class and successor classes.
func encodeSig(base int, next []int) string {
	buf := make([]byte, 0, 8*(len(next)+1))
	buf = appendInt(buf, base)
	buf = append(buf, '|')
	for _, v := range next {
		buf = appendInt(buf, v)
		buf = append(buf, ',')
	}
	return string(buf)
}

// appendInt appends the base-10 representation of v to buf.
func appendInt(buf []byte, v int) []byte {
	if v == 0 {
		return append(buf, '0')
	}
	neg := v < 0
	if neg {
		v = -v
	}
	var tmp [20]byte
	i := len(tmp)
	for v > 0 {
		i--
		tmp[i] = byte('0' + v%10)
		v /= 10
	}
	if neg {
		buf = append(buf, '-')
	}
	return append(buf, tmp[i:]...)
}
