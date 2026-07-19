package modelchecking

const initPredID = -1

// gpvwNode is a tableau node of the Gerth–Peled–Vardi–Wolper construction.
type gpvwNode struct {
	id       int
	incoming map[int]bool
	newF     []*LTL
	old      map[string]*LTL
	next     map[string]*LTL
}

// gpvw carries the mutable state of one tableau expansion.
type gpvw struct {
	nodes   []*gpvwNode
	counter int
}

func (g *gpvw) newID() int {
	id := g.counter
	g.counter++
	return id
}

func sameFormulaSet(a, b map[string]*LTL) bool {
	if len(a) != len(b) {
		return false
	}
	for k := range a {
		if _, ok := b[k]; !ok {
			return false
		}
	}
	return true
}

func copyIntSet(m map[int]bool) map[int]bool {
	c := make(map[int]bool, len(m))
	for k := range m {
		c[k] = true
	}
	return c
}

func copyFormulaSet(m map[string]*LTL) map[string]*LTL {
	c := make(map[string]*LTL, len(m))
	for k, v := range m {
		c[k] = v
	}
	return c
}

// negLiteralKey returns the ltlKey of the negation of a literal, and whether the
// argument was a literal at all.
func negLiteralKey(f *LTL) (string, bool) {
	switch f.Kind {
	case LTLAtomKind:
		return LTLNot(LTLVar(f.Atom)).ltlKey(), true
	case LTLNotKind:
		if f.L != nil && f.L.Kind == LTLAtomKind {
			return LTLVar(f.L.Atom).ltlKey(), true
		}
	case LTLTrueKind:
		return LTLBot().ltlKey(), true
	case LTLFalseKind:
		return LTLTop().ltlKey(), true
	}
	return "", false
}

// splitFormula returns the (new1, next1, new2) additions for a formula that
// forces a branch in the tableau: disjunction, until or release.
func splitFormula(eta *LTL) (new1 []*LTL, next1 []*LTL, new2 []*LTL) {
	switch eta.Kind {
	case LTLOrKind:
		return []*LTL{eta.L}, nil, []*LTL{eta.R}
	case LTLUntilKind:
		// a U b ≡ b ∨ (a ∧ X(a U b))
		return []*LTL{eta.L}, []*LTL{eta}, []*LTL{eta.R}
	case LTLReleaseKind:
		// a R b ≡ b ∧ (a ∨ X(a R b))
		return []*LTL{eta.R}, []*LTL{eta}, []*LTL{eta.L, eta.R}
	}
	return nil, nil, nil
}

func (g *gpvw) expand(node *gpvwNode) {
	if len(node.newF) == 0 {
		for _, r := range g.nodes {
			if sameFormulaSet(r.old, node.old) && sameFormulaSet(r.next, node.next) {
				for k := range node.incoming {
					r.incoming[k] = true
				}
				return
			}
		}
		g.nodes = append(g.nodes, node)
		succ := &gpvwNode{
			id:       g.newID(),
			incoming: map[int]bool{node.id: true},
			old:      map[string]*LTL{},
			next:     map[string]*LTL{},
		}
		for _, f := range node.next {
			succ.newF = append(succ.newF, f)
		}
		g.expand(succ)
		return
	}

	eta := node.newF[len(node.newF)-1]
	node.newF = node.newF[:len(node.newF)-1]
	key := eta.ltlKey()
	if _, ok := node.old[key]; ok {
		g.expand(node)
		return
	}

	switch {
	case eta.IsLiteral():
		if eta.Kind == LTLFalseKind {
			return // contradiction: discard node
		}
		if nk, ok := negLiteralKey(eta); ok {
			if _, clash := node.old[nk]; clash {
				return // contradiction
			}
		}
		node.old[key] = eta
		g.expand(node)
	case eta.Kind == LTLAndKind:
		node.old[key] = eta
		for _, f := range []*LTL{eta.L, eta.R} {
			if _, ok := node.old[f.ltlKey()]; !ok {
				node.newF = append(node.newF, f)
			}
		}
		g.expand(node)
	case eta.Kind == LTLNextKind:
		node.old[key] = eta
		node.next[eta.L.ltlKey()] = eta.L
		g.expand(node)
	default:
		// disjunction, until or release: branch into two nodes
		new1, next1, new2 := splitFormula(eta)
		node1 := &gpvwNode{
			id:       g.newID(),
			incoming: copyIntSet(node.incoming),
			newF:     append([]*LTL(nil), node.newF...),
			old:      copyFormulaSet(node.old),
			next:     copyFormulaSet(node.next),
		}
		node1.old[key] = eta
		for _, f := range new1 {
			if _, ok := node1.old[f.ltlKey()]; !ok {
				node1.newF = append(node1.newF, f)
			}
		}
		for _, f := range next1 {
			node1.next[f.ltlKey()] = f
		}
		node2 := &gpvwNode{
			id:       g.newID(),
			incoming: copyIntSet(node.incoming),
			newF:     append([]*LTL(nil), node.newF...),
			old:      copyFormulaSet(node.old),
			next:     copyFormulaSet(node.next),
		}
		node2.old[key] = eta
		for _, f := range new2 {
			if _, ok := node2.old[f.ltlKey()]; !ok {
				node2.newF = append(node2.newF, f)
			}
		}
		g.expand(node1)
		g.expand(node2)
	}
}

// guardOfOld builds the edge guard from the literals recorded in an old set.
func guardOfOld(old map[string]*LTL) Guard {
	var pos, neg []string
	for _, f := range old {
		switch {
		case f.Kind == LTLAtomKind:
			pos = append(pos, f.Atom)
		case f.Kind == LTLNotKind && f.L != nil && f.L.Kind == LTLAtomKind:
			neg = append(neg, f.L.Atom)
		}
	}
	return NewGuard(pos, neg)
}

// collectUntils returns the distinct Until subformulas of an NNF formula, keyed
// by their canonical string, in deterministic order.
func collectUntils(f *LTL) []*LTL {
	seen := map[string]bool{}
	var out []*LTL
	var rec func(*LTL)
	rec = func(n *LTL) {
		if n == nil {
			return
		}
		rec(n.L)
		rec(n.R)
		if n.Kind == LTLUntilKind {
			k := n.ltlKey()
			if !seen[k] {
				seen[k] = true
				out = append(out, n)
			}
		}
	}
	rec(f)
	return out
}

// LTLToGenBuchi builds a generalized Büchi automaton whose language is exactly
// the set of infinite words satisfying the LTL formula f, using the on-the-fly
// tableau construction of Gerth, Peled, Vardi and Wolper. The formula is first
// put into negation normal form by [LTL.NNF]. Edges are source-labelled by the
// literals asserted in a node, and there is one generalized acceptance set per
// Until subformula guaranteeing that every eventuality is eventually fulfilled.
func LTLToGenBuchi(f *LTL) *GenBuchi {
	nf := f.NNF()
	g := &gpvw{}
	initNode := &gpvwNode{
		id:       g.newID(),
		incoming: map[int]bool{initPredID: true},
		newF:     []*LTL{nf},
		old:      map[string]*LTL{},
		next:     map[string]*LTL{},
	}
	g.expand(initNode)

	m := len(g.nodes)
	gb := NewGenBuchi(m)
	gb.SetPropositions(f.Atoms())
	index := make(map[int]int, m)
	for j, node := range g.nodes {
		index[node.id] = j
	}
	for j, q := range g.nodes {
		if q.incoming[initPredID] {
			gb.SetInitial(j)
		}
		for pid := range q.incoming {
			if pid == initPredID {
				continue
			}
			if pi, ok := index[pid]; ok {
				// source-labelled edge p -> q with the source's literals
				gb.AddEdge(pi, j, guardOfOld(g.nodes[pi].old))
			}
		}
	}

	untils := collectUntils(nf)
	for _, u := range untils {
		uk := u.ltlKey()
		var nuKey string
		if u.R != nil {
			nuKey = u.R.ltlKey()
		}
		set := NewStateSet(m)
		for j, q := range g.nodes {
			_, hasU := q.old[uk]
			_, hasNu := q.old[nuKey]
			if !hasU || hasNu {
				set.Add(j)
			}
		}
		gb.AddAcceptanceSet(set)
	}
	return gb
}

// LTLToBuchi builds an ordinary Büchi automaton equivalent to the LTL formula f
// by degeneralizing the output of [LTLToGenBuchi].
func LTLToBuchi(f *LTL) *Buchi {
	return Degeneralize(LTLToGenBuchi(f))
}
