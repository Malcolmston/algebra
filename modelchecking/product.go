package modelchecking

import (
	"fmt"
	"strings"
)

// ProductBuchiKripke builds the synchronous product of a Kripke structure k and
// a Büchi automaton b whose edges are source-labelled by literal guards. The
// product is itself a Büchi automaton over the trivial alphabet: a product state
// (s, q) is encoded as the index s*b.NumStates() + q. There is a product edge
// (s, q) -> (s', q') whenever s -> s' is a transition of k and q -> q' is an edge
// of b whose guard is satisfied by the label of s. A product state is initial
// when both components are initial and accepting when the automaton component is
// accepting. An accepting run of the product projects to an infinite path of k
// whose trace is accepted by b.
func ProductBuchiKripke(k *Kripke, b *Buchi) *Buchi {
	nk, nb := k.n, b.n
	prod := NewBuchi(nk * nb)
	idx := func(s, q int) int { return s*nb + q }
	for s := 0; s < nk; s++ {
		label := k.labels[s]
		for q := 0; q < nb; q++ {
			for _, e := range b.edges[q] {
				if !e.Guard.Satisfies(label) {
					continue
				}
				for _, sp := range k.succ[s] {
					prod.AddEdge(idx(s, q), idx(sp, e.To), TrueGuard())
				}
			}
		}
	}
	for _, s := range k.InitialStates() {
		for _, q := range b.InitialStates() {
			prod.SetInitial(idx(s, q))
		}
	}
	for s := 0; s < nk; s++ {
		for q := 0; q < nb; q++ {
			if b.accept[q] {
				prod.SetAccepting(idx(s, q))
			}
		}
	}
	return prod
}

// projectLassoToKripke maps a product lasso back to the Kripke component by
// dividing each product state index by the number of Büchi states.
func projectLassoToKripke(l Lasso, nb int) Lasso {
	proj := func(states []int) []int {
		out := make([]int, len(states))
		for i, s := range states {
			out[i] = s / nb
		}
		return out
	}
	return Lasso{Prefix: proj(l.Prefix), Loop: proj(l.Loop)}
}

// Counterexample records a violating infinite path of a Kripke structure, found
// while model checking an LTL (or path-form CTL) property. The path is an
// ultimately periodic lasso over the states of Kripke.
type Counterexample struct {
	Kripke  *Kripke
	Lasso   Lasso
	Formula string
}

// StatePath returns one unrolling of the counterexample path (prefix followed by
// one copy of the loop).
func (c *Counterexample) StatePath() []int { return c.Lasso.States() }

// String renders the counterexample as a sequence of state names, marking the
// beginning of the periodic part with "(loop)".
func (c *Counterexample) String() string {
	if c == nil {
		return "<no counterexample>"
	}
	var b strings.Builder
	if c.Formula != "" {
		fmt.Fprintf(&b, "counterexample to %s:\n", c.Formula)
	}
	name := func(s int) string {
		if c.Kripke != nil {
			return c.Kripke.StateName(s)
		}
		return fmt.Sprintf("s%d", s)
	}
	for _, s := range c.Lasso.Prefix {
		fmt.Fprintf(&b, "  %s\n", name(s))
	}
	b.WriteString("  (loop)\n")
	for _, s := range c.Lasso.Loop {
		fmt.Fprintf(&b, "  %s\n", name(s))
	}
	return b.String()
}

// LTLModelCheck checks whether every initial path of k satisfies the LTL formula
// f. It builds a Büchi automaton for the negation of f, forms the product with
// k and tests emptiness. When the property holds it returns (true, nil, nil);
// otherwise it returns (false, cex, nil) with cex a violating lasso. The Kripke
// structure should be total; use [Kripke.MakeTotal] beforehand if it may have
// deadlocks.
func LTLModelCheck(k *Kripke, f *LTL) (bool, *Counterexample, error) {
	if f == nil {
		return false, nil, fmt.Errorf("modelchecking: nil LTL formula")
	}
	neg := LTLNot(f)
	b := LTLToBuchi(neg)
	prod := ProductBuchiKripke(k, b)
	lasso, ok := prod.AcceptingLasso()
	if !ok {
		return true, nil, nil
	}
	kl := projectLassoToKripke(lasso, b.NumStates())
	return false, &Counterexample{Kripke: k, Lasso: kl, Formula: f.String()}, nil
}

// LTLSatisfiable reports whether the LTL formula f is satisfiable, i.e. whether
// its Büchi automaton has a non-empty language. It also returns a satisfying
// ultimately periodic model as a lasso over the automaton states when one
// exists.
func LTLSatisfiable(f *LTL) (bool, Lasso) {
	b := LTLToBuchi(f)
	lasso, ok := b.AcceptingLasso()
	return ok, lasso
}

// LTLValid reports whether the LTL formula f holds on every infinite word, which
// is the case exactly when its negation is unsatisfiable.
func LTLValid(f *LTL) bool {
	sat, _ := LTLSatisfiable(LTLNot(f))
	return !sat
}

// LTLEquivalent reports whether two LTL formulas are satisfied by exactly the
// same infinite words, decided by checking validity of their biconditional.
func LTLEquivalent(f, g *LTL) bool {
	return LTLValid(LTLIff(f, g))
}
