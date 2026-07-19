package modelchecking

import "fmt"

// PreExists returns the existential predecessor image of target: the set of
// states having at least one successor in target. This is the semantics of the
// EX operator and the primitive step of every CTL fixpoint.
func PreExists(k *Kripke, target StateSet) StateSet {
	out := NewStateSet(k.n)
	for s := 0; s < k.n; s++ {
		for _, t := range k.succ[s] {
			if target.Contains(t) {
				out.Add(s)
				break
			}
		}
	}
	return out
}

// PreForall returns the universal predecessor image of target: the set of states
// all of whose successors lie in target. A state with no successors is included
// vacuously. This is the semantics of the AX operator.
func PreForall(k *Kripke, target StateSet) StateSet {
	return PreExists(k, target.Complement()).Complement()
}

// SatEX returns the set of states satisfying EX phi, given the set phi of states
// satisfying the argument.
func SatEX(k *Kripke, phi StateSet) StateSet { return PreExists(k, phi) }

// SatAX returns the set of states satisfying AX phi.
func SatAX(k *Kripke, phi StateSet) StateSet { return PreForall(k, phi) }

// SatEU returns the set of states satisfying E[a U b] by computing the least
// fixpoint Z = b ∪ (a ∩ EX Z).
func SatEU(k *Kripke, a, b StateSet) StateSet {
	z := b.Clone()
	for {
		next := b.Union(a.Intersect(PreExists(k, z)))
		if next.Equal(z) {
			return z
		}
		z = next
	}
}

// SatEF returns the set of states satisfying EF a ≡ E[true U a].
func SatEF(k *Kripke, a StateSet) StateSet {
	return SatEU(k, FullStateSet(k.n), a)
}

// SatEG returns the set of states satisfying EG a by computing the greatest
// fixpoint Z = a ∩ EX Z starting from a.
func SatEG(k *Kripke, a StateSet) StateSet {
	z := a.Clone()
	for {
		next := a.Intersect(PreExists(k, z))
		if next.Equal(z) {
			return z
		}
		z = next
	}
}

// SatAF returns the set of states satisfying AF a ≡ ¬EG ¬a.
func SatAF(k *Kripke, a StateSet) StateSet {
	return SatEG(k, a.Complement()).Complement()
}

// SatAG returns the set of states satisfying AG a ≡ ¬EF ¬a.
func SatAG(k *Kripke, a StateSet) StateSet {
	return SatEF(k, a.Complement()).Complement()
}

// SatAU returns the set of states satisfying A[a U b] by the least fixpoint
// Z = b ∪ (a ∩ AX Z).
func SatAU(k *Kripke, a, b StateSet) StateSet {
	z := b.Clone()
	for {
		next := b.Union(a.Intersect(PreForall(k, z)))
		if next.Equal(z) {
			return z
		}
		z = next
	}
}

// SatER returns the set of states satisfying E[a R b] by the greatest fixpoint
// Z = b ∩ (a ∪ EX Z).
func SatER(k *Kripke, a, b StateSet) StateSet {
	z := b.Clone()
	for {
		next := b.Intersect(a.Union(PreExists(k, z)))
		if next.Equal(z) {
			return z
		}
		z = next
	}
}

// SatAR returns the set of states satisfying A[a R b] by the greatest fixpoint
// Z = b ∩ (a ∪ AX Z).
func SatAR(k *Kripke, a, b StateSet) StateSet {
	z := b.Clone()
	for {
		next := b.Intersect(a.Union(PreForall(k, z)))
		if next.Equal(z) {
			return z
		}
		z = next
	}
}

// CTLCheck returns the set of states of k that satisfy the CTL formula f. It
// evaluates the formula bottom-up, computing a [StateSet] for every subformula
// via the fixpoint characterisations of the temporal operators.
func CTLCheck(k *Kripke, f *CTL) (StateSet, error) {
	if f == nil {
		return NewStateSet(k.n), fmt.Errorf("modelchecking: nil CTL formula")
	}
	switch f.Kind {
	case CTLTrueKind:
		return FullStateSet(k.n), nil
	case CTLFalseKind:
		return NewStateSet(k.n), nil
	case CTLAtomKind:
		return k.LabelStateSet(f.Atom), nil
	case CTLNotKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		return a.Complement(), nil
	case CTLAndKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		b, err := CTLCheck(k, f.R)
		if err != nil {
			return b, err
		}
		return a.Intersect(b), nil
	case CTLOrKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		b, err := CTLCheck(k, f.R)
		if err != nil {
			return b, err
		}
		return a.Union(b), nil
	case CTLImpliesKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		b, err := CTLCheck(k, f.R)
		if err != nil {
			return b, err
		}
		return a.Complement().Union(b), nil
	case CTLIffKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		b, err := CTLCheck(k, f.R)
		if err != nil {
			return b, err
		}
		both := a.Intersect(b)
		neither := a.Complement().Intersect(b.Complement())
		return both.Union(neither), nil
	case CTLEXKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		return SatEX(k, a), nil
	case CTLAXKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		return SatAX(k, a), nil
	case CTLEFKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		return SatEF(k, a), nil
	case CTLAFKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		return SatAF(k, a), nil
	case CTLEGKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		return SatEG(k, a), nil
	case CTLAGKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		return SatAG(k, a), nil
	case CTLEUKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		b, err := CTLCheck(k, f.R)
		if err != nil {
			return b, err
		}
		return SatEU(k, a, b), nil
	case CTLAUKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		b, err := CTLCheck(k, f.R)
		if err != nil {
			return b, err
		}
		return SatAU(k, a, b), nil
	case CTLERKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		b, err := CTLCheck(k, f.R)
		if err != nil {
			return b, err
		}
		return SatER(k, a, b), nil
	case CTLARKind:
		a, err := CTLCheck(k, f.L)
		if err != nil {
			return a, err
		}
		b, err := CTLCheck(k, f.R)
		if err != nil {
			return b, err
		}
		return SatAR(k, a, b), nil
	}
	return NewStateSet(k.n), fmt.Errorf("modelchecking: unknown CTL kind %d", f.Kind)
}

// CTLSatisfies reports whether state s of k satisfies f.
func CTLSatisfies(k *Kripke, s int, f *CTL) (bool, error) {
	set, err := CTLCheck(k, f)
	if err != nil {
		return false, err
	}
	return set.Contains(s), nil
}

// CTLModelCheck reports whether every initial state of k satisfies f. A Kripke
// structure with no initial states trivially satisfies every formula.
func CTLModelCheck(k *Kripke, f *CTL) (bool, error) {
	set, err := CTLCheck(k, f)
	if err != nil {
		return false, err
	}
	for _, s := range k.InitialStates() {
		if !set.Contains(s) {
			return false, nil
		}
	}
	return true, nil
}
