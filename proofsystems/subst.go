package proofsystems

import (
	"sort"
	"strings"
)

// Substitution maps variable names to replacement terms. It is the fundamental
// data structure of unification and resolution: applying a substitution
// simultaneously replaces every occurrence of each mapped variable by its
// image term. The zero value (a nil map wrapped by NewSubstitution) is the
// identity substitution.
type Substitution struct {
	m map[string]Term
}

// NewSubstitution returns the empty (identity) substitution.
func NewSubstitution() Substitution {
	return Substitution{m: map[string]Term{}}
}

// SingletonSubstitution returns the substitution {v -> t}.
func SingletonSubstitution(v string, t Term) Substitution {
	return Substitution{m: map[string]Term{v: t}}
}

// IsEmpty reports whether the substitution has no bindings (is the identity).
func (s Substitution) IsEmpty() bool { return len(s.m) == 0 }

// Len returns the number of variable bindings.
func (s Substitution) Len() int { return len(s.m) }

// Get returns the image of variable v and whether it is bound.
func (s Substitution) Get(v string) (Term, bool) {
	t, ok := s.m[v]
	return t, ok
}

// Bind returns a copy of the substitution with v additionally mapped to t. It
// does not mutate the receiver.
func (s Substitution) Bind(v string, t Term) Substitution {
	m := make(map[string]Term, len(s.m)+1)
	for k, val := range s.m {
		m[k] = val
	}
	m[v] = t
	return Substitution{m: m}
}

// Domain returns the sorted list of variable names bound by the substitution.
func (s Substitution) Domain() []string {
	out := make([]string, 0, len(s.m))
	for k := range s.m {
		out = append(out, k)
	}
	sort.Strings(out)
	return out
}

// ApplyTerm applies the substitution to a term, recursively replacing bound
// variables by their images. Images are themselves substituted so the result
// is idempotent when the substitution is already resolved.
func (s Substitution) ApplyTerm(t Term) Term {
	switch t.Kind {
	case KindVar:
		if img, ok := s.m[t.Name]; ok {
			return img
		}
		return t
	case KindConst:
		return t
	default:
		args := make([]Term, len(t.Args))
		for i, a := range t.Args {
			args[i] = s.ApplyTerm(a)
		}
		return Term{Kind: KindFunc, Name: t.Name, Args: args}
	}
}

// ApplyTerms applies the substitution to each term in a slice.
func (s Substitution) ApplyTerms(ts []Term) []Term {
	out := make([]Term, len(ts))
	for i, t := range ts {
		out[i] = s.ApplyTerm(t)
	}
	return out
}

// Compose returns the composition s then o (written o∘s): applying the result
// to a term is the same as applying s and then o. Bindings of s have o applied
// to their images, and any binding of o whose variable is not in the domain of
// s is appended.
func (s Substitution) Compose(o Substitution) Substitution {
	m := make(map[string]Term, len(s.m)+len(o.m))
	for v, t := range s.m {
		m[v] = o.ApplyTerm(t)
	}
	for v, t := range o.m {
		if _, ok := s.m[v]; !ok {
			m[v] = t
		}
	}
	// Drop trivial x -> x bindings that composition may introduce.
	for v, t := range m {
		if t.IsVar() && t.Name == v {
			delete(m, v)
		}
	}
	return Substitution{m: m}
}

// Restrict returns the substitution limited to the given set of variable
// names.
func (s Substitution) Restrict(vars []string) Substitution {
	keep := map[string]bool{}
	for _, v := range vars {
		keep[v] = true
	}
	m := map[string]Term{}
	for v, t := range s.m {
		if keep[v] {
			m[v] = t
		}
	}
	return Substitution{m: m}
}

// Equal reports whether two substitutions have identical bindings.
func (s Substitution) Equal(o Substitution) bool {
	if len(s.m) != len(o.m) {
		return false
	}
	for v, t := range s.m {
		u, ok := o.m[v]
		if !ok || !t.Equal(u) {
			return false
		}
	}
	return true
}

// String renders a substitution as {x -> t, y -> u} with the domain sorted.
func (s Substitution) String() string {
	if len(s.m) == 0 {
		return "{}"
	}
	parts := make([]string, 0, len(s.m))
	for _, v := range s.Domain() {
		parts = append(parts, v+" -> "+s.m[v].String())
	}
	return "{" + strings.Join(parts, ", ") + "}"
}

// ApplyLiteral applies the substitution to every argument term of a first-order
// literal.
func (s Substitution) ApplyLiteral(l FOLiteral) FOLiteral {
	return FOLiteral{Neg: l.Neg, Pred: l.Pred, Args: s.ApplyTerms(l.Args)}
}

// ApplyClause applies the substitution to every literal of a first-order
// clause.
func (s Substitution) ApplyClause(c FOClause) FOClause {
	lits := make([]FOLiteral, len(c.Lits))
	for i, l := range c.Lits {
		lits[i] = s.ApplyLiteral(l)
	}
	return FOClause{Lits: lits}
}
