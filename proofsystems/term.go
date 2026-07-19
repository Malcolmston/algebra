package proofsystems

import (
	"sort"
	"strings"
)

// TermKind classifies the three syntactic categories of a first-order term:
// logical variables, constant symbols, and compound function applications.
type TermKind int

const (
	// KindVar denotes a first-order (logical) variable such as x or y.
	KindVar TermKind = iota
	// KindConst denotes a constant symbol (a nullary function) such as a or 0.
	KindConst
	// KindFunc denotes a compound term f(t1,...,tn) with n >= 1 arguments.
	KindFunc
)

// String returns the human-readable name of a term kind.
func (k TermKind) String() string {
	switch k {
	case KindVar:
		return "Var"
	case KindConst:
		return "Const"
	case KindFunc:
		return "Func"
	default:
		return "Unknown"
	}
}

// Term is an immutable first-order term. Depending on Kind it is a variable
// (Name holds the variable name), a constant (Name holds the constant symbol),
// or a function application (Name holds the function symbol and Args its
// arguments). Terms are compared and hashed structurally.
type Term struct {
	Kind TermKind
	Name string
	Args []Term
}

// NewVar builds a first-order variable term with the given name.
func NewVar(name string) Term {
	return Term{Kind: KindVar, Name: name}
}

// NewConst builds a constant term (a nullary function symbol) with the given
// name.
func NewConst(name string) Term {
	return Term{Kind: KindConst, Name: name}
}

// NewFunc builds a compound term applying the function symbol name to the given
// argument terms. With zero arguments it degenerates to a constant.
func NewFunc(name string, args ...Term) Term {
	if len(args) == 0 {
		return NewConst(name)
	}
	cp := make([]Term, len(args))
	copy(cp, args)
	return Term{Kind: KindFunc, Name: name, Args: cp}
}

// IsVar reports whether the term is a variable.
func (t Term) IsVar() bool { return t.Kind == KindVar }

// IsConst reports whether the term is a constant symbol.
func (t Term) IsConst() bool { return t.Kind == KindConst }

// IsFunc reports whether the term is a compound function application.
func (t Term) IsFunc() bool { return t.Kind == KindFunc }

// Arity returns the number of argument sub-terms; it is zero for variables and
// constants.
func (t Term) Arity() int { return len(t.Args) }

// Equal reports whether two terms are structurally identical.
func (t Term) Equal(u Term) bool {
	if t.Kind != u.Kind || t.Name != u.Name || len(t.Args) != len(u.Args) {
		return false
	}
	for i := range t.Args {
		if !t.Args[i].Equal(u.Args[i]) {
			return false
		}
	}
	return true
}

// String renders a term in standard f(a,g(x)) notation.
func (t Term) String() string {
	if len(t.Args) == 0 {
		return t.Name
	}
	parts := make([]string, len(t.Args))
	for i, a := range t.Args {
		parts[i] = a.String()
	}
	return t.Name + "(" + strings.Join(parts, ",") + ")"
}

// IsGround reports whether the term contains no variables.
func (t Term) IsGround() bool {
	if t.Kind == KindVar {
		return false
	}
	for _, a := range t.Args {
		if !a.IsGround() {
			return false
		}
	}
	return true
}

// Occurs reports whether the variable named v occurs anywhere within the term.
// It underlies the occurs-check of Robinson unification.
func (t Term) Occurs(v string) bool {
	if t.Kind == KindVar {
		return t.Name == v
	}
	for _, a := range t.Args {
		if a.Occurs(v) {
			return true
		}
	}
	return false
}

// Vars returns the sorted, de-duplicated list of variable names occurring in
// the term.
func (t Term) Vars() []string {
	set := map[string]bool{}
	t.collectVars(set)
	return sortedKeys(set)
}

func (t Term) collectVars(set map[string]bool) {
	if t.Kind == KindVar {
		set[t.Name] = true
		return
	}
	for _, a := range t.Args {
		a.collectVars(set)
	}
}

// Size returns the total number of symbol occurrences (nodes) in the term.
func (t Term) Size() int {
	n := 1
	for _, a := range t.Args {
		n += a.Size()
	}
	return n
}

// Depth returns the nesting depth of the term: 0 for a variable or constant,
// and 1 plus the maximum argument depth for a function application.
func (t Term) Depth() int {
	if len(t.Args) == 0 {
		return 0
	}
	max := 0
	for _, a := range t.Args {
		if d := a.Depth(); d > max {
			max = d
		}
	}
	return max + 1
}

// Rename returns a copy of the term with every variable name x replaced by the
// name r(x) if r has an entry for it. Non-variable structure is preserved.
func (t Term) Rename(r map[string]string) Term {
	switch t.Kind {
	case KindVar:
		if nv, ok := r[t.Name]; ok {
			return NewVar(nv)
		}
		return t
	case KindConst:
		return t
	default:
		args := make([]Term, len(t.Args))
		for i, a := range t.Args {
			args[i] = a.Rename(r)
		}
		return Term{Kind: KindFunc, Name: t.Name, Args: args}
	}
}

// Subterms returns every subterm of t including t itself, in pre-order.
func (t Term) Subterms() []Term {
	out := []Term{t}
	for _, a := range t.Args {
		out = append(out, a.Subterms()...)
	}
	return out
}

// TermsEqual reports whether two slices of terms are pointwise structurally
// equal.
func TermsEqual(a, b []Term) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !a[i].Equal(b[i]) {
			return false
		}
	}
	return true
}

// CollectVars returns the sorted union of variable names occurring in all of
// the given terms.
func CollectVars(terms ...Term) []string {
	set := map[string]bool{}
	for _, t := range terms {
		t.collectVars(set)
	}
	return sortedKeys(set)
}

func sortedKeys(set map[string]bool) []string {
	out := make([]string, 0, len(set))
	for k := range set {
		out = append(out, k)
	}
	sort.Strings(out)
	return out
}
