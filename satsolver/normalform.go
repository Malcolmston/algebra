package satsolver

import "sort"

// literalClauses converts an NNF expression into a list of clauses, each clause
// a list of literal expressions understood disjunctively; the whole list is
// understood conjunctively (CNF form).
func cnfClauses(e Expr) [][]Expr {
	switch t := e.(type) {
	case And:
		return append(cnfClauses(t.X), cnfClauses(t.Y)...)
	case Or:
		left := cnfClauses(t.X)
		right := cnfClauses(t.Y)
		var out [][]Expr
		for _, a := range left {
			for _, b := range right {
				clause := make([]Expr, 0, len(a)+len(b))
				clause = append(clause, a...)
				clause = append(clause, b...)
				out = append(out, clause)
			}
		}
		return out
	default:
		return [][]Expr{{e}}
	}
}

// dnfTerms converts an NNF expression into a list of conjunctive terms
// understood disjunctively (DNF form).
func dnfTerms(e Expr) [][]Expr {
	switch t := e.(type) {
	case Or:
		return append(dnfTerms(t.X), dnfTerms(t.Y)...)
	case And:
		left := dnfTerms(t.X)
		right := dnfTerms(t.Y)
		var out [][]Expr
		for _, a := range left {
			for _, b := range right {
				term := make([]Expr, 0, len(a)+len(b))
				term = append(term, a...)
				term = append(term, b...)
				out = append(out, term)
			}
		}
		return out
	default:
		return [][]Expr{{e}}
	}
}

// ToCNFExpr converts an arbitrary Boolean expression into an equivalent
// expression in conjunctive normal form using De Morgan's laws and distribution
// of disjunction over conjunction. The result may be exponentially larger than
// the input; use [Tseitin] when only equisatisfiability is required.
func ToCNFExpr(e Expr) Expr {
	s := Simplify(e)
	if c, ok := s.(BoolConst); ok {
		return c
	}
	nnf := ToNNF(s)
	clauses := cnfClauses(nnf)
	ors := make([]Expr, len(clauses))
	for i, cl := range clauses {
		ors[i] = OrAll(cl...)
	}
	return AndAll(ors...)
}

// ToDNFExpr converts an arbitrary Boolean expression into an equivalent
// expression in disjunctive normal form.
func ToDNFExpr(e Expr) Expr {
	s := Simplify(e)
	if c, ok := s.(BoolConst); ok {
		return c
	}
	nnf := ToNNF(s)
	terms := dnfTerms(nnf)
	ands := make([]Expr, len(terms))
	for i, tm := range terms {
		ands[i] = AndAll(tm...)
	}
	return OrAll(ands...)
}

// IsCNFExpr reports whether e is syntactically in conjunctive normal form: a
// conjunction of disjunctions of literals.
func IsCNFExpr(e Expr) bool {
	var isLiteral func(Expr) bool
	isLiteral = func(x Expr) bool {
		switch t := x.(type) {
		case Variable, BoolConst:
			return true
		case Not:
			_, v := t.X.(Variable)
			return v
		}
		return false
	}
	var isClause func(Expr) bool
	isClause = func(x Expr) bool {
		if o, ok := x.(Or); ok {
			return isClause(o.X) && isClause(o.Y)
		}
		return isLiteral(x)
	}
	var walk func(Expr) bool
	walk = func(x Expr) bool {
		if a, ok := x.(And); ok {
			return walk(a.X) && walk(a.Y)
		}
		return isClause(x)
	}
	return walk(e)
}

// IsDNFExpr reports whether e is syntactically in disjunctive normal form: a
// disjunction of conjunctions of literals.
func IsDNFExpr(e Expr) bool {
	var isLiteral func(Expr) bool
	isLiteral = func(x Expr) bool {
		switch t := x.(type) {
		case Variable, BoolConst:
			return true
		case Not:
			_, v := t.X.(Variable)
			return v
		}
		return false
	}
	var isTerm func(Expr) bool
	isTerm = func(x Expr) bool {
		if a, ok := x.(And); ok {
			return isTerm(a.X) && isTerm(a.Y)
		}
		return isLiteral(x)
	}
	var walk func(Expr) bool
	walk = func(x Expr) bool {
		if o, ok := x.(Or); ok {
			return walk(o.X) && walk(o.Y)
		}
		return isTerm(x)
	}
	return walk(e)
}

// VarMap maps variable names to positive integer indices (1-based), the numeric
// encoding used by the clausal solver.
type VarMap struct {
	// Names lists the variable names in index order; Names[i-1] is the name of
	// variable i.
	Names []string
	index map[string]int
}

// NewVarMap builds a VarMap assigning indices 1..n to the given names in order.
func NewVarMap(names []string) *VarMap {
	m := &VarMap{index: map[string]int{}}
	for _, n := range names {
		m.Add(n)
	}
	return m
}

// Add returns the index of name, allocating a fresh one if it is new.
func (m *VarMap) Add(name string) int {
	if m.index == nil {
		m.index = map[string]int{}
	}
	if i, ok := m.index[name]; ok {
		return i
	}
	i := len(m.Names) + 1
	m.index[name] = i
	m.Names = append(m.Names, name)
	return i
}

// Index returns the index of name and whether it is present.
func (m *VarMap) Index(name string) (int, bool) {
	i, ok := m.index[name]
	return i, ok
}

// Name returns the variable name for index i (1-based), or the empty string if
// out of range.
func (m *VarMap) Name(i int) string {
	if i < 1 || i > len(m.Names) {
		return ""
	}
	return m.Names[i-1]
}

// Len returns the number of variables in the map.
func (m *VarMap) Len() int { return len(m.Names) }

// ToCNFFormula converts an expression to a clausal [CNF] together with the
// [VarMap] recording the name-to-index correspondence. The encoding is
// equivalence-preserving (via distribution), so it is exact rather than merely
// equisatisfiable.
func ToCNFFormula(e Expr) (CNF, *VarMap) {
	names := Vars(e)
	sort.Strings(names)
	vm := NewVarMap(names)
	s := Simplify(e)
	if c, ok := s.(BoolConst); ok {
		if bool(c) {
			return CNF{}, vm
		}
		// false: a single empty clause is unsatisfiable
		return CNF{Clauses: []Clause{{}}}, vm
	}
	nnf := ToNNF(s)
	clauses := cnfClauses(nnf)
	var out []Clause
	for _, cl := range clauses {
		clause, keep := literalsToClause(cl, vm)
		if keep {
			out = append(out, clause)
		}
	}
	return CNF{Clauses: out}, vm
}

// literalsToClause turns a slice of literal expressions into a numeric Clause.
// It returns keep=false when the clause is a tautology (contains a true
// constant or a literal and its complement).
func literalsToClause(lits []Expr, vm *VarMap) (Clause, bool) {
	var out Clause
	seen := map[Lit]bool{}
	for _, le := range lits {
		switch t := le.(type) {
		case BoolConst:
			if bool(t) {
				return nil, false // clause satisfied -> drop
			}
			// false literal contributes nothing
		case Variable:
			l := PosLit(vm.Add(string(t)))
			if seen[l.Negate()] {
				return nil, false
			}
			if !seen[l] {
				seen[l] = true
				out = append(out, l)
			}
		case Not:
			if v, ok := t.X.(Variable); ok {
				l := NegLit(vm.Add(string(v)))
				if seen[l.Negate()] {
					return nil, false
				}
				if !seen[l] {
					seen[l] = true
					out = append(out, l)
				}
			}
		}
	}
	return out, true
}

// ToDNFFormula converts an expression to a clausal [DNF] together with the
// [VarMap] recording the name-to-index correspondence.
func ToDNFFormula(e Expr) (DNF, *VarMap) {
	names := Vars(e)
	sort.Strings(names)
	vm := NewVarMap(names)
	s := Simplify(e)
	if c, ok := s.(BoolConst); ok {
		if bool(c) {
			return DNF{Terms: []Clause{{}}}, vm
		}
		return DNF{}, vm
	}
	nnf := ToNNF(s)
	terms := dnfTerms(nnf)
	var out []Clause
	for _, tm := range terms {
		term, keep := literalsToTerm(tm, vm)
		if keep {
			out = append(out, term)
		}
	}
	return DNF{Terms: out}, vm
}

// literalsToTerm turns literal expressions into a conjunctive numeric term. It
// returns keep=false when the term is contradictory (a false constant or a
// literal and its complement).
func literalsToTerm(lits []Expr, vm *VarMap) (Clause, bool) {
	var out Clause
	seen := map[Lit]bool{}
	for _, le := range lits {
		switch t := le.(type) {
		case BoolConst:
			if !bool(t) {
				return nil, false
			}
		case Variable:
			l := PosLit(vm.Add(string(t)))
			if seen[l.Negate()] {
				return nil, false
			}
			if !seen[l] {
				seen[l] = true
				out = append(out, l)
			}
		case Not:
			if v, ok := t.X.(Variable); ok {
				l := NegLit(vm.Add(string(v)))
				if seen[l.Negate()] {
					return nil, false
				}
				if !seen[l] {
					seen[l] = true
					out = append(out, l)
				}
			}
		}
	}
	return out, true
}

// CNFToExpr renders a clausal CNF as an [Expr] tree using the variable names in
// vm.
func CNFToExpr(f CNF, vm *VarMap) Expr {
	if len(f.Clauses) == 0 {
		return True
	}
	clauses := make([]Expr, len(f.Clauses))
	for i, c := range f.Clauses {
		if len(c) == 0 {
			clauses[i] = False
			continue
		}
		lits := make([]Expr, len(c))
		for j, l := range c {
			lits[j] = litToExpr(l, vm)
		}
		clauses[i] = OrAll(lits...)
	}
	return AndAll(clauses...)
}

// DNFToExpr renders a clausal DNF as an [Expr] tree using the names in vm.
func DNFToExpr(d DNF, vm *VarMap) Expr {
	if len(d.Terms) == 0 {
		return False
	}
	terms := make([]Expr, len(d.Terms))
	for i, t := range d.Terms {
		if len(t) == 0 {
			terms[i] = True
			continue
		}
		lits := make([]Expr, len(t))
		for j, l := range t {
			lits[j] = litToExpr(l, vm)
		}
		terms[i] = AndAll(lits...)
	}
	return OrAll(terms...)
}

func litToExpr(l Lit, vm *VarMap) Expr {
	name := vm.Name(l.Var())
	if name == "" {
		name = "x" + itoa(l.Var())
	}
	if l.IsNeg() {
		return Not{X: Variable(name)}
	}
	return Variable(name)
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}
