package satsolver

// Node identifies a vertex of a reduced ordered binary decision diagram managed
// by a [BDD]. The constants [BDDFalse] and [BDDTrue] are the two terminals.
type Node int

const (
	// BDDFalse is the false terminal node.
	BDDFalse Node = 0
	// BDDTrue is the true terminal node.
	BDDTrue Node = 1
)

const bddTerminalLevel = 1 << 30

type bddNode struct {
	level int // index into the variable ordering; terminals use bddTerminalLevel
	lo    Node
	hi    Node
}

// BDD is a reduced ordered binary decision diagram manager. It fixes a variable
// ordering at construction and maintains a shared, canonical set of nodes so
// that logically equal functions are represented by the identical [Node].
type BDD struct {
	order    []string
	varIndex map[string]int
	nodes    []bddNode
	unique   map[[3]int]Node
	iteCache map[[3]Node]Node
}

// NewBDD returns a manager whose variables are tested in the given order, most
// significant (closest to the root) first.
func NewBDD(order []string) *BDD {
	b := &BDD{
		order:    append([]string(nil), order...),
		varIndex: make(map[string]int, len(order)),
		unique:   map[[3]int]Node{},
		iteCache: map[[3]Node]Node{},
	}
	for i, name := range order {
		b.varIndex[name] = i
	}
	// Terminal nodes 0 (false) and 1 (true).
	b.nodes = []bddNode{
		{level: bddTerminalLevel},
		{level: bddTerminalLevel},
	}
	return b
}

// Order returns the variable ordering used by the manager.
func (b *BDD) Order() []string { return append([]string(nil), b.order...) }

// NumNodes returns the total number of distinct nodes allocated, including the
// two terminals.
func (b *BDD) NumNodes() int { return len(b.nodes) }

func (b *BDD) level(n Node) int { return b.nodes[n].level }

// mk returns the reduced node testing the variable at the given level with the
// given low and high children, applying the ROBDD reduction rules.
func (b *BDD) mk(level int, lo, hi Node) Node {
	if lo == hi {
		return lo
	}
	key := [3]int{level, int(lo), int(hi)}
	if n, ok := b.unique[key]; ok {
		return n
	}
	n := Node(len(b.nodes))
	b.nodes = append(b.nodes, bddNode{level: level, lo: lo, hi: hi})
	b.unique[key] = n
	return n
}

// Var returns the node representing the single variable name. The variable must
// be part of the manager's ordering.
func (b *BDD) Var(name string) Node {
	i, ok := b.varIndex[name]
	if !ok {
		i = len(b.order)
		b.order = append(b.order, name)
		b.varIndex[name] = i
	}
	return b.mk(i, BDDFalse, BDDTrue)
}

// cofactor returns the low or high child of n if n tests the variable at the
// given level, or n itself otherwise.
func (b *BDD) cofactor(n Node, level int, high bool) Node {
	if b.level(n) != level {
		return n
	}
	if high {
		return b.nodes[n].hi
	}
	return b.nodes[n].lo
}

// Ite computes the if-then-else operator ite(f, g, h) = (f AND g) OR (NOT f AND
// h), the universal ternary connective from which all binary Boolean operators
// are derived.
func (b *BDD) Ite(f, g, h Node) Node {
	// Terminal simplifications.
	if f == BDDTrue {
		return g
	}
	if f == BDDFalse {
		return h
	}
	if g == h {
		return g
	}
	if g == BDDTrue && h == BDDFalse {
		return f
	}
	key := [3]Node{f, g, h}
	if r, ok := b.iteCache[key]; ok {
		return r
	}
	lv := b.level(f)
	if l := b.level(g); l < lv {
		lv = l
	}
	if l := b.level(h); l < lv {
		lv = l
	}
	lo := b.Ite(b.cofactor(f, lv, false), b.cofactor(g, lv, false), b.cofactor(h, lv, false))
	hi := b.Ite(b.cofactor(f, lv, true), b.cofactor(g, lv, true), b.cofactor(h, lv, true))
	r := b.mk(lv, lo, hi)
	b.iteCache[key] = r
	return r
}

// Not returns the negation of n.
func (b *BDD) Not(n Node) Node { return b.Ite(n, BDDFalse, BDDTrue) }

// And returns the conjunction of f and g.
func (b *BDD) And(f, g Node) Node { return b.Ite(f, g, BDDFalse) }

// Or returns the disjunction of f and g.
func (b *BDD) Or(f, g Node) Node { return b.Ite(f, BDDTrue, g) }

// Xor returns the exclusive-or of f and g.
func (b *BDD) Xor(f, g Node) Node { return b.Ite(f, b.Not(g), g) }

// Nand returns the negated conjunction of f and g.
func (b *BDD) Nand(f, g Node) Node { return b.Not(b.And(f, g)) }

// Nor returns the negated disjunction of f and g.
func (b *BDD) Nor(f, g Node) Node { return b.Not(b.Or(f, g)) }

// Iff returns the biconditional (equivalence) of f and g.
func (b *BDD) Iff(f, g Node) Node { return b.Ite(f, g, b.Not(g)) }

// Implies returns the conditional f -> g.
func (b *BDD) Implies(f, g Node) Node { return b.Ite(f, g, BDDTrue) }

// Apply combines f and g with the named binary operator ("and", "or", "xor",
// "nand", "nor", "iff", "implies"). Unknown names default to conjunction.
func (b *BDD) Apply(op string, f, g Node) Node {
	switch op {
	case "or":
		return b.Or(f, g)
	case "xor":
		return b.Xor(f, g)
	case "nand":
		return b.Nand(f, g)
	case "nor":
		return b.Nor(f, g)
	case "iff":
		return b.Iff(f, g)
	case "implies":
		return b.Implies(f, g)
	default:
		return b.And(f, g)
	}
}

// Restrict returns the cofactor of n obtained by fixing the variable name to
// val (Shannon restriction).
func (b *BDD) Restrict(n Node, name string, val bool) Node {
	i, ok := b.varIndex[name]
	if !ok {
		return n
	}
	var rec func(Node) Node
	rec = func(m Node) Node {
		if m == BDDFalse || m == BDDTrue {
			return m
		}
		lvl := b.level(m)
		if lvl > i {
			return m
		}
		lo := b.nodes[m].lo
		hi := b.nodes[m].hi
		if lvl == i {
			if val {
				return hi
			}
			return lo
		}
		return b.mk(lvl, rec(lo), rec(hi))
	}
	return rec(n)
}

// IsTautology reports whether n is the true terminal, i.e. the function is
// identically true.
func (b *BDD) IsTautology(n Node) bool { return n == BDDTrue }

// IsContradiction reports whether n is the false terminal.
func (b *BDD) IsContradiction(n Node) bool { return n == BDDFalse }

// IsSatisfiable reports whether n is not the false terminal.
func (b *BDD) IsSatisfiable(n Node) bool { return n != BDDFalse }

// SatCount returns the number of assignments over all ordering variables under
// which the function represented by n is true.
func (b *BDD) SatCount(n Node) int {
	N := len(b.order)
	memo := map[Node]int{}
	var paths func(Node) int
	paths = func(m Node) int {
		if m == BDDFalse {
			return 0
		}
		if m == BDDTrue {
			return 1
		}
		if v, ok := memo[m]; ok {
			return v
		}
		lvl := b.level(m)
		lo := b.nodes[m].lo
		hi := b.nodes[m].hi
		gapLo := b.levelOrN(lo, N) - lvl - 1
		gapHi := b.levelOrN(hi, N) - lvl - 1
		res := (1<<gapLo)*paths(lo) + (1<<gapHi)*paths(hi)
		memo[m] = res
		return res
	}
	if n == BDDFalse {
		return 0
	}
	if n == BDDTrue {
		return 1 << N
	}
	top := b.level(n)
	return (1 << top) * paths(n)
}

func (b *BDD) levelOrN(n Node, N int) int {
	if n == BDDFalse || n == BDDTrue {
		return N
	}
	return b.level(n)
}

// AnySat returns one satisfying assignment of the function n as a name->value
// map and true, or nil and false when n is unsatisfiable. Variables not tested
// on the chosen path are set to false.
func (b *BDD) AnySat(n Node) (map[string]bool, bool) {
	if n == BDDFalse {
		return nil, false
	}
	env := map[string]bool{}
	for _, name := range b.order {
		env[name] = false
	}
	m := n
	for m != BDDTrue {
		lvl := b.level(m)
		name := b.order[lvl]
		if b.nodes[m].hi != BDDFalse {
			env[name] = true
			m = b.nodes[m].hi
		} else {
			env[name] = false
			m = b.nodes[m].lo
		}
	}
	return env, true
}

// Eval evaluates the function n under the given environment.
func (b *BDD) Eval(n Node, env map[string]bool) bool {
	m := n
	for m != BDDFalse && m != BDDTrue {
		name := b.order[b.level(m)]
		if env[name] {
			m = b.nodes[m].hi
		} else {
			m = b.nodes[m].lo
		}
	}
	return m == BDDTrue
}

// FromExpr builds the BDD node for the expression e, testing its variables in
// the manager's ordering. Variables absent from the ordering are appended.
func (b *BDD) FromExpr(e Expr) Node {
	switch t := e.(type) {
	case BoolConst:
		if bool(t) {
			return BDDTrue
		}
		return BDDFalse
	case Variable:
		return b.Var(string(t))
	case Not:
		return b.Not(b.FromExpr(t.X))
	case And:
		return b.And(b.FromExpr(t.X), b.FromExpr(t.Y))
	case Or:
		return b.Or(b.FromExpr(t.X), b.FromExpr(t.Y))
	case Xor:
		return b.Xor(b.FromExpr(t.X), b.FromExpr(t.Y))
	case Implies:
		return b.Implies(b.FromExpr(t.X), b.FromExpr(t.Y))
	case Iff:
		return b.Iff(b.FromExpr(t.X), b.FromExpr(t.Y))
	case Nand:
		return b.Nand(b.FromExpr(t.X), b.FromExpr(t.Y))
	case Nor:
		return b.Nor(b.FromExpr(t.X), b.FromExpr(t.Y))
	case Xnor:
		return b.Iff(b.FromExpr(t.X), b.FromExpr(t.Y))
	}
	return BDDFalse
}

// NewBDDFromExpr is a convenience that builds a manager whose ordering is the
// sorted variables of e and returns both the manager and the root node.
func NewBDDFromExpr(e Expr) (*BDD, Node) {
	b := NewBDD(Vars(e))
	return b, b.FromExpr(e)
}

// Equal reports whether two nodes denote the same Boolean function. Because the
// diagram is reduced and ordered, this is exactly node identity.
func (b *BDD) Equal(f, g Node) bool { return f == g }

// ToExpr reconstructs a Boolean [Expr] from the node n by Shannon expansion,
// using the manager's variable names.
func (b *BDD) ToExpr(n Node) Expr {
	switch n {
	case BDDFalse:
		return False
	case BDDTrue:
		return True
	}
	name := b.order[b.level(n)]
	lo := b.ToExpr(b.nodes[n].lo)
	hi := b.ToExpr(b.nodes[n].hi)
	// n = (v & hi) | (~v & lo)
	return Simplify(Or{
		X: And{X: Variable(name), Y: hi},
		Y: And{X: Not{X: Variable(name)}, Y: lo},
	})
}
