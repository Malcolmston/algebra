package satsolver

// TseitinEncoder builds an equisatisfiable CNF encoding of a Boolean expression
// by introducing a fresh auxiliary variable for every compound subexpression.
// The resulting formula grows only linearly in the size of the input, unlike
// the potentially exponential distribution performed by [ToCNFFormula].
type TseitinEncoder struct {
	vm      *VarMap
	clauses []Clause
	cache   map[string]Lit
}

// NewTseitinEncoder returns an empty encoder ready to encode expressions.
func NewTseitinEncoder() *TseitinEncoder {
	return &TseitinEncoder{
		vm:    NewVarMap(nil),
		cache: map[string]Lit{},
	}
}

// VarMap returns the variable map that records both original and auxiliary
// variable names introduced during encoding.
func (t *TseitinEncoder) VarMap() *VarMap { return t.vm }

// Clauses returns the clauses accumulated so far.
func (t *TseitinEncoder) Clauses() []Clause { return t.clauses }

func (t *TseitinEncoder) fresh() Lit {
	n := t.vm.Add("__aux" + itoa(t.vm.Len()+1))
	return PosLit(n)
}

func (t *TseitinEncoder) add(c ...Lit) {
	t.clauses = append(t.clauses, NewClause(c...))
}

// Encode registers the clauses that define e and returns a literal that is true
// exactly when e is true. Repeated identical subexpressions share their
// auxiliary variable.
func (t *TseitinEncoder) Encode(e Expr) Lit {
	key := e.String()
	if l, ok := t.cache[key]; ok {
		return l
	}
	l := t.encode(e)
	t.cache[key] = l
	return l
}

func (t *TseitinEncoder) encode(e Expr) Lit {
	switch n := e.(type) {
	case BoolConst:
		g := t.fresh()
		if bool(n) {
			t.add(g) // g must be true
		} else {
			t.add(g.Negate())
		}
		return g
	case Variable:
		return PosLit(t.vm.Add(string(n)))
	case Not:
		return t.Encode(n.X).Negate()
	case And:
		a := t.Encode(n.X)
		b := t.Encode(n.Y)
		g := t.fresh()
		// g <-> (a & b)
		t.add(a, g.Negate())
		t.add(b, g.Negate())
		t.add(a.Negate(), b.Negate(), g)
		return g
	case Or:
		a := t.Encode(n.X)
		b := t.Encode(n.Y)
		g := t.fresh()
		// g <-> (a | b)
		t.add(a.Negate(), g)
		t.add(b.Negate(), g)
		t.add(a, b, g.Negate())
		return g
	case Xor:
		a := t.Encode(n.X)
		b := t.Encode(n.Y)
		g := t.fresh()
		// g <-> (a xor b)
		t.add(a.Negate(), b.Negate(), g.Negate())
		t.add(a, b, g.Negate())
		t.add(a, b.Negate(), g)
		t.add(a.Negate(), b, g)
		return g
	case Iff, Xnor:
		var x, y Expr
		if v, ok := n.(Iff); ok {
			x, y = v.X, v.Y
		} else {
			v := n.(Xnor)
			x, y = v.X, v.Y
		}
		a := t.Encode(x)
		b := t.Encode(y)
		g := t.fresh()
		// g <-> (a == b)
		t.add(a.Negate(), b.Negate(), g)
		t.add(a, b, g)
		t.add(a, b.Negate(), g.Negate())
		t.add(a.Negate(), b, g.Negate())
		return g
	case Implies:
		a := t.Encode(n.X)
		b := t.Encode(n.Y)
		g := t.fresh()
		// g <-> (~a | b)
		t.add(a, g)
		t.add(b.Negate(), g)
		t.add(a.Negate(), b, g.Negate())
		return g
	case Nand:
		return t.Encode(Not{X: And{X: n.X, Y: n.Y}})
	case Nor:
		return t.Encode(Not{X: Or{X: n.X, Y: n.Y}})
	}
	return t.fresh()
}

// Assert records that the encoded expression must be true, by adding a unit
// clause on its top-level literal.
func (t *TseitinEncoder) Assert(e Expr) {
	l := t.Encode(e)
	t.add(l)
}

// CNF returns the accumulated clauses as a [CNF] formula.
func (t *TseitinEncoder) CNF() CNF {
	cs := make([]Clause, len(t.clauses))
	copy(cs, t.clauses)
	return CNF{Clauses: cs}
}

// Tseitin returns an equisatisfiable CNF encoding of e that asserts e is true,
// together with the variable map used. The original variables of e keep the
// lowest indices; auxiliary variables follow.
func Tseitin(e Expr) (CNF, *VarMap) {
	enc := NewTseitinEncoder()
	// Register the original variables first so they get the lowest indices.
	for _, name := range Vars(e) {
		enc.vm.Add(name)
	}
	enc.Assert(e)
	return enc.CNF(), enc.vm
}

// TseitinCNF returns just the equisatisfiable CNF from [Tseitin].
func TseitinCNF(e Expr) CNF {
	f, _ := Tseitin(e)
	return f
}
