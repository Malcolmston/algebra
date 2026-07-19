package proofsystems

import "fmt"

// TseitinResult is the outcome of a Tseitin transformation: the equisatisfiable
// clause set together with the name of the auxiliary variable representing the
// whole input formula. Asserting Top as a unit clause makes the CNF satisfiable
// exactly when the original formula is.
type TseitinResult struct {
	CNF PCNF
	Top string
}

// Tseitin performs the Tseitin transformation of a propositional formula,
// introducing a fresh auxiliary variable for each compound subformula and
// emitting the defining clauses that make the auxiliary equivalent to its
// subformula. The result is equisatisfiable with the input and grows linearly
// in the size of the formula. Auxiliary variables are named with the given
// prefix followed by an index.
func Tseitin(f Formula, prefix string) TseitinResult {
	t := &tseitinBuilder{prefix: prefix, memo: map[string]string{}}
	top := t.encode(ToNNF(f))
	return TseitinResult{CNF: PCNF{Clauses: t.clauses}, Top: top}
}

// TseitinCNF returns the equisatisfiable clause set of f with the top auxiliary
// variable asserted as a unit clause, so the returned CNF is satisfiable if and
// only if f is.
func TseitinCNF(f Formula, prefix string) PCNF {
	r := Tseitin(f, prefix)
	clauses := append([]PClause{}, r.CNF.Clauses...)
	clauses = append(clauses, NewPClause(PosPLit(r.Top)))
	return PCNF{Clauses: clauses}
}

type tseitinBuilder struct {
	prefix  string
	counter int
	clauses []PClause
	memo    map[string]string
}

func (t *tseitinBuilder) fresh() string {
	name := fmt.Sprintf("%s%d", t.prefix, t.counter)
	t.counter++
	return name
}

// encode returns the name of a literal-bearing variable equivalent to f,
// emitting definitional clauses as needed. The input is assumed to be in NNF.
func (t *tseitinBuilder) encode(f Formula) string {
	switch f.Conn {
	case ConnAtom:
		return f.Pred
	case ConnTrue:
		v := t.trueVar()
		return v
	case ConnFalse:
		v := t.falseVar()
		return v
	case ConnNot:
		// In NNF the operand is an atom.
		inner := f.Subs[0].Pred
		key := "!" + inner
		if v, ok := t.memo[key]; ok {
			return v
		}
		v := t.fresh()
		// v <-> !inner
		t.emit(NewPClause(NegPLit(v), NegPLit(inner)))
		t.emit(NewPClause(PosPLit(v), PosPLit(inner)))
		t.memo[key] = v
		return v
	case ConnAnd:
		a := t.encodeLit(f.Subs[0])
		b := t.encodeLit(f.Subs[1])
		v := t.fresh()
		// v <-> (a & b)
		t.emit(NewPClause(NegPLit(v), a))
		t.emit(NewPClause(NegPLit(v), b))
		t.emit(NewPClause(PosPLit(v), a.Negated(), b.Negated()))
		return v
	case ConnOr:
		a := t.encodeLit(f.Subs[0])
		b := t.encodeLit(f.Subs[1])
		v := t.fresh()
		// v <-> (a | b)
		t.emit(NewPClause(NegPLit(v), a, b))
		t.emit(NewPClause(PosPLit(v), a.Negated()))
		t.emit(NewPClause(PosPLit(v), b.Negated()))
		return v
	default:
		// ToNNF removed implications and biconditionals; nothing else remains.
		return t.encode(ToNNF(f))
	}
}

// encodeLit encodes a subformula and returns the literal standing for it,
// handling negated atoms directly so no auxiliary is needed for a bare literal.
func (t *tseitinBuilder) encodeLit(f Formula) PLiteral {
	if f.Conn == ConnAtom {
		return PosPLit(f.Pred)
	}
	if f.Conn == ConnNot && f.Subs[0].Conn == ConnAtom {
		return NegPLit(f.Subs[0].Pred)
	}
	return PosPLit(t.encode(f))
}

func (t *tseitinBuilder) trueVar() string {
	if v, ok := t.memo["__T"]; ok {
		return v
	}
	v := t.fresh()
	t.emit(NewPClause(PosPLit(v)))
	t.memo["__T"] = v
	return v
}

func (t *tseitinBuilder) falseVar() string {
	if v, ok := t.memo["__F"]; ok {
		return v
	}
	v := t.fresh()
	t.emit(NewPClause(NegPLit(v)))
	t.memo["__F"] = v
	return v
}

func (t *tseitinBuilder) emit(c PClause) {
	t.clauses = append(t.clauses, c)
}
