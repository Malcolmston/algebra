package algebra

import "sort"

// Equal reports structural equality of two expressions. Because the [Add],
// [Mul] and [Pow] constructors canonicalize their inputs, mathematically equal
// expressions built through the public API compare equal here.
func (s *Symbol) Equal(o Expr) bool { x, ok := o.(*Symbol); return ok && x.Name == s.Name }

// Equal reports whether o is the same integer.
func (i *Integer) Equal(o Expr) bool { x, ok := o.(*Integer); return ok && x.Val.Cmp(i.Val) == 0 }

// Equal reports whether o is the same rational.
func (r *Rational) Equal(o Expr) bool { x, ok := o.(*Rational); return ok && x.Val.Cmp(r.Val) == 0 }

// Equal reports whether o is the same float.
func (f *Float) Equal(o Expr) bool { x, ok := o.(*Float); return ok && x.Val == f.Val }

// Equal reports whether o is the same named constant.
func (c *Constant) Equal(o Expr) bool { x, ok := o.(*Constant); return ok && x.Name == c.Name }

// Equal reports whether o is a structurally identical sum.
func (s *sum) Equal(o Expr) bool {
	x, ok := o.(*sum)
	return ok && equalSlice(s.args, x.args)
}

// Equal reports whether o is a structurally identical product.
func (p *product) Equal(o Expr) bool {
	x, ok := o.(*product)
	return ok && equalSlice(p.factors, x.factors)
}

// Equal reports whether o is a structurally identical power.
func (p *power) Equal(o Expr) bool {
	x, ok := o.(*power)
	return ok && p.base.Equal(x.base) && p.exp.Equal(x.exp)
}

// Equal reports whether o is the same function applied to an equal argument.
func (f *fn) Equal(o Expr) bool {
	x, ok := o.(*fn)
	return ok && f.name == x.name && f.arg.Equal(x.arg)
}

// Equal reports whether o is the same two-argument function applied to equal
// arguments.
func (f *fn2) Equal(o Expr) bool {
	x, ok := o.(*fn2)
	return ok && f.name == x.name && f.arg1.Equal(x.arg1) && f.arg2.Equal(x.arg2)
}

// Equal reports whether o is an equal unevaluated integral.
func (n *integral) Equal(o Expr) bool {
	x, ok := o.(*integral)
	return ok && n.arg.Equal(x.arg) && n.v.Equal(x.v)
}

func equalSlice(a, b []Expr) bool {
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

// rank orders the different node kinds for canonical sorting.
func rank(e Expr) int {
	switch e.(type) {
	case *Integer, *Rational, *Float:
		return 0
	case *Constant:
		return 1
	case *Symbol:
		return 2
	case *fn:
		return 3
	case *power:
		return 4
	case *product:
		return 5
	case *sum:
		return 6
	case *fn2:
		return 7
	case *integral:
		return 8
	case *bigOp:
		return 9
	}
	return 10
}

// compareExpr defines a deterministic total order used to canonicalize the
// argument lists of sums and products.
func compareExpr(a, b Expr) int {
	ra, rb := rank(a), rank(b)
	if ra != rb {
		return ra - rb
	}
	switch x := a.(type) {
	case *Integer, *Rational, *Float:
		xr, xok := toRat(a)
		yr, yok := toRat(b)
		if xok && yok {
			return xr.Cmp(yr)
		}
		return cmpFloat(toFloat(a), toFloat(b))
	case *Constant:
		return cmpString(x.Name, b.(*Constant).Name)
	case *Symbol:
		return cmpString(x.Name, b.(*Symbol).Name)
	case *fn:
		y := b.(*fn)
		if c := cmpString(x.name, y.name); c != 0 {
			return c
		}
		return compareExpr(x.arg, y.arg)
	case *power:
		y := b.(*power)
		if c := compareExpr(x.base, y.base); c != 0 {
			return c
		}
		return compareExpr(x.exp, y.exp)
	case *product:
		return compareSlices(x.factors, b.(*product).factors)
	case *sum:
		return compareSlices(x.args, b.(*sum).args)
	case *fn2:
		y := b.(*fn2)
		if c := cmpString(x.name, y.name); c != 0 {
			return c
		}
		if c := compareExpr(x.arg1, y.arg1); c != 0 {
			return c
		}
		return compareExpr(x.arg2, y.arg2)
	case *integral:
		y := b.(*integral)
		if c := compareExpr(x.arg, y.arg); c != 0 {
			return c
		}
		return compareExpr(x.v, y.v)
	case *bigOp:
		y := b.(*bigOp)
		if c := cmpString(x.kind, y.kind); c != 0 {
			return c
		}
		return compareExpr(x.body, y.body)
	}
	return 0
}

func compareSlices(a, b []Expr) int {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		if c := compareExpr(a[i], b[i]); c != 0 {
			return c
		}
	}
	return len(a) - len(b)
}

func cmpString(a, b string) int {
	if a < b {
		return -1
	}
	if a > b {
		return 1
	}
	return 0
}

func cmpFloat(a, b float64) int {
	if a < b {
		return -1
	}
	if a > b {
		return 1
	}
	return 0
}

// sortExprs sorts a slice of expressions into canonical order in place.
func sortExprs(s []Expr) {
	sort.SliceStable(s, func(i, j int) bool { return compareExpr(s[i], s[j]) < 0 })
}
