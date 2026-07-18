package logic

// And returns the logical conjunction of a and b: true when both operands are
// true.
func And(a, b bool) bool { return a && b }

// Or returns the logical (inclusive) disjunction of a and b: true when at least
// one operand is true.
func Or(a, b bool) bool { return a || b }

// Not returns the logical negation of a.
func Not(a bool) bool { return !a }

// Xor returns the exclusive disjunction of a and b: true when the operands
// differ.
func Xor(a, b bool) bool { return a != b }

// Nand returns the negated conjunction of a and b, equivalent to Not(And(a, b)).
// Nand is functionally complete: every Boolean function can be built from it
// alone.
func Nand(a, b bool) bool { return !(a && b) }

// Nor returns the negated disjunction of a and b, equivalent to Not(Or(a, b)).
// Nor is functionally complete: every Boolean function can be built from it
// alone.
func Nor(a, b bool) bool { return !(a || b) }

// Xnor returns the negated exclusive disjunction (logical equivalence) of a and
// b: true when the operands are equal. It is identical to Iff.
func Xnor(a, b bool) bool { return a == b }

// Implies returns the material conditional a -> b, which is false in exactly one
// case: a true and b false.
func Implies(a, b bool) bool { return !a || b }

// Iff returns the biconditional a <-> b (logical equivalence): true when the
// operands share the same truth value. It is identical to Xnor.
func Iff(a, b bool) bool { return a == b }

// Majority returns the value taken by the majority of its three inputs: true
// when at least two of a, b and c are true.
func Majority(a, b, c bool) bool {
	return (a && b) || (a && c) || (b && c)
}

// Mux returns a two-to-one multiplexer output: b when sel is true, otherwise a.
func Mux(sel, a, b bool) bool {
	if sel {
		return b
	}
	return a
}

// AndAll returns the conjunction of every value in vals. The empty product is
// true, the identity element for And.
func AndAll(vals ...bool) bool {
	for _, v := range vals {
		if !v {
			return false
		}
	}
	return true
}

// OrAll returns the disjunction of every value in vals. The empty sum is false,
// the identity element for Or.
func OrAll(vals ...bool) bool {
	for _, v := range vals {
		if v {
			return true
		}
	}
	return false
}

// XorAll returns the parity of vals: true when an odd number of the values are
// true. The empty case is false.
func XorAll(vals ...bool) bool {
	r := false
	for _, v := range vals {
		r = r != v
	}
	return r
}
