package groebner

// Order is a monomial ordering: a function that compares two monomials and
// returns -1 if a < b, 0 if a == b, and +1 if a > b. A valid monomial order is
// a total order that is multiplicative and for which 1 is the smallest
// monomial. All monomials passed to an Order are assumed to have the same
// number of variables.
type Order func(a, b Monomial) int

func sign(x int) int {
	switch {
	case x < 0:
		return -1
	case x > 0:
		return 1
	default:
		return 0
	}
}

// CompareLex compares two monomials by the lexicographic order: the monomial
// with the larger exponent in the first variable where they differ is greater.
func CompareLex(a, b Monomial) int {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			return sign(a[i] - b[i])
		}
	}
	return 0
}

// CompareGrlex compares two monomials by the graded lexicographic order: first
// by total degree, breaking ties with the lexicographic order.
func CompareGrlex(a, b Monomial) int {
	da, db := a.Degree(), b.Degree()
	if da != db {
		return sign(da - db)
	}
	return CompareLex(a, b)
}

// CompareGrevlex compares two monomials by the graded reverse lexicographic
// order: first by total degree, breaking ties by the reverse lexicographic
// rule (the monomial with the smaller exponent in the last variable where they
// differ is greater).
func CompareGrevlex(a, b Monomial) int {
	da, db := a.Degree(), b.Degree()
	if da != db {
		return sign(da - db)
	}
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := n - 1; i >= 0; i-- {
		if a[i] != b[i] {
			return sign(b[i] - a[i])
		}
	}
	return 0
}

// LexOrder is the lexicographic monomial order. It is an elimination order:
// the first variables dominate.
func LexOrder() Order { return CompareLex }

// GrlexOrder is the graded lexicographic monomial order.
func GrlexOrder() Order { return CompareGrlex }

// GrevlexOrder is the graded reverse lexicographic monomial order, generally
// the most efficient order for Gröbner basis computation.
func GrevlexOrder() Order { return CompareGrevlex }

// ReverseOrder returns the order obtained by reversing the given order. Note
// that the reverse of a monomial order is generally not itself a monomial
// order; this helper is intended for sorting and diagnostics.
func ReverseOrder(o Order) Order {
	return func(a, b Monomial) int { return -o(a, b) }
}

// WeightOrder returns a monomial order that first compares the dot product of
// the exponent vectors with the given non-negative integer weights, breaking
// ties with the supplied tie order. For a strictly positive weight vector any
// tie order yields a genuine monomial order.
func WeightOrder(weights []int, tie Order) Order {
	w := append([]int(nil), weights...)
	return func(a, b Monomial) int {
		wa, wb := 0, 0
		for i := range a {
			if i < len(w) {
				wa += a[i] * w[i]
				wb += b[i] * w[i]
			}
		}
		if wa != wb {
			return sign(wa - wb)
		}
		return tie(a, b)
	}
}

// BlockOrder returns the product (block) order that compares the first sep
// variables using first and, on a tie, compares the remaining variables using
// second. Monomials involving any of the first sep variables are larger than
// those that do not, so BlockOrder eliminates the first sep variables.
func BlockOrder(sep int, first, second Order) Order {
	return func(a, b Monomial) int {
		ah := head(a, sep)
		bh := head(b, sep)
		if c := first(ah, bh); c != 0 {
			return c
		}
		return second(tail(a, sep), tail(b, sep))
	}
}

// EliminationOrder returns a monomial order suitable for eliminating the first
// sep variables from an ideal. It is a block order using graded reverse
// lexicographic within each block. The elements of a Gröbner basis with
// respect to this order that do not involve the first sep variables generate
// the corresponding elimination ideal.
func EliminationOrder(sep int) Order {
	return BlockOrder(sep, CompareGrevlex, CompareGrevlex)
}

func head(m Monomial, k int) Monomial {
	if k > len(m) {
		k = len(m)
	}
	return m[:k]
}

func tail(m Monomial, k int) Monomial {
	if k > len(m) {
		k = len(m)
	}
	return m[k:]
}

// MonomialLess reports whether a is strictly smaller than b under the order o.
func MonomialLess(o Order, a, b Monomial) bool { return o(a, b) < 0 }

// MonomialMax returns whichever of a and b is larger under the order o.
func MonomialMax(o Order, a, b Monomial) Monomial {
	if o(a, b) >= 0 {
		return a
	}
	return b
}
