package tropical

import (
	"math"
	"strconv"
)

// Kind enumerates the two tropical semirings supported by this package.
type Kind int

const (
	// MinPlus is the semiring (min, +) with tropical zero +Inf and tropical
	// one 0. Its addition is the minimum and its multiplication is the
	// ordinary sum.
	MinPlus Kind = iota
	// MaxPlus is the semiring (max, +) with tropical zero -Inf and tropical
	// one 0. Its addition is the maximum and its multiplication is the
	// ordinary sum.
	MaxPlus
)

// String returns the human-readable name of the semiring kind.
func (k Kind) String() string {
	switch k {
	case MinPlus:
		return "min-plus"
	case MaxPlus:
		return "max-plus"
	default:
		return "unknown"
	}
}

// Semiring is a value describing which tropical semiring a computation uses.
// It carries no mutable state and is safe to copy and share.
type Semiring struct {
	kind Kind
}

// NewSemiring returns the Semiring for the given kind. It panics if kind is
// neither MinPlus nor MaxPlus.
func NewSemiring(kind Kind) Semiring {
	if kind != MinPlus && kind != MaxPlus {
		panic("tropical: NewSemiring requires MinPlus or MaxPlus")
	}
	return Semiring{kind: kind}
}

// MinPlusSemiring returns the min-plus semiring.
func MinPlusSemiring() Semiring { return Semiring{kind: MinPlus} }

// MaxPlusSemiring returns the max-plus semiring.
func MaxPlusSemiring() Semiring { return Semiring{kind: MaxPlus} }

// Kind returns the kind of the semiring.
func (s Semiring) Kind() Kind { return s.kind }

// String returns the human-readable name of the semiring.
func (s Semiring) String() string { return s.kind.String() }

// IsMinPlus reports whether s is the min-plus semiring.
func (s Semiring) IsMinPlus() bool { return s.kind == MinPlus }

// IsMaxPlus reports whether s is the max-plus semiring.
func (s Semiring) IsMaxPlus() bool { return s.kind == MaxPlus }

// Dual returns the order-dual semiring: MinPlus for MaxPlus and vice versa.
// Negating every entry of a problem and switching to the dual semiring yields
// an equivalent problem.
func (s Semiring) Dual() Semiring {
	if s.kind == MinPlus {
		return Semiring{kind: MaxPlus}
	}
	return Semiring{kind: MinPlus}
}

// Zero returns the tropical zero (additive identity) of the semiring: +Inf for
// min-plus and -Inf for max-plus.
func (s Semiring) Zero() float64 {
	if s.kind == MinPlus {
		return math.Inf(1)
	}
	return math.Inf(-1)
}

// One returns the tropical one (multiplicative identity) of the semiring,
// which is 0 for both min-plus and max-plus.
func (s Semiring) One() float64 { return 0 }

// IsZero reports whether a is the tropical zero of the semiring.
func (s Semiring) IsZero(a float64) bool { return a == s.Zero() }

// IsOne reports whether a equals the tropical one (0).
func (s Semiring) IsOne(a float64) bool { return a == 0 }

// Add returns the tropical sum a (+) b: the minimum for min-plus and the
// maximum for max-plus.
func (s Semiring) Add(a, b float64) float64 {
	if s.kind == MinPlus {
		if a < b {
			return a
		}
		return b
	}
	if a > b {
		return a
	}
	return b
}

// Mul returns the tropical product a (*) b. On finite operands it is the
// ordinary sum a+b; if either operand is the tropical zero the result is the
// tropical zero (which absorbs under multiplication).
func (s Semiring) Mul(a, b float64) float64 {
	z := s.Zero()
	if a == z || b == z {
		return z
	}
	return a + b
}

// Div returns the tropical quotient a (/) b, the residual of multiplication:
// on finite operands it is the ordinary difference a-b. Dividing the tropical
// zero by anything yields the tropical zero, and dividing a finite value by
// the tropical zero yields the opposite infinity (the tropical "top").
func (s Semiring) Div(a, b float64) float64 {
	z := s.Zero()
	if a == z {
		return z
	}
	if b == z {
		if s.kind == MinPlus {
			return math.Inf(-1)
		}
		return math.Inf(1)
	}
	return a - b
}

// Pow returns the tropical power a raised to the integer exponent n, which
// equals n*a for finite a. The zeroth power is the tropical one (0). A
// positive power of the tropical zero is the tropical zero. Negative exponents
// are permitted only for finite a (the semiring is a semifield on finite
// values); a negative power of the tropical zero returns the tropical zero.
func (s Semiring) Pow(a float64, n int) float64 {
	if n == 0 {
		return s.One()
	}
	if s.IsZero(a) {
		return s.Zero()
	}
	return float64(n) * a
}

// Star returns the tropical Kleene star a* = 1 (+) a (+) a^2 (+) ... For
// min-plus this is 0 when a >= 0 and -Inf (divergent) when a < 0; for max-plus
// it is 0 when a <= 0 and +Inf (divergent) when a > 0.
func (s Semiring) Star(a float64) float64 {
	if s.kind == MinPlus {
		if a >= 0 {
			return 0
		}
		return math.Inf(-1)
	}
	if a <= 0 {
		return 0
	}
	return math.Inf(1)
}

// StarConverges reports whether the scalar Kleene star a* is finite: a >= 0 for
// min-plus and a <= 0 for max-plus.
func (s Semiring) StarConverges(a float64) bool {
	if s.kind == MinPlus {
		return a >= 0
	}
	return a <= 0
}

// Sum returns the tropical sum of all arguments, or the tropical zero when no
// arguments are given.
func (s Semiring) Sum(xs ...float64) float64 {
	r := s.Zero()
	for _, x := range xs {
		r = s.Add(r, x)
	}
	return r
}

// Prod returns the tropical product of all arguments, or the tropical one when
// no arguments are given.
func (s Semiring) Prod(xs ...float64) float64 {
	r := s.One()
	for _, x := range xs {
		r = s.Mul(r, x)
	}
	return r
}

// AtLeastAsGood reports whether a is at least as good as b in the sense of
// tropical addition, i.e. whether a (+) b == a. For min-plus this is a <= b and
// for max-plus it is a >= b.
func (s Semiring) AtLeastAsGood(a, b float64) bool {
	return s.Add(a, b) == a
}

// Better returns whichever of a and b tropical addition selects, identical to
// Add but named for readability at call sites that compare candidates.
func (s Semiring) Better(a, b float64) float64 { return s.Add(a, b) }

// FormatScalar renders a tropical scalar, printing the tropical zero as the
// mathematical symbol for infinity with the appropriate sign.
func (s Semiring) FormatScalar(a float64) string {
	if a == math.Inf(1) {
		return "+inf"
	}
	if a == math.Inf(-1) {
		return "-inf"
	}
	return strconv.FormatFloat(a, 'g', -1, 64)
}
