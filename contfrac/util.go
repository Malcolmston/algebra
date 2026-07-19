package contfrac

import "errors"

// Common sentinel errors returned by the package.
var (
	// ErrZeroDenominator is returned when a rational with denominator zero is
	// supplied to a routine that cannot make sense of it.
	ErrZeroDenominator = errors.New("contfrac: zero denominator")
	// ErrPerfectSquare is returned by routines that require a non-square
	// argument (for example the Pell solvers) when given a perfect square.
	ErrPerfectSquare = errors.New("contfrac: argument is a perfect square")
	// ErrNonPositive is returned when a strictly positive argument is required.
	ErrNonPositive = errors.New("contfrac: argument must be positive")
	// ErrNotProper is returned by Egyptian-fraction routines when the input is
	// not a proper fraction in the open interval (0, 1).
	ErrNotProper = errors.New("contfrac: value must lie strictly between 0 and 1")
)

// GCD returns the non-negative greatest common divisor of a and b using the
// Euclidean algorithm. GCD(0, 0) is 0.
func GCD(a, b int64) int64 {
	if a < 0 {
		a = -a
	}
	if b < 0 {
		b = -b
	}
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// LCM returns the least common multiple of a and b as a non-negative value.
// LCM(a, 0) and LCM(0, b) are 0.
func LCM(a, b int64) int64 {
	if a == 0 || b == 0 {
		return 0
	}
	g := GCD(a, b)
	r := (a / g) * b
	if r < 0 {
		r = -r
	}
	return r
}

// ReduceFraction returns p/q in lowest terms with a positive denominator.
// It panics if q == 0.
func ReduceFraction(p, q int64) (int64, int64) {
	if q == 0 {
		panic("contfrac: ReduceFraction requires q != 0")
	}
	if q < 0 {
		p, q = -p, -q
	}
	g := GCD(p, q)
	if g == 0 {
		return 0, 1
	}
	return p / g, q / g
}

// Isqrt returns the integer square root floor(sqrt(n)) for n >= 0 using a
// Newton iteration that is exact for every non-negative int64. It panics for
// negative n.
func Isqrt(n int64) int64 {
	if n < 0 {
		panic("contfrac: Isqrt of negative number")
	}
	if n < 2 {
		return n
	}
	x := int64(1) << ((bitLen(uint64(n)) + 1) / 2)
	for {
		y := (x + n/x) / 2
		if y >= x {
			break
		}
		x = y
	}
	return x
}

// bitLen returns the number of bits needed to represent v (0 for v == 0).
func bitLen(v uint64) int {
	n := 0
	for v > 0 {
		v >>= 1
		n++
	}
	return n
}

// IsPerfectSquare reports whether n is a perfect square. Negative numbers are
// never perfect squares.
func IsPerfectSquare(n int64) bool {
	if n < 0 {
		return false
	}
	r := Isqrt(n)
	return r*r == n
}

// EulerPhi returns Euler's totient of n for n >= 1: the count of integers in
// [1, n] coprime to n. EulerPhi of a non-positive number is 0.
func EulerPhi(n int64) int64 {
	if n <= 0 {
		return 0
	}
	result := n
	m := n
	for p := int64(2); p*p <= m; p++ {
		if m%p == 0 {
			for m%p == 0 {
				m /= p
			}
			result -= result / p
		}
	}
	if m > 1 {
		result -= result / m
	}
	return result
}

// floorDiv returns floor(a/b) for b > 0. Go's built-in division truncates
// toward zero, which differs for negative a.
func floorDiv(a, b int64) int64 {
	q := a / b
	if a%b != 0 && (a < 0) != (b < 0) {
		q--
	}
	return q
}

// absInt returns the absolute value of x.
func absInt(x int64) int64 {
	if x < 0 {
		return -x
	}
	return x
}
