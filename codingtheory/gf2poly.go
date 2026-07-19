package codingtheory

// This file provides carry-less (GF(2)) polynomial arithmetic on plain ints.
// A polynomial is encoded so that bit i is the coefficient of x^i; addition
// and subtraction are therefore the exclusive-or of the encodings.

// GF2PolyAdd returns the sum a+b of two GF(2) polynomials, which equals their
// bitwise exclusive-or.
func GF2PolyAdd(a, b int) int { return a ^ b }

// GF2PolySub returns the difference a-b of two GF(2) polynomials. Over GF(2)
// subtraction and addition coincide, so this equals GF2PolyAdd(a, b).
func GF2PolySub(a, b int) int { return a ^ b }

// GF2PolyDegree returns the degree of the GF(2) polynomial p, i.e. the index
// of its highest set bit. The zero polynomial has degree -1.
func GF2PolyDegree(p int) int {
	deg := -1
	for p != 0 {
		p >>= 1
		deg++
	}
	return deg
}

// GF2PolyMul returns the carry-less product a*b of two GF(2) polynomials.
func GF2PolyMul(a, b int) int {
	var result int
	for b != 0 {
		if b&1 != 0 {
			result ^= a
		}
		a <<= 1
		b >>= 1
	}
	return result
}

// GF2PolyDivMod divides the GF(2) polynomial a by the non-zero polynomial b
// and returns the quotient and remainder, so that a = quotient*b + remainder
// with GF2PolyDegree(remainder) < GF2PolyDegree(b). It panics if b is zero.
func GF2PolyDivMod(a, b int) (quotient, remainder int) {
	if b == 0 {
		panic("codingtheory: GF2 polynomial division by zero")
	}
	db := GF2PolyDegree(b)
	remainder = a
	for {
		dr := GF2PolyDegree(remainder)
		if dr < db {
			break
		}
		shift := dr - db
		quotient ^= 1 << shift
		remainder ^= b << shift
	}
	return quotient, remainder
}

// GF2PolyMod returns the remainder of the GF(2) polynomial a modulo the
// non-zero polynomial b.
func GF2PolyMod(a, b int) int {
	_, r := GF2PolyDivMod(a, b)
	return r
}

// GF2PolyDiv returns the quotient of the GF(2) polynomial a divided by the
// non-zero polynomial b.
func GF2PolyDiv(a, b int) int {
	q, _ := GF2PolyDivMod(a, b)
	return q
}

// GF2PolyGCD returns a greatest common divisor of the GF(2) polynomials a and
// b using the Euclidean algorithm. GF2PolyGCD(0,0) is 0.
func GF2PolyGCD(a, b int) int {
	for b != 0 {
		a, b = b, GF2PolyMod(a, b)
	}
	return a
}

// GF2PolyMulMod returns (a*b) mod m for GF(2) polynomials, with m non-zero.
func GF2PolyMulMod(a, b, m int) int {
	return GF2PolyMod(GF2PolyMul(a, b), m)
}

// GF2PolyPowMod returns (base^exp) mod m for GF(2) polynomials using
// square-and-multiply. exp must be non-negative and m non-zero.
func GF2PolyPowMod(base, exp, m int) int {
	result := GF2PolyMod(1, m)
	base = GF2PolyMod(base, m)
	for exp > 0 {
		if exp&1 != 0 {
			result = GF2PolyMulMod(result, base, m)
		}
		base = GF2PolyMulMod(base, base, m)
		exp >>= 1
	}
	return result
}

// GF2PolyEval evaluates the GF(2) polynomial p at the bit value x (0 or 1),
// returning 0 or 1. For x=0 it is the constant term; for x=1 it is the parity
// of the number of terms.
func GF2PolyEval(p, x int) int {
	if x&1 == 0 {
		return p & 1
	}
	// popcount parity
	parity := 0
	for p != 0 {
		parity ^= p & 1
		p >>= 1
	}
	return parity
}

// GF2PolyString renders a GF(2) polynomial in descending powers of x, e.g.
// "x^3 + x + 1". The zero polynomial renders as "0".
func GF2PolyString(p int) string {
	if p == 0 {
		return "0"
	}
	deg := GF2PolyDegree(p)
	var parts []string
	for i := deg; i >= 0; i-- {
		if p&(1<<uint(i)) == 0 {
			continue
		}
		switch i {
		case 0:
			parts = append(parts, "1")
		case 1:
			parts = append(parts, "x")
		default:
			parts = append(parts, "x^"+itoa(i))
		}
	}
	out := parts[0]
	for _, s := range parts[1:] {
		out += " + " + s
	}
	return out
}

// itoa is a tiny helper avoiding a strconv import churn in hot paths.
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
