package seq

import "math"

// Triangular returns the n-th triangular number Tₙ = n(n+1)/2, the number of
// dots in a triangular arrangement with n dots on a side: 0, 1, 3, 6, 10, 15,
// … n must be non-negative.
func Triangular(n int) uint64 {
	if n < 0 {
		panic("seq: Triangular requires n >= 0")
	}
	u := uint64(n)
	return u * (u + 1) / 2
}

// IsTriangular reports whether x is a triangular number, that is whether
// 8x+1 is a perfect square.
func IsTriangular(x uint64) bool {
	return IsPerfectSquare(8*x + 1)
}

// Square returns the n-th square number n². n must be non-negative.
func Square(n int) uint64 {
	if n < 0 {
		panic("seq: Square requires n >= 0")
	}
	u := uint64(n)
	return u * u
}

// IsSquare reports whether x is a perfect square. It is an alias of
// IsPerfectSquare provided for symmetry with the other figurate predicates.
func IsSquare(x uint64) bool {
	return IsPerfectSquare(x)
}

// Pentagonal returns the n-th pentagonal number Pₙ = n(3n−1)/2: 0, 1, 5, 12,
// 22, 35, 51, … n must be non-negative. For the signed indexing used in
// Euler's pentagonal number theorem see GeneralizedPentagonal.
func Pentagonal(n int) uint64 {
	if n < 0 {
		panic("seq: Pentagonal requires n >= 0")
	}
	u := uint64(n)
	return u * (3*u - 1) / 2
}

// GeneralizedPentagonal returns n(3n−1)/2 evaluated at the signed index n, so
// that ranging n over 0, 1, −1, 2, −2, … yields the generalized pentagonal
// numbers 0, 1, 2, 5, 7, 12, 15, … that appear in Euler's pentagonal number
// theorem. The result is int64 because the index may be negative.
func GeneralizedPentagonal(n int) int64 {
	m := int64(n)
	return m * (3*m - 1) / 2
}

// IsPentagonal reports whether x is a (non-negative-index) pentagonal number,
// that is whether (1+√(24x+1))/6 is a positive integer.
func IsPentagonal(x uint64) bool {
	return IsPolygonal(5, x)
}

// Hexagonal returns the n-th hexagonal number Hₙ = n(2n−1): 0, 1, 6, 15, 28,
// 45, 66, … n must be non-negative.
func Hexagonal(n int) uint64 {
	if n < 0 {
		panic("seq: Hexagonal requires n >= 0")
	}
	u := uint64(n)
	return u * (2*u - 1)
}

// IsHexagonal reports whether x is a hexagonal number, that is whether
// (1+√(8x+1))/4 is a positive integer.
func IsHexagonal(x uint64) bool {
	return IsPolygonal(6, x)
}

// Heptagonal returns the n-th heptagonal number n(5n−3)/2: 0, 1, 7, 18, 34,
// 55, 81, … n must be non-negative.
func Heptagonal(n int) uint64 {
	if n < 0 {
		panic("seq: Heptagonal requires n >= 0")
	}
	u := uint64(n)
	return u * (5*u - 3) / 2
}

// Octagonal returns the n-th octagonal number n(3n−2): 0, 1, 8, 21, 40, 65,
// 96, … n must be non-negative.
func Octagonal(n int) uint64 {
	if n < 0 {
		panic("seq: Octagonal requires n >= 0")
	}
	u := uint64(n)
	return u * (3*u - 2)
}

// Polygonal returns the n-th s-gonal number
//
//	P(s, n) = ((s−2)·n² − (s−4)·n) / 2,
//
// the generalization of the triangular (s=3), square (s=4), pentagonal (s=5)
// and higher polygonal numbers. s must be at least 3 and n must be
// non-negative.
func Polygonal(s, n int) uint64 {
	if s < 3 {
		panic("seq: Polygonal requires s >= 3")
	}
	if n < 0 {
		panic("seq: Polygonal requires n >= 0")
	}
	u := uint64(n)
	return ((uint64(s)-2)*u*u - (uint64(s)-4)*u) / 2
}

// IsPolygonal reports whether x is an s-gonal number for some non-negative
// index, i.e. whether the equation P(s, n) = x has a non-negative integer
// solution n. s must be at least 3.
func IsPolygonal(s int, x uint64) bool {
	if s < 3 {
		panic("seq: IsPolygonal requires s >= 3")
	}
	if x == 0 {
		return true // n = 0
	}
	// Solve (s-2)n^2 - (s-4)n - 2x = 0 for a floating estimate, then verify
	// the neighbouring integer candidates exactly.
	a := float64(s - 2)
	b := -float64(s - 4)
	c := -2 * float64(x)
	disc := b*b - 4*a*c
	if disc < 0 {
		return false
	}
	est := (-b + math.Sqrt(disc)) / (2 * a)
	cand := uint64(math.Round(est))
	for _, cc := range []uint64{cand, cand + 1} {
		if cc >= 1 && Polygonal(s, int(cc)) == x {
			return true
		}
	}
	if cand >= 1 && Polygonal(s, int(cand-1)) == x {
		return true
	}
	return false
}

// CenteredPolygonal returns the n-th centered s-gonal number
// s·n(n+1)/2 + 1, the count of dots in a centered polygonal arrangement with
// n rings around a central dot: for s=4 this is 1, 5, 13, 25, 41, … s must be
// at least 3 and n must be non-negative.
func CenteredPolygonal(s, n int) uint64 {
	if s < 3 {
		panic("seq: CenteredPolygonal requires s >= 3")
	}
	if n < 0 {
		panic("seq: CenteredPolygonal requires n >= 0")
	}
	u := uint64(n)
	return uint64(s)*u*(u+1)/2 + 1
}

// CenteredTriangular returns the n-th centered triangular number
// 3n(n+1)/2 + 1: 1, 4, 10, 19, 31, … n must be non-negative.
func CenteredTriangular(n int) uint64 {
	return CenteredPolygonal(3, n)
}

// CenteredSquare returns the n-th centered square number 2n(n+1) + 1:
// 1, 5, 13, 25, 41, … n must be non-negative.
func CenteredSquare(n int) uint64 {
	return CenteredPolygonal(4, n)
}

// CenteredHexagonal returns the n-th centered hexagonal number 3n(n+1) + 1,
// also known as a hex number: 1, 7, 19, 37, 61, … n must be non-negative.
func CenteredHexagonal(n int) uint64 {
	return CenteredPolygonal(6, n)
}

// Tetrahedral returns the n-th tetrahedral (triangular-pyramidal) number
// n(n+1)(n+2)/6, the number of spheres in a tetrahedral stack of n layers:
// 0, 1, 4, 10, 20, 35, 56, … n must be non-negative.
func Tetrahedral(n int) uint64 {
	if n < 0 {
		panic("seq: Tetrahedral requires n >= 0")
	}
	u := uint64(n)
	return u * (u + 1) * (u + 2) / 6
}

// SquarePyramidal returns the n-th square-pyramidal number
// n(n+1)(2n+1)/6, the number of spheres in a pyramid with a square base of
// side n: 0, 1, 5, 14, 30, 55, 91, … n must be non-negative.
func SquarePyramidal(n int) uint64 {
	if n < 0 {
		panic("seq: SquarePyramidal requires n >= 0")
	}
	u := uint64(n)
	return u * (u + 1) * (2*u + 1) / 6
}

// Pentatope returns the n-th pentatope (4-simplex) number
// n(n+1)(n+2)(n+3)/24, the four-dimensional analogue of the triangular and
// tetrahedral numbers: 0, 1, 5, 15, 35, 70, 126, … n must be non-negative.
func Pentatope(n int) uint64 {
	if n < 0 {
		panic("seq: Pentatope requires n >= 0")
	}
	u := uint64(n)
	return u * (u + 1) * (u + 2) * (u + 3) / 24
}

// Pronic returns the n-th pronic (oblong) number n(n+1), the product of two
// consecutive integers: 0, 2, 6, 12, 20, 30, … n must be non-negative.
func Pronic(n int) uint64 {
	if n < 0 {
		panic("seq: Pronic requires n >= 0")
	}
	u := uint64(n)
	return u * (u + 1)
}

// IsPronic reports whether x is a pronic number, that is whether x = n(n+1)
// for some non-negative integer n.
func IsPronic(x uint64) bool {
	// x = n(n+1) => 4x+1 = (2n+1)^2 is an odd perfect square.
	r := IntSqrt(4*x + 1)
	return r*r == 4*x+1
}

// StarNumber returns the n-th star number (centered 12-gonal number)
// 6n(n−1) + 1: 1, 13, 37, 73, 121, … n must be at least 1.
func StarNumber(n int) uint64 {
	if n < 1 {
		panic("seq: StarNumber requires n >= 1")
	}
	u := uint64(n)
	return 6*u*(u-1) + 1
}

// Cube returns the n-th cube number n³. n must be non-negative.
func Cube(n int) uint64 {
	if n < 0 {
		panic("seq: Cube requires n >= 0")
	}
	u := uint64(n)
	return u * u * u
}

// Gnomonic returns the n-th gnomonic number 2n−1, the n-th odd number and the
// difference between consecutive square numbers: 1, 3, 5, 7, 9, … n must be at
// least 1.
func Gnomonic(n int) uint64 {
	if n < 1 {
		panic("seq: Gnomonic requires n >= 1")
	}
	return 2*uint64(n) - 1
}
