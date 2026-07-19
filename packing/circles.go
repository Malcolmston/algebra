package packing

import "math"

// ----------------------------------------------------------------------------
// Planar (2-D) packing and covering constants.
// ----------------------------------------------------------------------------

// HexagonalPackingDensity returns the density of the densest packing of equal
// circles in the plane, pi/sqrt(12) = pi/(2*sqrt(3)) ~ 0.9069, achieved by the
// hexagonal (triangular) lattice. This is optimal (Thue, Fejes Toth).
func HexagonalPackingDensity() float64 { return math.Pi / math.Sqrt(12) }

// SquarePackingDensity returns the density of equal circles packed on the
// square lattice, pi/4 ~ 0.7854.
func SquarePackingDensity() float64 { return math.Pi / 4 }

// HexagonalCoveringThickness returns the thickness of the thinnest covering of
// the plane by equal circles, 2*pi/sqrt(27) ~ 1.2092, achieved by the
// hexagonal lattice (Kershner). Thickness 1 would be an exact tiling.
func HexagonalCoveringThickness() float64 { return 2 * math.Pi / math.Sqrt(27) }

// SquareCoveringThickness returns the thickness of the covering of the plane by
// equal circles centered on the square lattice, pi/2 ~ 1.5708.
func SquareCoveringThickness() float64 { return math.Pi / 2 }

// HexagonalKissingNumber returns the kissing number of the hexagonal circle
// packing, which is 6.
func HexagonalKissingNumber() int { return 6 }

// SquareKissingNumber returns the kissing number of the square circle packing,
// which is 4.
func SquareKissingNumber() int { return 4 }

// ----------------------------------------------------------------------------
// n equal circles packed inside a unit square.
//
// The problem "pack n points in a unit square maximizing the minimum pairwise
// distance" is equivalent to packing n equal circles in a square: if the
// optimal point spread (maximum minimum distance) is m(n), then n equal circles
// of radius r fit in the unit square with r = m / (2(1+m)), because the centers
// live in the concentric square of side 1-2r and must be pairwise 2r apart.
// ----------------------------------------------------------------------------

// squareSpread holds best known point spreads m(n): the maximum over all
// placements of n points in the closed unit square of the minimum pairwise
// distance. Entries up to n = 9 are proven optimal; n = 10 is the best known
// value. Grid values for perfect squares are handled analytically.
var squareSpread = map[int]float64{
	2:  math.Sqrt2,                      // sqrt(2)
	3:  math.Sqrt(6) - math.Sqrt2,       // sqrt(6)-sqrt(2)
	4:  1,                               // unit square corners
	5:  math.Sqrt2 / 2,                  // 1/sqrt(2)
	6:  math.Sqrt(13) / 6,               // sqrt(13)/6
	7:  2 * (2 - math.Sqrt(3)),          // 4-2sqrt(3)
	8:  (math.Sqrt(6) - math.Sqrt2) / 2, // (sqrt(6)-sqrt(2))/2
	9:  0.5,                             // 1/2
	10: 0.421279543983903,               // best known
}

// CircleInSquareSpread returns the best known optimal point spread m(n) for n
// points in the unit square (the largest achievable minimum pairwise distance),
// together with ok reporting whether a value is available. For a perfect square
// n = k^2 the k-by-k grid gives the proven optimum 1/(k-1); tabulated values
// cover 2 <= n <= 10. For n = 1 the spread is undefined and ok is false.
func CircleInSquareSpread(n int) (m float64, ok bool) {
	if n < 2 {
		return 0, false
	}
	if v, found := squareSpread[n]; found {
		return v, true
	}
	if k := isqrt(n); k*k == n && k >= 2 {
		return 1 / float64(k-1), true
	}
	return 0, false
}

// GridSpread returns the point spread of the k-by-k square grid of k^2 points
// in the unit square, 1/(k-1). This is the proven optimum when n = k^2.
func GridSpread(k int) float64 {
	if k < 2 {
		return math.NaN()
	}
	return 1 / float64(k-1)
}

// CircleInSquareRadius returns the radius of n equal circles packed in a unit
// square using the best known spread from [CircleInSquareSpread]. For n = 1 the
// single circle has radius 1/2. It returns ok = false when no packing is known.
func CircleInSquareRadius(n int) (r float64, ok bool) {
	if n == 1 {
		return 0.5, true
	}
	m, found := CircleInSquareSpread(n)
	if !found {
		return 0, false
	}
	return m / (2 * (1 + m)), true
}

// CircleInSquareDensity returns the fraction of a unit square covered by n
// equal circles packed with the best known radius, n*pi*r^2. It returns ok =
// false when no packing is known.
func CircleInSquareDensity(n int) (density float64, ok bool) {
	r, found := CircleInSquareRadius(n)
	if !found {
		return 0, false
	}
	return float64(n) * math.Pi * r * r, true
}

// CircleInSquareKnown reports whether a best known packing of n equal circles
// in a square is tabulated (including the analytic grid and single-circle
// cases).
func CircleInSquareKnown(n int) bool {
	if n == 1 {
		return true
	}
	_, ok := CircleInSquareSpread(n)
	return ok
}

// ----------------------------------------------------------------------------
// n equal circles packed inside a larger circle.
// ----------------------------------------------------------------------------

// circleRatio holds best known ratios R/r of the enclosing radius R to the
// small-circle radius r for n unit circles packed in a circle. Values up to
// n = 9 (with n = 6 and n = 7 both equal to 3) are proven optimal.
var circleRatio = map[int]float64{
	1: 1,
	2: 2,
	3: 1 + 2/math.Sqrt(3),
	4: 1 + math.Sqrt2,
	5: 1 + 1/math.Sin(math.Pi/5),
	6: 3,
	7: 3,
	8: 1 + 1/math.Sin(math.Pi/7),
	9: 1 + 1/math.Sin(math.Pi/8),
}

// CircleInCircleRatio returns the best known ratio R/r for packing n equal
// circles of radius r inside a circle of radius R, together with ok reporting
// whether a value is tabulated. A smaller ratio is a tighter packing.
func CircleInCircleRatio(n int) (ratio float64, ok bool) {
	v, found := circleRatio[n]
	return v, found
}

// CircleInCircleDensity returns the fraction of the enclosing circle covered by
// n equal circles packed with the best known ratio, n / (R/r)^2. It returns
// ok = false when no packing is tabulated.
func CircleInCircleDensity(n int) (density float64, ok bool) {
	ratio, found := CircleInCircleRatio(n)
	if !found {
		return 0, false
	}
	return float64(n) / (ratio * ratio), true
}

// CircleInCircleKnown reports whether a best known packing of n equal circles
// in a circle is tabulated.
func CircleInCircleKnown(n int) bool {
	_, ok := circleRatio[n]
	return ok
}

// RingRatio returns the enclosing ratio R/r for k equal circles arranged in a
// single symmetric ring touching the boundary of the outer circle,
// 1 + 1/sin(pi/k). This is the optimal packing for 2 <= k <= 6 (and, together
// with a central circle, for k+1 = 7). For k = 1 the ratio is 1.
func RingRatio(k int) float64 {
	if k <= 1 {
		return 1
	}
	return 1 + 1/math.Sin(math.Pi/float64(k))
}

// RingDensity returns the packing density k / RingRatio(k)^2 of a single ring
// of k circles inside the enclosing circle.
func RingDensity(k int) float64 {
	ratio := RingRatio(k)
	return float64(k) / (ratio * ratio)
}

// isqrt returns the integer square root of a non-negative integer n (the floor
// of sqrt(n)).
func isqrt(n int) int {
	if n < 0 {
		return 0
	}
	k := int(math.Sqrt(float64(n)))
	for (k+1)*(k+1) <= n {
		k++
	}
	for k*k > n {
		k--
	}
	return k
}
