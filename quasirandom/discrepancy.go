package quasirandom

import (
	"errors"
	"math"
	"sort"
)

// ErrEmptyPointSet is returned when a discrepancy routine is given no points.
var ErrEmptyPointSet = errors.New("quasirandom: empty point set")

// ErrRaggedPointSet is returned when the points of a set do not all share the
// same dimension.
var ErrRaggedPointSet = errors.New("quasirandom: points have differing dimensions")

// ErrOutOfUnitCube is returned when a point coordinate lies outside [0,1].
var ErrOutOfUnitCube = errors.New("quasirandom: coordinate outside the unit cube")

// validatePoints checks that points is a non-empty, rectangular set of vectors
// whose coordinates lie in [0,1], returning the point count and dimension.
func validatePoints(points [][]float64) (int, int, error) {
	n := len(points)
	if n == 0 {
		return 0, 0, ErrEmptyPointSet
	}
	d := len(points[0])
	if d == 0 {
		return 0, 0, ErrDimension
	}
	for _, p := range points {
		if len(p) != d {
			return 0, 0, ErrRaggedPointSet
		}
		for _, x := range p {
			if x < 0 || x > 1 || math.IsNaN(x) {
				return 0, 0, ErrOutOfUnitCube
			}
		}
	}
	return n, d, nil
}

// L2StarDiscrepancy returns the L2 star discrepancy of the point set using
// Warnock's closed-form formula, an O(d N^2) exact evaluation of the L2 norm of
// the local discrepancy over all anchored boxes [0,t). It returns an error when
// the point set is empty, ragged or leaves the unit cube.
func L2StarDiscrepancy(points [][]float64) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	N := float64(n)
	term1 := math.Pow(1.0/3.0, float64(d))
	var s2 float64
	for _, p := range points {
		prod := 1.0
		for _, x := range p {
			prod *= (1 - x*x)
		}
		s2 += prod
	}
	term2 := math.Pow(2, 1-float64(d)) / N * s2
	var s3 float64
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			prod := 1.0
			for k := 0; k < d; k++ {
				prod *= (1 - math.Max(points[i][k], points[j][k]))
			}
			s3 += prod
		}
	}
	term3 := s3 / (N * N)
	val := term1 - term2 + term3
	return math.Sqrt(math.Max(val, 0)), nil
}

// CenteredL2Discrepancy returns Hickernell's centered L2-discrepancy, which is
// invariant under reflection of the point set about the center of the cube and
// accounts for all lower-dimensional projections. It returns an error for an
// invalid point set.
func CenteredL2Discrepancy(points [][]float64) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	N := float64(n)
	term1 := math.Pow(13.0/12.0, float64(d))
	var s2 float64
	for _, p := range points {
		prod := 1.0
		for _, x := range p {
			a := math.Abs(x - 0.5)
			prod *= 1 + 0.5*a - 0.5*a*a
		}
		s2 += prod
	}
	term2 := 2.0 / N * s2
	var s3 float64
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			prod := 1.0
			for k := 0; k < d; k++ {
				ai := math.Abs(points[i][k] - 0.5)
				aj := math.Abs(points[j][k] - 0.5)
				aij := math.Abs(points[i][k] - points[j][k])
				prod *= 1 + 0.5*ai + 0.5*aj - 0.5*aij
			}
			s3 += prod
		}
	}
	term3 := s3 / (N * N)
	val := term1 - term2 + term3
	return math.Sqrt(math.Max(val, 0)), nil
}

// WrapAroundL2Discrepancy returns Hickernell's wrap-around L2-discrepancy, which
// treats the unit cube as a torus and so is invariant under coordinate shifts
// modulo one. It returns an error for an invalid point set.
func WrapAroundL2Discrepancy(points [][]float64) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	N := float64(n)
	term1 := -math.Pow(4.0/3.0, float64(d))
	var s float64
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			prod := 1.0
			for k := 0; k < d; k++ {
				a := math.Abs(points[i][k] - points[j][k])
				prod *= 1.5 - a*(1-a)
			}
			s += prod
		}
	}
	val := term1 + s/(N*N)
	return math.Sqrt(math.Max(val, 0)), nil
}

// SymmetricL2Discrepancy returns Hickernell's symmetric L2-discrepancy, which is
// invariant under reflection of any single coordinate. It returns an error for
// an invalid point set.
func SymmetricL2Discrepancy(points [][]float64) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	N := float64(n)
	term1 := math.Pow(4.0/3.0, float64(d))
	var s2 float64
	for _, p := range points {
		prod := 1.0
		for _, x := range p {
			prod *= 1 + 2*x - 2*x*x
		}
		s2 += prod
	}
	term2 := 2.0 / N * s2
	var s3 float64
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			prod := 1.0
			for k := 0; k < d; k++ {
				prod *= 1 - math.Abs(points[i][k]-points[j][k])
			}
			s3 += prod
		}
	}
	term3 := math.Pow(2, float64(d)) / (N * N) * s3
	val := term1 - term2 + term3
	return math.Sqrt(math.Max(val, 0)), nil
}

// ModifiedL2StarDiscrepancy returns Hickernell's modified L2-star discrepancy,
// which weights the anchored boxes so that low-dimensional projections
// contribute. It returns an error for an invalid point set.
func ModifiedL2StarDiscrepancy(points [][]float64) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	N := float64(n)
	term1 := math.Pow(4.0/3.0, float64(d))
	var s2 float64
	for _, p := range points {
		prod := 1.0
		for _, x := range p {
			prod *= (3 - x*x) / 2
		}
		s2 += prod
	}
	term2 := 2.0 / N * s2
	var s3 float64
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			prod := 1.0
			for k := 0; k < d; k++ {
				prod *= 2 - math.Max(points[i][k], points[j][k])
			}
			s3 += prod
		}
	}
	term3 := s3 / (N * N)
	val := term1 - term2 + term3
	return math.Sqrt(math.Max(val, 0)), nil
}

// StarDiscrepancy1D returns the exact star discrepancy of a one-dimensional
// point set (each point a length-one slice) using the closed form
// 1/(2N) + max_i |x_(i) - (2i-1)/(2N)| over the sorted samples. It returns an
// error for an invalid point set or one whose points are not one-dimensional.
func StarDiscrepancy1D(points [][]float64) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	if d != 1 {
		return 0, ErrDimension
	}
	xs := make([]float64, n)
	for i := range points {
		xs[i] = points[i][0]
	}
	sort.Float64s(xs)
	N := float64(n)
	maxDev := 0.0
	for i := 0; i < n; i++ {
		dev := math.Abs(xs[i] - float64(2*i+1)/(2*N))
		if dev > maxDev {
			maxDev = dev
		}
	}
	return 1/(2*N) + maxDev, nil
}

// StarDiscrepancy returns the exact star discrepancy D*_N of the point set,
// computed by the coordinate-grid method: the supremum of the local
// discrepancy over anchored boxes is attained at a corner whose coordinates are
// drawn from the sample coordinates. The cost grows like N^(d+1), so the
// routine is intended for small sets in low dimension. It returns an error for
// an invalid point set.
func StarDiscrepancy(points [][]float64) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	N := float64(n)
	// Candidate coordinate values per dimension.
	coords := make([][]float64, d)
	coordsWithOne := make([][]float64, d)
	for k := 0; k < d; k++ {
		col := make([]float64, n)
		for i := range points {
			col[i] = points[i][k]
		}
		u := sortedUnique(col)
		coords[k] = u
		withOne := append(append([]float64(nil), u...), 1)
		coordsWithOne[k] = sortedUnique(withOne)
	}
	best := 0.0
	// Upper (closed box): maximize count/N - volume over grid of coords.
	idx := make([]int, d)
	for {
		vol := 1.0
		for k := 0; k < d; k++ {
			vol *= coords[k][idx[k]]
		}
		cnt := 0
		for i := 0; i < n; i++ {
			inside := true
			for k := 0; k < d; k++ {
				if points[i][k] > coords[k][idx[k]] {
					inside = false
					break
				}
			}
			if inside {
				cnt++
			}
		}
		if dev := float64(cnt)/N - vol; dev > best {
			best = dev
		}
		if !odometer(idx, coords) {
			break
		}
	}
	// Lower (open box): maximize volume - count/N over grid including 1.
	for k := range idx {
		idx[k] = 0
	}
	for {
		vol := 1.0
		for k := 0; k < d; k++ {
			vol *= coordsWithOne[k][idx[k]]
		}
		cnt := 0
		for i := 0; i < n; i++ {
			inside := true
			for k := 0; k < d; k++ {
				if points[i][k] >= coordsWithOne[k][idx[k]] {
					inside = false
					break
				}
			}
			if inside {
				cnt++
			}
		}
		if dev := vol - float64(cnt)/N; dev > best {
			best = dev
		}
		if !odometer(idx, coordsWithOne) {
			break
		}
	}
	return best, nil
}

// odometer advances the mixed-radix index vector idx over the ragged grid
// defined by the lengths of grid[k], returning false once it wraps around.
func odometer(idx []int, grid [][]float64) bool {
	for k := 0; k < len(idx); k++ {
		idx[k]++
		if idx[k] < len(grid[k]) {
			return true
		}
		idx[k] = 0
	}
	return false
}

// StarDiscrepancyGrid returns an approximation to the star discrepancy obtained
// by sampling the local discrepancy on a regular grid of res divisions per
// dimension. It is a cheap lower bound on the true star discrepancy that is
// useful in higher dimensions where the exact routine is impractical. It
// returns an error for an invalid point set or res < 1.
func StarDiscrepancyGrid(points [][]float64, res int) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	if res < 1 {
		return 0, ErrNonPositive
	}
	N := float64(n)
	best := 0.0
	idx := make([]int, d)
	t := make([]float64, d)
	for {
		vol := 1.0
		for k := 0; k < d; k++ {
			t[k] = float64(idx[k]+1) / float64(res)
			vol *= t[k]
		}
		cnt := 0
		for i := 0; i < n; i++ {
			inside := true
			for k := 0; k < d; k++ {
				if points[i][k] >= t[k] {
					inside = false
					break
				}
			}
			if inside {
				cnt++
			}
		}
		if dev := math.Abs(float64(cnt)/N - vol); dev > best {
			best = dev
		}
		// advance idx over res^d grid
		carry := true
		for k := 0; k < d && carry; k++ {
			idx[k]++
			if idx[k] < res {
				carry = false
			} else {
				idx[k] = 0
			}
		}
		if carry {
			break
		}
	}
	return best, nil
}

// Bernoulli2 returns the second Bernoulli polynomial B2(x) = x^2 - x + 1/6,
// the building block of the diaphony's reproducing kernel.
func Bernoulli2(x float64) float64 { return x*x - x + 1.0/6.0 }

// Diaphony returns the classical (Zinterhof) diaphony of a point set, a
// Fourier-based uniformity measure with the closed form
// F_N^2 = (1/N^2) sum_{i,j} [ prod_k (1 + 2 pi^2 B2({x_ik - x_jk})) - 1 ],
// where {} denotes the fractional part. Like the L2 discrepancies it vanishes
// only in the limit of perfect equidistribution. It returns an error for an
// invalid point set.
func Diaphony(points [][]float64) (float64, error) {
	n, d, err := validatePoints(points)
	if err != nil {
		return 0, err
	}
	N := float64(n)
	twoPi2 := 2 * math.Pi * math.Pi
	var sum float64
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			prod := 1.0
			for k := 0; k < d; k++ {
				t := Frac(points[i][k] - points[j][k])
				prod *= 1 + twoPi2*Bernoulli2(t)
			}
			sum += prod - 1
		}
	}
	val := sum / (N * N)
	return math.Sqrt(math.Max(val, 0)), nil
}

// L2StarDiscrepancySquared returns the square of the L2 star discrepancy,
// avoiding the final square root when only the squared quantity is needed. It
// returns an error for an invalid point set.
func L2StarDiscrepancySquared(points [][]float64) (float64, error) {
	v, err := L2StarDiscrepancy(points)
	if err != nil {
		return 0, err
	}
	return v * v, nil
}

// CenteredL2DiscrepancySquared returns the square of the centered
// L2-discrepancy. It returns an error for an invalid point set.
func CenteredL2DiscrepancySquared(points [][]float64) (float64, error) {
	v, err := CenteredL2Discrepancy(points)
	if err != nil {
		return 0, err
	}
	return v * v, nil
}

// KoksmaHlawkaBound returns the Koksma–Hlawka upper bound on the absolute error
// of an equal-weight cubature rule, namely variation * discrepancy, given the
// Hardy–Krause variation of the integrand and the star discrepancy of the
// nodes. Both arguments must be non-negative.
func KoksmaHlawkaBound(variation, starDiscrepancy float64) float64 {
	return variation * starDiscrepancy
}
