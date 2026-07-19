package ellipticcurves

import (
	"math"
	"math/big"
)

// NaiveHeightQ returns the naive logarithmic height of the point pt, defined as
// log(max(|num|, |den|)) of its x-coordinate in lowest terms. The point at
// infinity has height 0.
func (c *CurveQ) NaiveHeightQ(pt PointQ) float64 {
	if pt.Infinity {
		return 0
	}
	num := new(big.Int).Abs(pt.X.Num())
	den := new(big.Int).Abs(pt.X.Denom())
	m := num
	if den.Cmp(num) > 0 {
		m = den
	}
	if m.Sign() == 0 {
		return 0
	}
	return bigLog(m)
}

// bigLog returns the natural logarithm of a positive integer, accurate for very
// large values via its bit length and leading mantissa.
func bigLog(n *big.Int) float64 {
	if n.Sign() <= 0 {
		return 0
	}
	bits := n.BitLen()
	if bits <= 62 {
		return math.Log(float64(n.Int64()))
	}
	shift := uint(bits - 52)
	top := new(big.Int).Rsh(n, shift)
	return math.Log(float64(top.Int64())) + float64(shift)*math.Ln2
}

// CanonicalHeightApprox returns an approximation to the Neron-Tate canonical
// height of pt, computed as the limit 4^{-n} * h(x(2^n * pt)) of the naive
// x-coordinate height under repeated doubling. The iterations parameter sets n;
// values around 10-20 converge for small points. The point at infinity has
// canonical height 0.
func (c *CurveQ) CanonicalHeightApprox(pt PointQ, iterations int) float64 {
	if pt.Infinity {
		return 0
	}
	p := clonePointQ(pt)
	scale := 1.0
	h := 0.0
	for n := 0; n < iterations; n++ {
		if p.Infinity {
			return 0
		}
		h = c.NaiveHeightQ(p) * scale
		p = c.Double(p)
		scale /= 4.0
	}
	return h
}

// HeightPairing returns the Neron-Tate height pairing
// <P, Q> = (h(P+Q) - h(P) - h(Q)) / 2 using canonical heights approximated with
// the given number of doubling iterations.
func (c *CurveQ) HeightPairing(p, q PointQ, iterations int) float64 {
	hp := c.CanonicalHeightApprox(p, iterations)
	hq := c.CanonicalHeightApprox(q, iterations)
	hpq := c.CanonicalHeightApprox(c.Add(p, q), iterations)
	return (hpq - hp - hq) / 2
}

// HeightPairingMatrix returns the Gram matrix of the height pairing of the given
// points, an n-by-n symmetric matrix of canonical height pairings.
func (c *CurveQ) HeightPairingMatrix(points []PointQ, iterations int) [][]float64 {
	n := len(points)
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			v := c.HeightPairing(points[i], points[j], iterations)
			m[i][j] = v
			m[j][i] = v
		}
	}
	return m
}

// Regulator returns the determinant of the height-pairing Gram matrix of the
// given points, an approximation to the elliptic regulator. A value bounded away
// from zero indicates that the points are independent in E(Q)/torsion.
func (c *CurveQ) Regulator(points []PointQ, iterations int) float64 {
	return matrixDeterminant(c.HeightPairingMatrix(points, iterations))
}

// AreIndependent reports whether the given points are linearly independent in
// E(Q) modulo torsion, judged by the height-pairing regulator exceeding tol in
// absolute value.
func (c *CurveQ) AreIndependent(points []PointQ, iterations int, tol float64) bool {
	if len(points) == 0 {
		return false
	}
	return math.Abs(c.Regulator(points, iterations)) > tol
}

// RankLowerBound returns a heuristic lower bound on the Mordell-Weil rank of the
// curve given a list of rational points: the numerical rank of the
// height-pairing Gram matrix, i.e. the size of a maximal independent subset. It
// uses Gaussian elimination with the given tolerance.
func (c *CurveQ) RankLowerBound(points []PointQ, iterations int, tol float64) int {
	// Drop torsion points (canonical height ~ 0) up front.
	var nonTorsion []PointQ
	for _, pt := range points {
		if c.CanonicalHeightApprox(pt, iterations) > tol {
			nonTorsion = append(nonTorsion, pt)
		}
	}
	if len(nonTorsion) == 0 {
		return 0
	}
	m := c.HeightPairingMatrix(nonTorsion, iterations)
	return matrixRank(m, tol)
}

// matrixDeterminant returns the determinant of a square float matrix via LU
// decomposition with partial pivoting.
func matrixDeterminant(a [][]float64) float64 {
	n := len(a)
	if n == 0 {
		return 1
	}
	m := make([][]float64, n)
	for i := range a {
		m[i] = append([]float64(nil), a[i]...)
	}
	det := 1.0
	for col := 0; col < n; col++ {
		pivot := col
		best := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if math.Abs(m[r][col]) > best {
				best = math.Abs(m[r][col])
				pivot = r
			}
		}
		if best == 0 {
			return 0
		}
		if pivot != col {
			m[col], m[pivot] = m[pivot], m[col]
			det = -det
		}
		det *= m[col][col]
		for r := col + 1; r < n; r++ {
			f := m[r][col] / m[col][col]
			for cc := col; cc < n; cc++ {
				m[r][cc] -= f * m[col][cc]
			}
		}
	}
	return det
}

// matrixRank returns the numerical rank of a square float matrix using Gaussian
// elimination with the given pivot tolerance.
func matrixRank(a [][]float64, tol float64) int {
	n := len(a)
	if n == 0 {
		return 0
	}
	m := make([][]float64, n)
	maxAbs := 1.0
	for i := range a {
		m[i] = append([]float64(nil), a[i]...)
		for _, v := range a[i] {
			if math.Abs(v) > maxAbs {
				maxAbs = math.Abs(v)
			}
		}
	}
	// Use a threshold relative to the matrix scale so that the numerical noise
	// of the approximate canonical heights does not inflate the rank.
	threshold := tol * maxAbs
	rank := 0
	rowUsed := make([]bool, n)
	for col := 0; col < n; col++ {
		pivot := -1
		best := threshold
		for r := 0; r < n; r++ {
			if !rowUsed[r] && math.Abs(m[r][col]) > best {
				best = math.Abs(m[r][col])
				pivot = r
			}
		}
		if pivot == -1 {
			continue
		}
		rowUsed[pivot] = true
		rank++
		for r := 0; r < n; r++ {
			if r != pivot && math.Abs(m[r][col]) > 0 {
				f := m[r][col] / m[pivot][col]
				for cc := 0; cc < n; cc++ {
					m[r][cc] -= f * m[pivot][cc]
				}
			}
		}
	}
	return rank
}
