package lattices

import (
	"math"
	"sort"
)

// BallVolume returns the volume of the n-dimensional Euclidean ball of radius r,
// namely pi^(n/2)/Gamma(n/2+1) * r^n.
func BallVolume(n int, r float64) float64 {
	if n < 0 {
		return 0
	}
	return UnitBallVolume(n) * math.Pow(r, float64(n))
}

// UnitBallVolume returns the volume of the n-dimensional unit ball,
// pi^(n/2)/Gamma(n/2+1).
func UnitBallVolume(n int) float64 {
	if n < 0 {
		return 0
	}
	return math.Pow(math.Pi, float64(n)/2) / math.Gamma(float64(n)/2+1)
}

// exactHermite holds the exact Hermite constants gamma_n for 1 <= n <= 8.
var exactHermite = map[int]float64{
	1: 1,
	2: 2.0 / math.Sqrt(3),
	3: math.Cbrt(2),
	4: math.Sqrt2,
	5: math.Pow(8, 1.0/5),
	6: math.Pow(64.0/3, 1.0/6),
	7: math.Pow(64, 1.0/7),
	8: 2,
}

// HermiteConstant returns Hermite's constant gamma_n. For 1 <= n <= 8 the exact
// value is returned; for larger n it returns Minkowski's upper bound (see
// HermiteConstantUpperBound), which is an over-estimate of the true constant.
func HermiteConstant(n int) float64 {
	if v, ok := exactHermite[n]; ok {
		return v
	}
	return HermiteConstantUpperBound(n)
}

// HermiteConstantUpperBound returns Minkowski's upper bound on Hermite's
// constant, gamma_n <= (2/pi) * Gamma(2 + n/2)^(2/n).
func HermiteConstantUpperBound(n int) float64 {
	if n <= 0 {
		return 0
	}
	g := math.Gamma(2 + float64(n)/2)
	return (2 / math.Pi) * math.Pow(g, 2.0/float64(n))
}

// HermiteConstantLowerBound returns a Minkowski-Hlawka style lower bound on
// Hermite's constant, gamma_n >= (n/(2*pi*e)) for large n; it is only a rough
// asymptotic estimate for small n.
func HermiteConstantLowerBound(n int) float64 {
	if n <= 0 {
		return 0
	}
	return float64(n) / (2 * math.Pi * math.E)
}

// MinkowskiBound returns the upper bound on the first minimum lambda_1 given by
// Minkowski's convex body theorem for a lattice of the given rank and
// determinant: lambda_1 <= 2 * (det / V_n)^(1/n), where V_n is the unit ball
// volume.
func MinkowskiBound(n int, det float64) float64 {
	if n <= 0 || det <= 0 {
		return 0
	}
	return 2 * math.Pow(det/UnitBallVolume(n), 1.0/float64(n))
}

// MinkowskiBoundBasis returns MinkowskiBound evaluated at the rank and
// determinant of the basis.
func (b Basis) MinkowskiBoundBasis() float64 {
	return MinkowskiBound(len(b), b.Determinant())
}

// HermiteBound returns the Hermite bound sqrt(gamma_n) * det^(1/n) on the first
// minimum lambda_1 of a rank-n lattice with the given determinant.
func HermiteBound(n int, det float64) float64 {
	if n <= 0 || det <= 0 {
		return 0
	}
	return math.Sqrt(HermiteConstant(n)) * math.Pow(det, 1.0/float64(n))
}

// HermiteBoundBasis returns HermiteBound evaluated at the rank and determinant
// of the basis.
func (b Basis) HermiteBoundBasis() float64 {
	return HermiteBound(len(b), b.Determinant())
}

// GaussianHeuristic returns the Gaussian heuristic estimate of the first
// minimum, (det / V_n)^(1/n), the expected shortest length for a random lattice
// of the given rank and determinant.
func GaussianHeuristic(n int, det float64) float64 {
	if n <= 0 || det <= 0 {
		return 0
	}
	return math.Pow(det/UnitBallVolume(n), 1.0/float64(n))
}

// GaussianHeuristicBasis returns the Gaussian heuristic evaluated for the
// basis.
func (b Basis) GaussianHeuristicBasis() float64 {
	return GaussianHeuristic(len(b), b.Determinant())
}

// MinkowskiSecondUpper returns the upper bound of Minkowski's second theorem on
// the product of the successive minima: prod_i lambda_i <= 2^n * det / V_n.
func MinkowskiSecondUpper(n int, det float64) float64 {
	if n <= 0 || det <= 0 {
		return 0
	}
	return math.Pow(2, float64(n)) * det / UnitBallVolume(n)
}

// MinkowskiSecondLower returns the lower bound of Minkowski's second theorem on
// the product of the successive minima: prod_i lambda_i >= 2^n * det / (n! V_n).
func MinkowskiSecondLower(n int, det float64) float64 {
	if n <= 0 || det <= 0 {
		return 0
	}
	fact := 1.0
	for i := 2; i <= n; i++ {
		fact *= float64(i)
	}
	return math.Pow(2, float64(n)) * det / (fact * UnitBallVolume(n))
}

// RankinBound returns Rankin's bound gamma_{n,m} <= gamma_n^m used to bound the
// minimum determinant of an m-dimensional sublattice; here it returns
// gamma_n^m as a convenience.
func RankinBound(n, m int) float64 {
	return math.Pow(HermiteConstant(n), float64(m))
}

// SuccessiveMinima returns the successive minima lambda_1 <= ... <= lambda_n of
// the lattice, computed exactly by enumerating short vectors and greedily
// collecting linearly independent ones. It returns ErrEmpty or ErrNotFullRank
// as appropriate and ErrNoSolution if the enumeration limit is exhausted before
// n independent vectors are found.
func (b Basis) SuccessiveMinima() ([]float64, error) {
	if len(b) == 0 {
		return nil, ErrEmpty
	}
	red := b.LLLDefault()
	if !red.IsFullRank(1e-9) {
		return nil, ErrNotFullRank
	}
	radius := red.MaxNorm()*(1+1e-7) + 1e-12
	vecs, complete, err := red.EnumerateVectors(radius, 1<<20)
	if err != nil {
		return nil, err
	}
	if !complete {
		return nil, ErrNoSolution
	}
	sort.Slice(vecs, func(i, j int) bool { return vecs[i].Norm2 < vecs[j].Norm2 })
	n := len(b)
	var chosen []Vec
	minima := make([]float64, 0, n)
	for _, v := range vecs {
		cand := append(append([]Vec{}, chosen...), v.Vector)
		if Basis(cand).Matrix().Rank(1e-9) == len(cand) {
			chosen = append(chosen, v.Vector)
			minima = append(minima, v.Norm())
			if len(minima) == n {
				break
			}
		}
	}
	if len(minima) < n {
		return nil, ErrNoSolution
	}
	return minima, nil
}
