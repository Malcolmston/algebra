package chaos

import (
	"fmt"
	"math"
	"sort"
	"strings"
)

// FitLine returns the least-squares slope and intercept of the line y = m x + b
// through the paired samples xs, ys.
func FitLine(xs, ys []float64) (m, b float64) {
	n := float64(len(xs))
	if n == 0 {
		return 0, 0
	}
	var sx, sy, sxx, sxy float64
	for i := range xs {
		sx += xs[i]
		sy += ys[i]
		sxx += xs[i] * xs[i]
		sxy += xs[i] * ys[i]
	}
	den := n*sxx - sx*sx
	if den == 0 {
		return 0, sy / n
	}
	m = (n*sxy - sx*sy) / den
	b = (sy - m*sx) / n
	return m, b
}

// FitSlope returns just the least-squares slope of ys against xs.
func FitSlope(xs, ys []float64) float64 {
	m, _ := FitLine(xs, ys)
	return m
}

// BoxCount returns the number of occupied boxes of side eps needed to cover the
// set of points, where each point is a Vec of common dimension. Points are
// binned by integer box index in each coordinate.
func BoxCount(points []Vec, eps float64) int {
	if eps <= 0 || len(points) == 0 {
		return 0
	}
	seen := make(map[string]struct{})
	var sb strings.Builder
	for _, p := range points {
		sb.Reset()
		for _, c := range p {
			fmt.Fprintf(&sb, "%d,", int(math.Floor(c/eps)))
		}
		seen[sb.String()] = struct{}{}
	}
	return len(seen)
}

// BoxCountingDimension estimates the box-counting (capacity) dimension of a set
// of points by counting occupied boxes over a geometric sequence of box sizes
// from epsMax down to epsMin and fitting log N(eps) against log(1/eps). It
// returns the estimated dimension and the (log(1/eps), log N) samples used.
func BoxCountingDimension(points []Vec, epsMin, epsMax float64, scales int) (dim float64, logInvEps, logN []float64) {
	if scales < 2 || epsMin <= 0 || epsMax <= epsMin {
		return 0, nil, nil
	}
	logInvEps = make([]float64, 0, scales)
	logN = make([]float64, 0, scales)
	ratio := math.Pow(epsMax/epsMin, 1/float64(scales-1))
	eps := epsMax
	for i := 0; i < scales; i++ {
		n := BoxCount(points, eps)
		if n > 0 {
			logInvEps = append(logInvEps, math.Log(1/eps))
			logN = append(logN, math.Log(float64(n)))
		}
		eps /= ratio
	}
	dim = FitSlope(logInvEps, logN)
	return dim, logInvEps, logN
}

// CorrelationSum returns the Grassberger-Procaccia correlation sum C(r): the
// fraction of distinct point pairs whose Euclidean distance is below r.
func CorrelationSum(points []Vec, r float64) float64 {
	n := len(points)
	if n < 2 {
		return 0
	}
	var count int
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if points[i].Distance(points[j]) < r {
				count++
			}
		}
	}
	pairs := n * (n - 1) / 2
	return float64(count) / float64(pairs)
}

// CorrelationDimension estimates the correlation dimension D2 by the
// Grassberger-Procaccia method: it evaluates the correlation sum over a
// geometric sequence of radii and fits log C(r) against log r. It returns the
// estimated dimension and the (log r, log C) samples.
func CorrelationDimension(points []Vec, rMin, rMax float64, scales int) (dim float64, logR, logC []float64) {
	if scales < 2 || rMin <= 0 || rMax <= rMin {
		return 0, nil, nil
	}
	pairDist := allPairDistances(points)
	sort.Float64s(pairDist)
	ratio := math.Pow(rMax/rMin, 1/float64(scales-1))
	r := rMin
	total := float64(len(pairDist))
	for i := 0; i < scales; i++ {
		// Number of pair distances below r via binary search.
		k := sort.SearchFloat64s(pairDist, r)
		if k > 0 && total > 0 {
			c := float64(k) / total
			logR = append(logR, math.Log(r))
			logC = append(logC, math.Log(c))
		}
		r *= ratio
	}
	dim = FitSlope(logR, logC)
	return dim, logR, logC
}

// allPairDistances returns all n(n-1)/2 pairwise Euclidean distances.
func allPairDistances(points []Vec) []float64 {
	n := len(points)
	out := make([]float64, 0, n*(n-1)/2)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			out = append(out, points[i].Distance(points[j]))
		}
	}
	return out
}

// InformationDimension estimates the information (D1) dimension of a set of
// points using box occupancy probabilities: for each box size it forms the
// Shannon entropy of the occupancy distribution and fits it against
// log(1/eps).
func InformationDimension(points []Vec, epsMin, epsMax float64, scales int) (dim float64, logInvEps, entropy []float64) {
	if scales < 2 || epsMin <= 0 || epsMax <= epsMin {
		return 0, nil, nil
	}
	ratio := math.Pow(epsMax/epsMin, 1/float64(scales-1))
	eps := epsMax
	for i := 0; i < scales; i++ {
		counts := boxCounts(points, eps)
		total := 0
		for _, c := range counts {
			total += c
		}
		if total > 0 && len(counts) > 0 {
			var h float64
			for _, c := range counts {
				p := float64(c) / float64(total)
				h -= p * math.Log(p)
			}
			logInvEps = append(logInvEps, math.Log(1/eps))
			entropy = append(entropy, h)
		}
		eps /= ratio
	}
	dim = FitSlope(logInvEps, entropy)
	return dim, logInvEps, entropy
}

// boxCounts returns the occupancy counts of every non-empty box of side eps.
func boxCounts(points []Vec, eps float64) map[string]int {
	m := make(map[string]int)
	var sb strings.Builder
	for _, p := range points {
		sb.Reset()
		for _, c := range p {
			fmt.Fprintf(&sb, "%d,", int(math.Floor(c/eps)))
		}
		m[sb.String()]++
	}
	return m
}

// GeneralizedDimension estimates the Renyi generalized dimension D_q of order q
// (q != 1) from box occupancy probabilities. For q=0 it reduces to the
// box-counting dimension and for q=2 to the correlation dimension.
func GeneralizedDimension(points []Vec, q, epsMin, epsMax float64, scales int) float64 {
	if math.Abs(q-1) < 1e-9 {
		d, _, _ := InformationDimension(points, epsMin, epsMax, scales)
		return d
	}
	if scales < 2 || epsMin <= 0 || epsMax <= epsMin {
		return 0
	}
	ratio := math.Pow(epsMax/epsMin, 1/float64(scales-1))
	eps := epsMax
	var logInvEps, logSum []float64
	for i := 0; i < scales; i++ {
		counts := boxCounts(points, eps)
		total := 0
		for _, c := range counts {
			total += c
		}
		if total > 0 && len(counts) > 0 {
			var s float64
			for _, c := range counts {
				p := float64(c) / float64(total)
				s += math.Pow(p, q)
			}
			logInvEps = append(logInvEps, math.Log(1/eps))
			logSum = append(logSum, math.Log(s))
		}
		eps /= ratio
	}
	// log S(eps) ~ (q-1) D_q log eps, and logInvEps = -log eps, so the slope
	// of log S against log(1/eps) equals -(q-1) D_q, giving D_q below.
	slope := FitSlope(logInvEps, logSum)
	return -slope / (q - 1)
}

// TakensEmbedding reconstructs a delay-coordinate embedding of a scalar time
// series with embedding dimension m and delay tau, returning the sequence of
// m-dimensional delay vectors.
func TakensEmbedding(series []float64, m, tau int) []Vec {
	if m < 1 || tau < 1 {
		return nil
	}
	span := (m - 1) * tau
	if len(series) <= span {
		return nil
	}
	out := make([]Vec, 0, len(series)-span)
	for i := 0; i+span < len(series); i++ {
		v := make(Vec, m)
		for k := 0; k < m; k++ {
			v[k] = series[i+k*tau]
		}
		out = append(out, v)
	}
	return out
}
