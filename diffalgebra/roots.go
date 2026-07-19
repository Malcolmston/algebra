package diffalgebra

import (
	"math"
	"math/cmplx"
	"math/rand"
	"sort"
)

// durandKerner finds all complex roots of the polynomial whose coefficients are
// given in ascending degree order (coeffs[len-1] is the leading coefficient).
// The iteration is seeded so results are reproducible.
func durandKerner(coeffs []complex128, seed int64) []complex128 {
	n := len(coeffs) - 1
	if n < 1 {
		return nil
	}
	// Normalise to monic form.
	lead := coeffs[n]
	mon := make([]complex128, n+1)
	for i := range coeffs {
		mon[i] = coeffs[i] / lead
	}
	eval := func(z complex128) complex128 {
		acc := complex(0, 0)
		for i := n; i >= 0; i-- {
			acc = acc*z + mon[i]
		}
		return acc
	}
	rng := rand.New(rand.NewSource(seed))
	roots := make([]complex128, n)
	// Spread initial guesses around a circle of a data-dependent radius.
	radius := 1.0 + rootUpperBound(mon)
	for i := 0; i < n; i++ {
		ang := 2*math.Pi*float64(i)/float64(n) + 0.4*rng.Float64()
		r := radius * (0.5 + 0.5*rng.Float64())
		roots[i] = cmplx.Rect(r, ang)
	}
	for iter := 0; iter < 2000; iter++ {
		maxDelta := 0.0
		for i := 0; i < n; i++ {
			denom := complex(1, 0)
			for j := 0; j < n; j++ {
				if j != i {
					denom *= roots[i] - roots[j]
				}
			}
			if denom == 0 {
				denom = complex(1e-14, 0)
			}
			delta := eval(roots[i]) / denom
			roots[i] -= delta
			if d := cmplx.Abs(delta); d > maxDelta {
				maxDelta = d
			}
		}
		if maxDelta < 1e-14 {
			break
		}
	}
	// Clean tiny imaginary/real parts.
	for i := range roots {
		re, im := real(roots[i]), imag(roots[i])
		if math.Abs(im) < 1e-9 {
			im = 0
		}
		if math.Abs(re) < 1e-12 {
			re = 0
		}
		roots[i] = complex(re, im)
	}
	sort.Slice(roots, func(a, b int) bool {
		if math.Abs(real(roots[a])-real(roots[b])) > 1e-9 {
			return real(roots[a]) < real(roots[b])
		}
		return imag(roots[a]) < imag(roots[b])
	})
	return roots
}

// rootUpperBound returns a Cauchy-style upper bound on the magnitude of the
// roots of the monic polynomial mon.
func rootUpperBound(mon []complex128) float64 {
	n := len(mon) - 1
	max := 0.0
	for i := 0; i < n; i++ {
		if a := cmplx.Abs(mon[i]); a > max {
			max = a
		}
	}
	return max
}

// RootCluster is a complex root together with the multiplicity inferred by
// grouping numerically coincident roots.
type RootCluster struct {
	Value complex128
	Mult  int
}

// clusterRoots groups roots that agree within tol into multiplicity clusters,
// averaging each cluster. The result is sorted by increasing real then
// imaginary part.
func clusterRoots(roots []complex128, tol float64) []RootCluster {
	used := make([]bool, len(roots))
	var out []RootCluster
	for i := range roots {
		if used[i] {
			continue
		}
		used[i] = true
		sum := roots[i]
		count := 1
		for j := i + 1; j < len(roots); j++ {
			if used[j] {
				continue
			}
			if cmplx.Abs(roots[j]-roots[i]) < tol {
				used[j] = true
				sum += roots[j]
				count++
			}
		}
		out = append(out, RootCluster{Value: sum / complex(float64(count), 0), Mult: count})
	}
	sort.Slice(out, func(a, b int) bool {
		if math.Abs(real(out[a].Value)-real(out[b].Value)) > tol {
			return real(out[a].Value) < real(out[b].Value)
		}
		return imag(out[a].Value) < imag(out[b].Value)
	})
	return out
}

// PolyComplexRootClusters returns the clustered roots of p with inferred
// multiplicities, using tol to merge coincident roots.
func (p Poly) PolyComplexRootClusters(seed int64, tol float64) []RootCluster {
	return clusterRoots(p.ComplexRootsFloat(seed), tol)
}
