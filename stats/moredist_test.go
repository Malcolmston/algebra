package stats

import (
	"testing"
)

func TestGeometricDist(t *testing.T) {
	g := Geometric{P: 0.5}
	approx(t, g.PMF(1), 0.5, tol, "geom pmf1")
	approx(t, g.PMF(2), 0.25, tol, "geom pmf2")
	approx(t, g.PMF(3), 0.125, tol, "geom pmf3")
	approx(t, g.PMF(0), 0, 0, "geom pmf0")
	approx(t, g.CDF(2), 0.75, tol, "geom cdf2")
	approx(t, g.Mean(), 2, tol, "geom mean")
	approx(t, g.Variance(), 2, tol, "geom var")
	// Quantile is the smallest k with CDF(k) >= p.
	approx(t, Geometric{P: 0.5}.Quantile(0.75), 2, tol, "geom quantile")
}

func TestNegativeBinomialDist(t *testing.T) {
	nb := NegativeBinomial{R: 3, P: 0.5}
	approx(t, nb.PMF(0), 0.125, tol, "nbinom pmf0")
	// PMF should sum with CDF consistently.
	approx(t, nb.CDF(0), nb.PMF(0), tol, "nbinom cdf0")
	var sum float64
	for k := 0; k <= 3; k++ {
		sum += nb.PMF(k)
	}
	approx(t, nb.CDF(3), sum, 1e-9, "nbinom cdf3")
	approx(t, nb.Mean(), 3, tol, "nbinom mean")
	approx(t, nb.Variance(), 6, tol, "nbinom var")
}
