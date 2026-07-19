package randommatrix

import (
	"math"
	"sort"
)

// PoissonSpacingDensity returns the density exp(-s) of nearest-neighbour level
// spacings for an uncorrelated (Poisson) spectrum with unit mean spacing.
func PoissonSpacingDensity(s float64) float64 {
	if s < 0 {
		return 0
	}
	return math.Exp(-s)
}

// PoissonSpacingCDF returns the cumulative distribution 1 - exp(-s) of Poisson
// level spacings.
func PoissonSpacingCDF(s float64) float64 {
	if s < 0 {
		return 0
	}
	return 1 - math.Exp(-s)
}

// wignerAB returns the constants a and b of the Wigner surmise density
// p(s) = a s^beta exp(-b s^2) normalised to unit total mass and unit mean.
func wignerAB(beta float64) (a, b float64) {
	g1 := math.Gamma((beta + 1) / 2)
	g2 := math.Gamma((beta + 2) / 2)
	b = (g2 / g1) * (g2 / g1)
	a = 2 * math.Pow(b, (beta+1)/2) / g1
	return a, b
}

// WignerSurmise returns the Wigner surmise density for Dyson index beta at
// spacing s, normalised to unit total mass and unit mean spacing. The classical
// cases are beta = 1 (GOE), beta = 2 (GUE) and beta = 4 (GSE).
func WignerSurmise(s, beta float64) float64 {
	if s < 0 {
		return 0
	}
	a, b := wignerAB(beta)
	return a * math.Pow(s, beta) * math.Exp(-b*s*s)
}

// WignerSurmiseCDF returns the cumulative distribution of the Wigner surmise for
// Dyson index beta at spacing s, integrated numerically.
func WignerSurmiseCDF(s, beta float64) float64 {
	if s <= 0 {
		return 0
	}
	return simpson(func(t float64) float64 { return WignerSurmise(t, beta) }, 0, s, 400)
}

// GOESpacingDensity returns the GOE (beta = 1) Wigner surmise density
// (pi/2) s exp(-pi s^2/4).
func GOESpacingDensity(s float64) float64 {
	if s < 0 {
		return 0
	}
	return (math.Pi / 2) * s * math.Exp(-math.Pi*s*s/4)
}

// GOESpacingCDF returns the GOE spacing distribution 1 - exp(-pi s^2/4).
func GOESpacingCDF(s float64) float64 {
	if s < 0 {
		return 0
	}
	return 1 - math.Exp(-math.Pi*s*s/4)
}

// GUESpacingDensity returns the GUE (beta = 2) Wigner surmise density
// (32/pi^2) s^2 exp(-4 s^2/pi).
func GUESpacingDensity(s float64) float64 {
	if s < 0 {
		return 0
	}
	return (32 / (math.Pi * math.Pi)) * s * s * math.Exp(-4*s*s/math.Pi)
}

// GUESpacingCDF returns the GUE spacing distribution, integrated numerically.
func GUESpacingCDF(s float64) float64 {
	if s <= 0 {
		return 0
	}
	return simpson(GUESpacingDensity, 0, s, 400)
}

// GSESpacingDensity returns the GSE (beta = 4) Wigner surmise density
// (2^18/(3^6 pi^3)) s^4 exp(-64 s^2/(9 pi)).
func GSESpacingDensity(s float64) float64 {
	if s < 0 {
		return 0
	}
	coef := math.Pow(2, 18) / (math.Pow(3, 6) * math.Pow(math.Pi, 3))
	return coef * math.Pow(s, 4) * math.Exp(-64*s*s/(9*math.Pi))
}

// GSESpacingCDF returns the GSE spacing distribution, integrated numerically.
func GSESpacingCDF(s float64) float64 {
	if s <= 0 {
		return 0
	}
	return simpson(GSESpacingDensity, 0, s, 400)
}

// NearestNeighbourSpacings returns the consecutive gaps of the sorted eigenvalue
// sample.
func NearestNeighbourSpacings(eigs []float64) []float64 {
	return SpectralGaps(eigs)
}

// MeanSpacing returns the mean nearest-neighbour spacing of the eigenvalue
// sample.
func MeanSpacing(eigs []float64) float64 {
	g := SpectralGaps(eigs)
	return Mean(g)
}

// NormalizedSpacings returns the nearest-neighbour spacings rescaled by their
// mean, a crude global unfolding that produces unit mean spacing.
func NormalizedSpacings(eigs []float64) []float64 {
	g := SpectralGaps(eigs)
	m := Mean(g)
	if m == 0 || math.IsNaN(m) {
		return g
	}
	out := make([]float64, len(g))
	for i, v := range g {
		out[i] = v / m
	}
	return out
}

// Unfold maps each eigenvalue lambda to n*cdf(lambda), where n is the number of
// eigenvalues and cdf is the integrated limiting spectral density. The
// transformed levels have unit mean spacing.
func Unfold(eigs []float64, cdf func(float64) float64) []float64 {
	s := SortedCopy(eigs)
	n := float64(len(s))
	out := make([]float64, len(s))
	for i, l := range s {
		out[i] = n * cdf(l)
	}
	return out
}

// UnfoldedSpacings returns the nearest-neighbour spacings of the eigenvalues
// after unfolding with the given integrated spectral density.
func UnfoldedSpacings(eigs []float64, cdf func(float64) float64) []float64 {
	return SpectralGaps(Unfold(eigs, cdf))
}

// SpacingRatios returns the ratios s_{i+1}/s_i of consecutive nearest-neighbour
// spacings of the eigenvalue sample.
func SpacingRatios(eigs []float64) []float64 {
	g := SpectralGaps(eigs)
	if len(g) < 2 {
		return []float64{}
	}
	out := make([]float64, len(g)-1)
	for i := 1; i < len(g); i++ {
		if g[i-1] == 0 {
			out[i-1] = math.Inf(1)
		} else {
			out[i-1] = g[i] / g[i-1]
		}
	}
	return out
}

// ConsecutiveSpacingRatios returns the restricted ratios
// r_i = min(s_i, s_{i-1}) / max(s_i, s_{i-1}) in [0, 1], which need no spectral
// unfolding.
func ConsecutiveSpacingRatios(eigs []float64) []float64 {
	g := SpectralGaps(eigs)
	if len(g) < 2 {
		return []float64{}
	}
	out := make([]float64, len(g)-1)
	for i := 1; i < len(g); i++ {
		lo, hi := g[i-1], g[i]
		if lo > hi {
			lo, hi = hi, lo
		}
		if hi == 0 {
			out[i-1] = 0
		} else {
			out[i-1] = lo / hi
		}
	}
	return out
}

// MeanConsecutiveRatio returns the mean of the restricted spacing ratios of the
// eigenvalue sample.
func MeanConsecutiveRatio(eigs []float64) float64 {
	return Mean(ConsecutiveSpacingRatios(eigs))
}

// PoissonMeanRatio returns the mean restricted spacing ratio 2 ln 2 - 1 for a
// Poisson spectrum.
func PoissonMeanRatio() float64 { return 2*math.Ln2 - 1 }

// GOEMeanRatio returns the surmise value 0.5307 of the mean restricted spacing
// ratio for the GOE.
func GOEMeanRatio() float64 { return 0.5307 }

// GUEMeanRatio returns the surmise value 0.5996 of the mean restricted spacing
// ratio for the GUE.
func GUEMeanRatio() float64 { return 0.5996 }

// GSEMeanRatio returns the surmise value 0.6744 of the mean restricted spacing
// ratio for the GSE.
func GSEMeanRatio() float64 { return 0.6744 }

// RatioSurmiseDensity returns the Atas et al. surmise density for the restricted
// spacing ratio r in [0, 1] at Dyson index beta,
// P(r) = (1/Z) (r + r^2)^beta / (1 + r + r^2)^(1 + 3 beta/2), with Z fixed by
// numerical normalisation.
func RatioSurmiseDensity(r, beta float64) float64 {
	if r < 0 || r > 1 {
		return 0
	}
	z := ratioSurmiseNorm(beta)
	return ratioSurmiseUnnorm(r, beta) / z
}

func ratioSurmiseUnnorm(r, beta float64) float64 {
	return math.Pow(r+r*r, beta) / math.Pow(1+r+r*r, 1+1.5*beta)
}

func ratioSurmiseNorm(beta float64) float64 {
	return simpson(func(t float64) float64 { return ratioSurmiseUnnorm(t, beta) }, 0, 1, 2000)
}

// PoissonRatioDensity returns the surmise density 1/(1+r)^2 of the unrestricted
// consecutive-spacing ratio for a Poisson spectrum.
func PoissonRatioDensity(r float64) float64 {
	if r < 0 {
		return 0
	}
	return 1 / ((1 + r) * (1 + r))
}

// SortSpacings returns the ascending-sorted spacings, convenient for building
// empirical spacing distributions.
func SortSpacings(spacings []float64) []float64 {
	out := make([]float64, len(spacings))
	copy(out, spacings)
	sort.Float64s(out)
	return out
}
