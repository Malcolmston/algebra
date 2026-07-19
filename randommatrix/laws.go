package randommatrix

import (
	"math"
	"math/cmplx"
)

// SemicircleRadius returns the radius R = 2*sqrt(variance) of the Wigner
// semicircle law whose second moment (variance) equals the given value.
func SemicircleRadius(variance float64) float64 {
	return 2 * math.Sqrt(variance)
}

// SemicircleVariance returns the variance R^2/4 of the semicircle law of radius
// R.
func SemicircleVariance(radius float64) float64 {
	return radius * radius / 4
}

// SemicircleSupport returns the endpoints (-R, R) of the semicircle law of
// radius R.
func SemicircleSupport(radius float64) (lo, hi float64) {
	return -radius, radius
}

// SemicircleDensity returns the Wigner semicircle density
// (2/(pi R^2)) sqrt(R^2 - x^2) at x for radius R; it is zero outside [-R, R].
func SemicircleDensity(x, radius float64) float64 {
	if radius <= 0 {
		return math.NaN()
	}
	if x < -radius || x > radius {
		return 0
	}
	return 2 * math.Sqrt(radius*radius-x*x) / (math.Pi * radius * radius)
}

// SemicircleCDF returns the cumulative distribution function of the semicircle
// law of radius R evaluated at x.
func SemicircleCDF(x, radius float64) float64 {
	if radius <= 0 {
		return math.NaN()
	}
	if x <= -radius {
		return 0
	}
	if x >= radius {
		return 1
	}
	r2 := radius * radius
	return 0.5 + x*math.Sqrt(r2-x*x)/(math.Pi*r2) + math.Asin(x/radius)/math.Pi
}

// SemicircleMoment returns the k-th moment of the semicircle law of radius R.
// Odd moments vanish; the (2m)-th moment equals the m-th Catalan number times
// (R/2)^(2m).
func SemicircleMoment(k int, radius float64) float64 {
	if k < 0 {
		return math.NaN()
	}
	if k%2 == 1 {
		return 0
	}
	m := k / 2
	return CatalanNumber(m) * math.Pow(radius/2, float64(k))
}

// SemicirclePeak returns the maximum value 2/(pi R) of the semicircle density,
// attained at x = 0.
func SemicirclePeak(radius float64) float64 {
	return 2 / (math.Pi * radius)
}

// SemicircleStieltjes returns the Stieltjes transform m(z) of the semicircle
// law of radius R, the value of the integral of 1/(x-z) against the density.
// The branch of the square root is chosen so that m(z) ~ -1/z as z -> infinity.
func SemicircleStieltjes(z complex128, radius float64) complex128 {
	r2 := complex(radius*radius, 0)
	root := cmplx.Sqrt(z*z - r2)
	// Choose the branch with the same imaginary-part sign as z for Im(z) != 0.
	if imag(z) > 0 && imag(root) < 0 || imag(z) < 0 && imag(root) > 0 {
		root = -root
	}
	sig2 := complex(radius*radius/4, 0)
	return (-z + root) / (2 * sig2)
}

// MarchenkoPasturSupport returns the edges lambda_minus and lambda_plus of the
// Marchenko-Pastur law with aspect ratio c and per-entry variance sigma2.
func MarchenkoPasturSupport(c, sigma2 float64) (lo, hi float64) {
	s := math.Sqrt(c)
	lo = sigma2 * (1 - s) * (1 - s)
	hi = sigma2 * (1 + s) * (1 + s)
	return lo, hi
}

// MarchenkoPasturDensity returns the continuous part of the Marchenko-Pastur
// density with aspect ratio c and variance sigma2 at x. It is zero outside the
// bulk support and does not include the possible atom at zero when c > 1.
func MarchenkoPasturDensity(x, c, sigma2 float64) float64 {
	if c <= 0 || sigma2 <= 0 {
		return math.NaN()
	}
	lo, hi := MarchenkoPasturSupport(c, sigma2)
	if x < lo || x > hi || x == 0 {
		return 0
	}
	return math.Sqrt((hi-x)*(x-lo)) / (2 * math.Pi * sigma2 * c * x)
}

// MarchenkoPasturHasAtom reports whether the Marchenko-Pastur law with aspect
// ratio c places an atom at zero (which happens exactly when c > 1).
func MarchenkoPasturHasAtom(c float64) bool { return c > 1 }

// MarchenkoPasturAtomMass returns the mass max(0, 1 - 1/c) of the atom at zero
// of the Marchenko-Pastur law with aspect ratio c.
func MarchenkoPasturAtomMass(c float64) float64 {
	if c <= 1 {
		return 0
	}
	return 1 - 1/c
}

// MarchenkoPasturCDF returns the cumulative distribution function of the
// Marchenko-Pastur law with aspect ratio c and variance sigma2 at x, including
// the atom at zero when c > 1. The continuous part is integrated numerically.
func MarchenkoPasturCDF(x, c, sigma2 float64) float64 {
	atom := MarchenkoPasturAtomMass(c)
	lo, hi := MarchenkoPasturSupport(c, sigma2)
	if x < 0 {
		return 0
	}
	acc := atom
	if x <= lo {
		return acc
	}
	upper := x
	if upper > hi {
		upper = hi
	}
	// Composite Simpson integration of the continuous density on [lo, upper].
	const steps = 2000
	f := func(t float64) float64 { return MarchenkoPasturDensity(t, c, sigma2) }
	acc += simpson(f, lo, upper, steps)
	if acc > 1 {
		acc = 1
	}
	return acc
}

// MarchenkoPasturMoment returns the k-th moment of the Marchenko-Pastur law
// with aspect ratio c and variance sigma2. It uses the Narayana expansion
// m_k = sigma2^k * sum_{r=0}^{k-1} 1/(r+1) * C(k,r) * C(k-1,r) * c^r.
func MarchenkoPasturMoment(k int, c, sigma2 float64) float64 {
	if k < 0 {
		return math.NaN()
	}
	if k == 0 {
		return 1
	}
	var s float64
	for r := 0; r <= k-1; r++ {
		s += BinomialCoefficient(k, r) * BinomialCoefficient(k-1, r) * math.Pow(c, float64(r)) / float64(r+1)
	}
	return math.Pow(sigma2, float64(k)) * s
}

// MarchenkoPasturMean returns the mean sigma2 of the Marchenko-Pastur law.
func MarchenkoPasturMean(c, sigma2 float64) float64 { return sigma2 }

// MarchenkoPasturVariance returns the variance c*sigma2^2 of the
// Marchenko-Pastur law.
func MarchenkoPasturVariance(c, sigma2 float64) float64 { return c * sigma2 * sigma2 }

// MarchenkoPasturStieltjes returns the Stieltjes transform of the standard
// (sigma2 = 1) Marchenko-Pastur law with ratio c at complex argument z, using
// the branch that behaves like -1/z at infinity.
func MarchenkoPasturStieltjes(z complex128, c float64) complex128 {
	cc := complex(c, 0)
	one := complex(1, 0)
	disc := cmplx.Sqrt((z-one-cc)*(z-one-cc) - 4*cc)
	if imag(z) > 0 && imag(disc) < 0 || imag(z) < 0 && imag(disc) > 0 {
		disc = -disc
	}
	return (-(z - one + cc) + disc) / (2 * cc * z)
}

// simpson approximates the integral of f over [a, b] with n subintervals
// (n is rounded up to an even number) using the composite Simpson rule.
func simpson(f func(float64) float64, a, b float64, n int) float64 {
	if b <= a {
		return 0
	}
	if n < 2 {
		n = 2
	}
	if n%2 == 1 {
		n++
	}
	h := (b - a) / float64(n)
	sum := f(a) + f(b)
	for i := 1; i < n; i++ {
		x := a + float64(i)*h
		if i%2 == 1 {
			sum += 4 * f(x)
		} else {
			sum += 2 * f(x)
		}
	}
	return sum * h / 3
}
