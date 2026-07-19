package randommatrix

import (
	"math"
	"math/big"
)

// FactorialBig returns n! as an exact big integer for n >= 0 and nil for n < 0.
func FactorialBig(n int) *big.Int {
	if n < 0 {
		return nil
	}
	r := big.NewInt(1)
	for i := 2; i <= n; i++ {
		r.Mul(r, big.NewInt(int64(i)))
	}
	return r
}

// BinomialBig returns the exact binomial coefficient C(n, k) as a big integer.
// It returns zero when k < 0 or k > n.
func BinomialBig(n, k int) *big.Int {
	if k < 0 || k > n || n < 0 {
		return big.NewInt(0)
	}
	r := big.NewInt(1)
	r.Binomial(int64(n), int64(k))
	return r
}

// CatalanBig returns the n-th Catalan number C_n = C(2n, n)/(n+1) as an exact
// big integer for n >= 0.
func CatalanBig(n int) *big.Int {
	if n < 0 {
		return big.NewInt(0)
	}
	c := BinomialBig(2*n, n)
	c.Div(c, big.NewInt(int64(n+1)))
	return c
}

// BinomialCoefficient returns C(n, k) as a float64, computed exactly with big
// integers and then converted. It returns zero for k < 0 or k > n.
func BinomialCoefficient(n, k int) float64 {
	b := BinomialBig(n, k)
	f := new(big.Float).SetInt(b)
	v, _ := f.Float64()
	return v
}

// CatalanNumber returns the n-th Catalan number as a float64.
func CatalanNumber(n int) float64 {
	c := CatalanBig(n)
	f := new(big.Float).SetInt(c)
	v, _ := f.Float64()
	return v
}

// NonCrossingPartitionsCount returns the number of non-crossing partitions of an
// n-element set, which equals the n-th Catalan number.
func NonCrossingPartitionsCount(n int) float64 { return CatalanNumber(n) }

// SpectralMoment returns the k-th empirical spectral moment (1/n) sum lambda_i^k
// of the eigenvalue sample eigs.
func SpectralMoment(eigs []float64, k int) float64 {
	if len(eigs) == 0 {
		return math.NaN()
	}
	var s float64
	for _, l := range eigs {
		s += math.Pow(l, float64(k))
	}
	return s / float64(len(eigs))
}

// SpectralMoments returns the empirical spectral moments of orders 1..kmax.
func SpectralMoments(eigs []float64, kmax int) []float64 {
	out := make([]float64, kmax)
	for k := 1; k <= kmax; k++ {
		out[k-1] = SpectralMoment(eigs, k)
	}
	return out
}

// SpectralMean returns the mean of the eigenvalue sample.
func SpectralMean(eigs []float64) float64 { return SpectralMoment(eigs, 1) }

// SpectralVariance returns the population variance of the eigenvalue sample.
func SpectralVariance(eigs []float64) float64 {
	if len(eigs) == 0 {
		return math.NaN()
	}
	mean := SpectralMean(eigs)
	var s float64
	for _, l := range eigs {
		d := l - mean
		s += d * d
	}
	return s / float64(len(eigs))
}

// CentralSpectralMoment returns the k-th central moment (1/n) sum (lambda_i -
// mean)^k of the eigenvalue sample.
func CentralSpectralMoment(eigs []float64, k int) float64 {
	if len(eigs) == 0 {
		return math.NaN()
	}
	mean := SpectralMean(eigs)
	var s float64
	for _, l := range eigs {
		s += math.Pow(l-mean, float64(k))
	}
	return s / float64(len(eigs))
}

// EmpiricalStieltjes returns the Stieltjes transform (1/n) sum 1/(lambda_i - z)
// of the empirical spectral distribution at complex argument z.
func EmpiricalStieltjes(eigs []float64, z complex128) complex128 {
	if len(eigs) == 0 {
		return 0
	}
	var s complex128
	for _, l := range eigs {
		s += 1 / (complex(l, 0) - z)
	}
	return s / complex(float64(len(eigs)), 0)
}

// EmpiricalCauchy returns the Cauchy transform (1/n) sum 1/(z - lambda_i), the
// negative of the Stieltjes transform.
func EmpiricalCauchy(eigs []float64, z complex128) complex128 {
	return -EmpiricalStieltjes(eigs, z)
}

// MomentsFromFreeCumulants returns the moments m_1..m_N implied by the free
// cumulants kappa_1..kappa_N via the free moment-cumulant relation
// m_n = sum over non-crossing partitions of products of cumulants.
func MomentsFromFreeCumulants(kappa []float64) []float64 {
	N := len(kappa)
	m := make([]float64, N+1)
	m[0] = 1
	for n := 1; n <= N; n++ {
		pow := make([]float64, n)
		copy(pow, m[:n])
		var sum float64
		for s := 1; s <= n; s++ {
			idx := n - s
			if idx >= 0 && idx < len(pow) {
				sum += kappa[s-1] * pow[idx]
			}
			if s < n {
				pow = convolveTrunc(pow, m[:n], n)
			}
		}
		m[n] = sum
	}
	out := make([]float64, N)
	copy(out, m[1:])
	return out
}

// FreeCumulantsFromMoments returns the free cumulants kappa_1..kappa_N implied
// by the moments m_1..m_N by inverting the free moment-cumulant relation.
func FreeCumulantsFromMoments(moments []float64) []float64 {
	N := len(moments)
	m := make([]float64, N+1)
	m[0] = 1
	copy(m[1:], moments)
	kappa := make([]float64, N)
	for n := 1; n <= N; n++ {
		pow := make([]float64, n)
		copy(pow, m[:n])
		var sum float64
		for s := 1; s <= n-1; s++ {
			idx := n - s
			if idx >= 0 && idx < len(pow) {
				sum += kappa[s-1] * pow[idx]
			}
			pow = convolveTrunc(pow, m[:n], n)
		}
		kappa[n-1] = m[n] - sum
	}
	return kappa
}

// convolveTrunc returns the polynomial product of a and b truncated to length n.
func convolveTrunc(a, b []float64, n int) []float64 {
	out := make([]float64, n)
	for i := 0; i < len(a) && i < n; i++ {
		if a[i] == 0 {
			continue
		}
		for j := 0; j < len(b) && i+j < n; j++ {
			out[i+j] += a[i] * b[j]
		}
	}
	return out
}

// FreeConvolutionMoments returns the moments of the free additive convolution of
// two distributions given their moment sequences. Free cumulants add, so the
// result is obtained by adding cumulants and converting back. The shorter of the
// two inputs sets the number of returned moments.
func FreeConvolutionMoments(momentsA, momentsB []float64) []float64 {
	n := len(momentsA)
	if len(momentsB) < n {
		n = len(momentsB)
	}
	ka := FreeCumulantsFromMoments(momentsA[:n])
	kb := FreeCumulantsFromMoments(momentsB[:n])
	sum := make([]float64, n)
	for i := range sum {
		sum[i] = ka[i] + kb[i]
	}
	return MomentsFromFreeCumulants(sum)
}

// RTransformSeries returns the coefficients of the R-transform power series
// R(z) = sum_{n>=1} kappa_n z^(n-1) from the free cumulants kappa_1..kappa_N.
// The returned slice holds the coefficient of z^0, z^1, ... .
func RTransformSeries(kappa []float64) []float64 {
	out := make([]float64, len(kappa))
	copy(out, kappa)
	return out
}

// SemicircleFreeCumulants returns the free cumulants of the semicircle law of
// the given variance up to order N. Only the second cumulant is non-zero and
// equals the variance.
func SemicircleFreeCumulants(variance float64, N int) []float64 {
	k := make([]float64, N)
	if N >= 2 {
		k[1] = variance
	}
	return k
}

// FreeConvolutionSemicircleVariance returns the variance of the free additive
// convolution of two semicircle laws, which is the sum of their variances.
func FreeConvolutionSemicircleVariance(var1, var2 float64) float64 {
	return var1 + var2
}

// FreeConvolutionSemicircleRadius returns the radius of the semicircle law
// obtained as the free additive convolution of two semicircle laws of the given
// radii, namely sqrt(R1^2 + R2^2).
func FreeConvolutionSemicircleRadius(r1, r2 float64) float64 {
	return math.Sqrt(r1*r1 + r2*r2)
}
