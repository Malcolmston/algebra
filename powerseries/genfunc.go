package powerseries

// powerseriesFactorial returns i! as a float64.
func powerseriesFactorial(i int) float64 {
	f := 1.0
	for k := 2; k <= i; k++ {
		f *= float64(k)
	}
	return f
}

// powerseriesBinom returns the binomial coefficient C(n, k) as a float64, using
// a product form that avoids forming the individual factorials.
func powerseriesBinom(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	if k > n-k {
		k = n - k
	}
	res := 1.0
	for i := 0; i < k; i++ {
		res = res * float64(n-i) / float64(i+1)
	}
	return res
}

// OGFtoEGF reinterprets the ordinary generating function s, whose degree-n
// coefficient is the sequence value a_n, as the exponential generating function
// of the same sequence, dividing the degree-n coefficient by n!. It is the
// Borel transform.
func OGFtoEGF(s Series) Series {
	out := make([]float64, len(s.coeffs))
	for i, c := range s.coeffs {
		out[i] = c / powerseriesFactorial(i)
	}
	return Series{coeffs: out}
}

// EGFtoOGF reinterprets the exponential generating function s, whose degree-n
// coefficient is a_n/n!, as the ordinary generating function of the same
// sequence, multiplying the degree-n coefficient by n!. It is the inverse Borel
// transform.
func EGFtoOGF(s Series) Series {
	out := make([]float64, len(s.coeffs))
	for i, c := range s.coeffs {
		out[i] = c * powerseriesFactorial(i)
	}
	return Series{coeffs: out}
}

// SequenceFromOGF returns the underlying sequence a_0, a_1, … of the series read
// as an ordinary generating function, that is a plain copy of the coefficients.
func (s Series) SequenceFromOGF() []float64 {
	return s.Coeffs()
}

// SequenceFromEGF returns the underlying sequence a_0, a_1, … of the series read
// as an exponential generating function, multiplying the degree-n coefficient by
// n!.
func (s Series) SequenceFromEGF() []float64 {
	out := make([]float64, len(s.coeffs))
	for i, c := range s.coeffs {
		out[i] = c * powerseriesFactorial(i)
	}
	return out
}

// OGFFromSequence returns the ordinary generating function whose degree-n
// coefficient is a[n].
func OGFFromSequence(a []float64) Series {
	return FromSlice(a)
}

// EGFFromSequence returns the exponential generating function whose degree-n
// coefficient is a[n]/n!.
func EGFFromSequence(a []float64) Series {
	out := make([]float64, len(a))
	for i, v := range a {
		out[i] = v / powerseriesFactorial(i)
	}
	if len(out) == 0 {
		out = []float64{0}
	}
	return Series{coeffs: out}
}

// GeometricGF returns the length-n ordinary generating function 1/(1−x) =
// Σ x^k, whose coefficients are all one.
func GeometricGF(n int) Series {
	s := Zero(n)
	for i := range s.coeffs {
		s.coeffs[i] = 1
	}
	return s
}

// GeometricRatioGF returns the length-n ordinary generating function
// 1/(1−r·x) = Σ r^k·x^k.
func GeometricRatioGF(r float64, n int) Series {
	s := Zero(n)
	p := 1.0
	for i := range s.coeffs {
		s.coeffs[i] = p
		p *= r
	}
	return s
}

// ExpGF returns the length-n exponential series e^x = Σ x^k/k!, the exponential
// generating function of the all-ones sequence.
func ExpGF(n int) Series {
	s := Zero(n)
	for i := range s.coeffs {
		s.coeffs[i] = 1 / powerseriesFactorial(i)
	}
	return s
}

// BinomialGF returns the length-n binomial series (1+x)^alpha, whose degree-k
// coefficient is the generalized binomial coefficient C(alpha, k). For a
// non-negative integer alpha it terminates as a polynomial.
func BinomialGF(alpha float64, n int) Series {
	s := Zero(n)
	c := 1.0
	for k := 0; k < n; k++ {
		s.coeffs[k] = c
		c = c * (alpha - float64(k)) / float64(k+1)
	}
	return s
}

// CatalanGF returns the length-n ordinary generating function of the Catalan
// numbers, C(x) = (1 − √(1−4x))/(2x), computed with the series square root so
// that its degree-k coefficient is the k-th Catalan number.
func CatalanGF(n int) Series {
	// u = 1 − 4x to precision n+1 so the square root delivers n+1 coefficients.
	u := Zero(n + 1)
	u.coeffs[0] = 1
	if n+1 > 1 {
		u.coeffs[1] = -4
	}
	root := u.Sqrt()
	// C(x) = (1 − root)/(2x): coefficient k equals −root[k+1]/2.
	s := Zero(n)
	for k := 0; k < n; k++ {
		s.coeffs[k] = -root.coeffs[k+1] / 2
	}
	return s
}

// FibonacciGF returns the length-n ordinary generating function of the Fibonacci
// numbers, x/(1 − x − x²), whose degree-k coefficient is F_k with F_0 = 0 and
// F_1 = 1.
func FibonacciGF(n int) Series {
	den := Zero(n)
	den.coeffs[0] = 1
	if n > 1 {
		den.coeffs[1] = -1
	}
	if n > 2 {
		den.coeffs[2] = -1
	}
	return den.Inverse().Shift(1)
}

// BernoulliGF returns the length-n exponential generating function of the
// Bernoulli numbers, x/(e^x − 1) = Σ B_k·x^k/k!, using the convention B_1 =
// −1/2. Its degree-k coefficient is B_k/k!.
func BernoulliGF(n int) Series {
	// d(x) = (e^x − 1)/x has coefficient 1/(k+1)! and constant term 1.
	d := Zero(n)
	for k := 0; k < n; k++ {
		d.coeffs[k] = 1 / powerseriesFactorial(k+1)
	}
	return d.Inverse()
}

// PartitionGF returns the length-n ordinary generating function of the integer
// partition numbers, Π_{k≥1} 1/(1−x^k), whose degree-m coefficient is p(m).
func PartitionGF(n int) Series {
	res := One(n)
	for k := 1; k < n; k++ {
		factor := Zero(n)
		for j := 0; j < n; j += k {
			factor.coeffs[j] = 1
		}
		res = res.Mul(factor)
	}
	return res
}

// BellGF returns the length-n exponential generating function of the Bell
// numbers, exp(e^x − 1), whose degree-k coefficient is B_k/k! where B_k counts
// the partitions of a k-element set.
func BellGF(n int) Series {
	inner := ExpGF(n)
	inner.coeffs[0] = 0 // e^x − 1
	return inner.Exp()
}

// DerangementGF returns the length-n exponential generating function of the
// derangement numbers, e^{−x}/(1−x), whose degree-k coefficient is D_k/k!.
func DerangementGF(n int) Series {
	expNeg := Zero(n)
	for k := 0; k < n; k++ {
		expNeg.coeffs[k] = powerseriesAltFactInv(k)
	}
	return expNeg.Mul(GeometricGF(n))
}

// powerseriesAltFactInv returns (−1)^k / k!.
func powerseriesAltFactInv(k int) float64 {
	v := 1 / powerseriesFactorial(k)
	if k%2 == 1 {
		return -v
	}
	return v
}

// HarmonicGF returns the length-n ordinary generating function of the harmonic
// numbers, −log(1−x)/(1−x), whose degree-k coefficient is H_k = Σ_{j=1}^{k} 1/j
// with H_0 = 0.
func HarmonicGF(n int) Series {
	oneMinusX := Zero(n)
	oneMinusX.coeffs[0] = 1
	if n > 1 {
		oneMinusX.coeffs[1] = -1
	}
	negLog := oneMinusX.Log().Neg()
	return negLog.Mul(GeometricGF(n))
}

// MotzkinGF returns the length-n ordinary generating function of the Motzkin
// numbers, (1 − x − √(1−2x−3x²))/(2x²), whose degree-k coefficient is the k-th
// Motzkin number.
func MotzkinGF(n int) Series {
	return FromSlice(MotzkinNumbers(n))
}

// CatalanNumbers returns the first n Catalan numbers C_0, C_1, … as float64
// values, where C_k = C(2k, k)/(k+1).
func CatalanNumbers(n int) []float64 {
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	out[0] = 1
	for k := 1; k < n; k++ {
		out[k] = out[k-1] * 2 * float64(2*k-1) / float64(k+1)
	}
	return out
}

// Catalan returns the n-th Catalan number as a float64.
func Catalan(n int) float64 {
	if n < 0 {
		return 0
	}
	return CatalanNumbers(n + 1)[n]
}

// BernoulliNumbers returns the first n Bernoulli numbers B_0, B_1, … as float64
// values, using the convention B_1 = −1/2. They are generated from the identity
// Σ_{j=0}^{m} C(m+1, j)·B_j = 0.
func BernoulliNumbers(n int) []float64 {
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	out[0] = 1
	for m := 1; m < n; m++ {
		var acc float64
		for j := 0; j < m; j++ {
			acc += powerseriesBinom(m+1, j) * out[j]
		}
		out[m] = -acc / float64(m+1)
	}
	return out
}

// Bernoulli returns the n-th Bernoulli number as a float64, using B_1 = −1/2.
func Bernoulli(n int) float64 {
	if n < 0 {
		return 0
	}
	return BernoulliNumbers(n + 1)[n]
}

// PartitionNumbers returns the first n integer-partition counts p(0), p(1), …
// as float64 values, computed with Euler's pentagonal-number recurrence.
func PartitionNumbers(n int) []float64 {
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	out[0] = 1
	for m := 1; m < n; m++ {
		var acc float64
		for k := 1; ; k++ {
			g1 := k * (3*k - 1) / 2
			g2 := k * (3*k + 1) / 2
			if g1 > m && g2 > m {
				break
			}
			sign := 1.0
			if k%2 == 0 {
				sign = -1.0
			}
			if g1 <= m {
				acc += sign * out[m-g1]
			}
			if g2 <= m {
				acc += sign * out[m-g2]
			}
		}
		out[m] = acc
	}
	return out
}

// PartitionCount returns p(n), the number of integer partitions of n, as a
// float64.
func PartitionCount(n int) float64 {
	if n < 0 {
		return 0
	}
	return PartitionNumbers(n + 1)[n]
}

// FibonacciNumbers returns the first n Fibonacci numbers F_0, F_1, … as float64
// values with F_0 = 0 and F_1 = 1.
func FibonacciNumbers(n int) []float64 {
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	out[0] = 0
	if n > 1 {
		out[1] = 1
	}
	for k := 2; k < n; k++ {
		out[k] = out[k-1] + out[k-2]
	}
	return out
}

// Fibonacci returns the n-th Fibonacci number as a float64 with F_0 = 0.
func Fibonacci(n int) float64 {
	if n < 0 {
		return 0
	}
	return FibonacciNumbers(n + 1)[n]
}

// BellNumbers returns the first n Bell numbers B_0, B_1, … as float64 values,
// computed with the Bell triangle.
func BellNumbers(n int) []float64 {
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	row := []float64{1}
	out[0] = 1
	for i := 1; i < n; i++ {
		next := make([]float64, i+1)
		next[0] = row[i-1]
		for j := 1; j <= i; j++ {
			next[j] = next[j-1] + row[j-1]
		}
		out[i] = next[0]
		row = next
	}
	return out
}

// Bell returns the n-th Bell number as a float64.
func Bell(n int) float64 {
	if n < 0 {
		return 0
	}
	return BellNumbers(n + 1)[n]
}

// DerangementNumbers returns the first n derangement counts D_0, D_1, … as
// float64 values, where D_k is the number of fixed-point-free permutations of k
// elements and D_k = (k−1)(D_{k−1} + D_{k−2}).
func DerangementNumbers(n int) []float64 {
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	out[0] = 1
	if n > 1 {
		out[1] = 0
	}
	for k := 2; k < n; k++ {
		out[k] = float64(k-1) * (out[k-1] + out[k-2])
	}
	return out
}

// Derangement returns the n-th derangement number as a float64.
func Derangement(n int) float64 {
	if n < 0 {
		return 0
	}
	return DerangementNumbers(n + 1)[n]
}

// MotzkinNumbers returns the first n Motzkin numbers M_0, M_1, … as float64
// values, using the recurrence (n+2)·M_n = (2n+1)·M_{n−1} + (3n−3)·M_{n−2}.
func MotzkinNumbers(n int) []float64 {
	out := make([]float64, n)
	if n == 0 {
		return out
	}
	out[0] = 1
	if n > 1 {
		out[1] = 1
	}
	for k := 2; k < n; k++ {
		out[k] = (float64(2*k+1)*out[k-1] + float64(3*k-3)*out[k-2]) / float64(k+2)
	}
	return out
}

// LagrangeInversion returns the length-n power series w(t) with zero constant
// term that solves the functional equation w = t·φ(w), where phi is the given
// series. The Lagrange inversion formula gives [t^m] w = (1/m)·[x^{m−1}] φ(x)^m.
// It panics if phi has a zero constant term.
func LagrangeInversion(phi Series, n int) Series {
	if phi.Coeff(0) == 0 {
		panic("powerseries: LagrangeInversion requires a non-zero constant term in phi")
	}
	return powerseriesLagrange(phi.Extend(n), n)
}

// powerseriesLagrange computes the Lagrange-inversion series w with w = t·φ(w)
// to length n from the coefficients of phi.
func powerseriesLagrange(phi Series, n int) Series {
	w := make([]float64, n)
	if n <= 1 {
		return Series{coeffs: w}
	}
	phipow := phi.Clone()
	for m := 1; m < n; m++ {
		w[m] = phipow.Coeff(m-1) / float64(m)
		if m+1 < n {
			phipow = phipow.Mul(phi)
		}
	}
	return Series{coeffs: w}
}

// LagrangeInversionApply returns the length-n series H(w(t)), where w is the
// solution of w = t·φ(w) with zero constant term. By the Lagrange–Bürmann
// formula [t^m] H(w) = (1/m)·[x^{m−1}](H′(x)·φ(x)^m) for m ≥ 1, and the constant
// term is H(0). It panics if phi has a zero constant term.
func LagrangeInversionApply(phi, h Series, n int) Series {
	if phi.Coeff(0) == 0 {
		panic("powerseries: LagrangeInversionApply requires a non-zero constant term in phi")
	}
	phi = phi.Extend(n)
	hp := h.Derivative()
	out := make([]float64, n)
	out[0] = h.Coeff(0)
	if n <= 1 {
		return Series{coeffs: out}
	}
	phipow := phi.Clone()
	for m := 1; m < n; m++ {
		prod := hp.Mul(phipow)
		out[m] = prod.Coeff(m-1) / float64(m)
		if m+1 < n {
			phipow = phipow.Mul(phi)
		}
	}
	return Series{coeffs: out}
}
