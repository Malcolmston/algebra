package powerseries

import "math"

// Add returns the sum s + other. The result has the precision of the longer
// operand; the shorter is treated as zero-padded.
func (s Series) Add(other Series) Series {
	n := powerseriesMaxLen(len(s.coeffs), len(other.coeffs))
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = s.Coeff(i) + other.Coeff(i)
	}
	return Series{coeffs: out}
}

// Sub returns the difference s − other. The result has the precision of the
// longer operand.
func (s Series) Sub(other Series) Series {
	n := powerseriesMaxLen(len(s.coeffs), len(other.coeffs))
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = s.Coeff(i) - other.Coeff(i)
	}
	return Series{coeffs: out}
}

// Neg returns −s.
func (s Series) Neg() Series {
	out := make([]float64, len(s.coeffs))
	for i, c := range s.coeffs {
		out[i] = -c
	}
	return Series{coeffs: out}
}

// Scale returns c·s, multiplying every coefficient by the scalar c.
func (s Series) Scale(c float64) Series {
	out := make([]float64, len(s.coeffs))
	for i, v := range s.coeffs {
		out[i] = c * v
	}
	return Series{coeffs: out}
}

// Mul returns the Cauchy product s · other, truncated to the precision of the
// longer operand. Coefficients beyond the stored range of either operand are
// taken to be zero.
func (s Series) Mul(other Series) Series {
	n := powerseriesMaxLen(len(s.coeffs), len(other.coeffs))
	out := make([]float64, n)
	for i, a := range s.coeffs {
		if a == 0 {
			continue
		}
		for j, b := range other.coeffs {
			if i+j < n {
				out[i+j] += a * b
			}
		}
	}
	return Series{coeffs: out}
}

// Hadamard returns the coefficient-wise (Hadamard) product of s and other, whose
// degree-i coefficient is s.Coeff(i)·other.Coeff(i). The result has the
// precision of the longer operand.
func (s Series) Hadamard(other Series) Series {
	n := powerseriesMaxLen(len(s.coeffs), len(other.coeffs))
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = s.Coeff(i) * other.Coeff(i)
	}
	return Series{coeffs: out}
}

// Shift returns x^k · s, moving every coefficient up by k degrees and dropping
// terms pushed beyond the current precision. A negative k divides by x^|k|,
// discarding the low-order coefficients (which must be understood as an exact
// division only when those coefficients vanish). The result keeps the receiver's
// precision.
func (s Series) Shift(k int) Series {
	n := len(s.coeffs)
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		j := i - k
		if j >= 0 && j < n {
			out[i] = s.coeffs[j]
		}
	}
	return Series{coeffs: out}
}

// Pow returns s raised to the non-negative integer power k, using
// exponentiation by squaring. Pow(0) is the constant series 1 at the receiver's
// precision. It panics if k is negative; use Inverse together with Pow, or
// PowReal, for negative or fractional powers.
func (s Series) Pow(k int) Series {
	if k < 0 {
		panic("powerseries: Pow requires a non-negative exponent")
	}
	result := One(len(s.coeffs))
	base := s.Clone()
	for k > 0 {
		if k&1 == 1 {
			result = result.Mul(base)
		}
		k >>= 1
		if k > 0 {
			base = base.Mul(base)
		}
	}
	return result
}

// Derivative returns the formal derivative d/dx of the series. Differentiating a
// length-n series yields a length-(n−1) series (length one for a constant),
// reflecting the coefficient that is genuinely known.
func (s Series) Derivative() Series {
	n := len(s.coeffs)
	if n <= 1 {
		return Series{coeffs: []float64{0}}
	}
	out := make([]float64, n-1)
	for i := 1; i < n; i++ {
		out[i-1] = float64(i) * s.coeffs[i]
	}
	return Series{coeffs: out}
}

// IntegralConst returns the antiderivative of the series whose constant term is
// c. Integrating a length-n series yields a length-(n+1) series.
func (s Series) IntegralConst(c float64) Series {
	n := len(s.coeffs)
	out := make([]float64, n+1)
	out[0] = c
	for i := 0; i < n; i++ {
		out[i+1] = s.coeffs[i] / float64(i+1)
	}
	return Series{coeffs: out}
}

// Integral returns the antiderivative of the series whose constant term is zero.
func (s Series) Integral() Series {
	return s.IntegralConst(0)
}

// Inverse returns the multiplicative inverse 1/s as a truncated power series at
// the receiver's precision. It panics if the constant term is zero, since only
// series with an invertible constant term have a power-series reciprocal.
func (s Series) Inverse() Series {
	n := len(s.coeffs)
	if s.coeffs[0] == 0 {
		panic("powerseries: Inverse requires a non-zero constant term")
	}
	b := make([]float64, n)
	b[0] = 1 / s.coeffs[0]
	for k := 1; k < n; k++ {
		var acc float64
		for i := 1; i <= k; i++ {
			acc += s.coeffs[i] * b[k-i]
		}
		b[k] = -acc / s.coeffs[0]
	}
	return Series{coeffs: b}
}

// Div returns the quotient s/other as a truncated power series. It panics if the
// constant term of other is zero.
func (s Series) Div(other Series) Series {
	return s.Mul(other.Inverse())
}

// Compose returns the composition s(other(x)), truncated to the precision of
// other. The inner series must have a zero constant term, other.Coeff(0) == 0,
// so that the composition is a well-defined formal power series; the method
// panics otherwise.
func (s Series) Compose(other Series) Series {
	if other.Coeff(0) != 0 {
		panic("powerseries: Compose requires the inner series to have zero constant term")
	}
	n := len(other.coeffs)
	top := len(s.coeffs)
	if top > n {
		top = n
	}
	result := Zero(n)
	for i := top - 1; i >= 0; i-- {
		result = result.Mul(other)
		result.coeffs[0] += s.coeffs[i]
	}
	return result
}

// Reversion returns the compositional inverse g of the series, the unique series
// with zero constant term satisfying s(g(x)) = g(s(x)) = x. The series must have
// a zero constant term and a non-zero linear term; the method panics otherwise.
// The result is obtained from the Lagrange inversion formula and has the
// receiver's precision.
func (s Series) Reversion() Series {
	n := len(s.coeffs)
	if s.Coeff(0) != 0 {
		panic("powerseries: Reversion requires a zero constant term")
	}
	if s.Coeff(1) == 0 {
		panic("powerseries: Reversion requires a non-zero linear term")
	}
	// Write s(x) = x·psi(x) with psi(0) = s.Coeff(1) != 0. The inverse g solves
	// g = t·phi(g) with phi = 1/psi, so Lagrange inversion applies with this phi.
	if n == 1 {
		return Zero(1)
	}
	psi := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		psi[i] = s.coeffs[i+1]
	}
	phi := Series{coeffs: psi}.Inverse()
	return powerseriesLagrange(phi, n)
}

// Exp returns exp(s) as a truncated power series. The constant term may be
// arbitrary: the result has constant term exp(s.Coeff(0)). The computation uses
// the recurrence n·b[n] = Σ k·a[k]·b[n−k] derived from b′ = a′·b.
func (s Series) Exp() Series {
	n := len(s.coeffs)
	b := make([]float64, n)
	b[0] = math.Exp(s.coeffs[0])
	for m := 1; m < n; m++ {
		var acc float64
		for k := 1; k <= m; k++ {
			acc += float64(k) * s.coeffs[k] * b[m-k]
		}
		b[m] = acc / float64(m)
	}
	return Series{coeffs: b}
}

// Log returns the natural logarithm log(s) as a truncated power series. The
// constant term must be strictly positive; the result has constant term
// log(s.Coeff(0)). The method panics if the constant term is not positive.
func (s Series) Log() Series {
	n := len(s.coeffs)
	a0 := s.coeffs[0]
	if a0 <= 0 {
		panic("powerseries: Log requires a positive constant term")
	}
	b := make([]float64, n)
	b[0] = math.Log(a0)
	for m := 1; m < n; m++ {
		acc := s.coeffs[m]
		for k := 1; k < m; k++ {
			acc -= s.coeffs[k] * float64(m-k) * b[m-k] / float64(m)
		}
		b[m] = acc / a0
	}
	return Series{coeffs: b}
}

// PowReal returns s raised to the real power p as a truncated power series. The
// constant term must be strictly positive; the result has constant term
// s.Coeff(0)^p. The method panics if the constant term is not positive.
func (s Series) PowReal(p float64) Series {
	n := len(s.coeffs)
	a0 := s.coeffs[0]
	if a0 <= 0 {
		panic("powerseries: PowReal requires a positive constant term")
	}
	b := make([]float64, n)
	b[0] = math.Pow(a0, p)
	for m := 1; m < n; m++ {
		acc := p * float64(m) * s.coeffs[m] * b[0]
		for k := 1; k < m; k++ {
			acc += (p*float64(k) - float64(m-k)) * s.coeffs[k] * b[m-k]
		}
		b[m] = acc / (a0 * float64(m))
	}
	return Series{coeffs: b}
}

// Sqrt returns the principal square root of the series as a truncated power
// series. It is defined for a strictly positive constant term and equals
// PowReal(0.5); the method panics if the constant term is not positive.
func (s Series) Sqrt() Series {
	return s.PowReal(0.5)
}

// Sin returns sin(s) as a truncated power series, valid for an arbitrary
// constant term. It is computed together with Cos from the coupled recurrence
// for s′ = c·a′, c′ = −s·a′.
func (s Series) Sin() Series {
	sin, _ := powerseriesSinCos(s)
	return sin
}

// Cos returns cos(s) as a truncated power series, valid for an arbitrary
// constant term.
func (s Series) Cos() Series {
	_, cos := powerseriesSinCos(s)
	return cos
}

// Tan returns tan(s) as a truncated power series, computed as Sin/Cos. It panics
// if the leading term of cos(s) vanishes, that is when cos(s.Coeff(0)) == 0.
func (s Series) Tan() Series {
	sin, cos := powerseriesSinCos(s)
	return sin.Div(cos)
}

// Sinh returns sinh(s) as a truncated power series, valid for an arbitrary
// constant term.
func (s Series) Sinh() Series {
	sinh, _ := powerseriesSinhCosh(s)
	return sinh
}

// Cosh returns cosh(s) as a truncated power series, valid for an arbitrary
// constant term.
func (s Series) Cosh() Series {
	_, cosh := powerseriesSinhCosh(s)
	return cosh
}

// Atan returns arctan(s) as a truncated power series, valid for an arbitrary
// constant term. It integrates a′/(1 + a²) using the identity (1 + a²)·b′ = a′.
func (s Series) Atan() Series {
	n := len(s.coeffs)
	// D = 1 + s².
	d := s.Mul(s)
	d.coeffs[0] += 1
	bp := make([]float64, n) // coefficients of b′
	for m := 0; m < n; m++ {
		var ap float64
		if m+1 < n {
			ap = float64(m+1) * s.coeffs[m+1]
		}
		acc := ap
		for k := 1; k <= m; k++ {
			acc -= d.coeffs[k] * bp[m-k]
		}
		bp[m] = acc / d.coeffs[0]
	}
	b := make([]float64, n)
	b[0] = math.Atan(s.coeffs[0])
	for m := 1; m < n; m++ {
		b[m] = bp[m-1] / float64(m)
	}
	return Series{coeffs: b}
}

// powerseriesSinCos returns sin(s) and cos(s) via the coupled recurrence
// n·sin[n] = Σ k·a[k]·cos[n−k], n·cos[n] = −Σ k·a[k]·sin[n−k].
func powerseriesSinCos(s Series) (Series, Series) {
	n := len(s.coeffs)
	sin := make([]float64, n)
	cos := make([]float64, n)
	sin[0] = math.Sin(s.coeffs[0])
	cos[0] = math.Cos(s.coeffs[0])
	for m := 1; m < n; m++ {
		var sAcc, cAcc float64
		for k := 1; k <= m; k++ {
			ka := float64(k) * s.coeffs[k]
			sAcc += ka * cos[m-k]
			cAcc += ka * sin[m-k]
		}
		sin[m] = sAcc / float64(m)
		cos[m] = -cAcc / float64(m)
	}
	return Series{coeffs: sin}, Series{coeffs: cos}
}

// powerseriesSinhCosh returns sinh(s) and cosh(s) via the coupled recurrence
// n·sinh[n] = Σ k·a[k]·cosh[n−k], n·cosh[n] = Σ k·a[k]·sinh[n−k].
func powerseriesSinhCosh(s Series) (Series, Series) {
	n := len(s.coeffs)
	sinh := make([]float64, n)
	cosh := make([]float64, n)
	sinh[0] = math.Sinh(s.coeffs[0])
	cosh[0] = math.Cosh(s.coeffs[0])
	for m := 1; m < n; m++ {
		var sAcc, cAcc float64
		for k := 1; k <= m; k++ {
			ka := float64(k) * s.coeffs[k]
			sAcc += ka * cosh[m-k]
			cAcc += ka * sinh[m-k]
		}
		sinh[m] = sAcc / float64(m)
		cosh[m] = cAcc / float64(m)
	}
	return Series{coeffs: sinh}, Series{coeffs: cosh}
}
