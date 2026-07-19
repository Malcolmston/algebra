package approxtheory

// PadeResult holds a Pade approximant P(x)/Q(x) with the numerator and
// denominator given as monomial coefficients in ascending order. By
// convention Q[0] = 1.
type PadeResult struct {
	Num []float64 // numerator coefficients, degree M
	Den []float64 // denominator coefficients, degree N, Den[0] == 1
}

// Eval evaluates the Pade approximant at x.
func (p *PadeResult) Eval(x float64) float64 {
	return Polyval(p.Num, x) / Polyval(p.Den, x)
}

// PadeApprox computes the [m/n] Pade approximant of the power series whose
// Taylor coefficients (about 0) are given in taylor, taylor[k] being the
// coefficient of x**k. It requires len(taylor) >= m+n+1. The numerator has
// degree m and the denominator degree n with Q(0)=1.
func PadeApprox(taylor []float64, m, n int) (*PadeResult, error) {
	if m < 0 || n < 0 {
		return nil, ErrDimensionMismatch
	}
	if len(taylor) < m+n+1 {
		return nil, ErrDimensionMismatch
	}
	a := taylor
	// Solve for denominator coefficients q_1..q_n from the n equations
	//   sum_{j=1}^{n} q_j a_{m+i-j} = -a_{m+i},   i = 1..n.
	den := make([]float64, n+1)
	den[0] = 1
	if n > 0 {
		A := make([][]float64, n)
		rhs := make([]float64, n)
		for i := 1; i <= n; i++ {
			row := make([]float64, n)
			for j := 1; j <= n; j++ {
				idx := m + i - j
				if idx >= 0 && idx < len(a) {
					row[j-1] = a[idx]
				}
			}
			A[i-1] = row
			rhs[i-1] = -aAt(a, m+i)
		}
		q, err := solveLinear(A, rhs)
		if err != nil {
			return nil, err
		}
		copy(den[1:], q)
	}
	// Numerator: p_i = sum_{j=0}^{i} q_j a_{i-j}, i=0..m.
	num := make([]float64, m+1)
	for i := 0; i <= m; i++ {
		var s float64
		for j := 0; j <= n && j <= i; j++ {
			s += den[j] * aAt(a, i-j)
		}
		num[i] = s
	}
	return &PadeResult{Num: num, Den: den}, nil
}

// aAt returns taylor[i] or 0 when the index is out of range.
func aAt(a []float64, i int) float64 {
	if i < 0 || i >= len(a) {
		return 0
	}
	return a[i]
}

// TaylorExp returns the first n+1 Taylor coefficients of exp about 0.
func TaylorExp(n int) []float64 {
	out := make([]float64, n+1)
	fact := 1.0
	for k := 0; k <= n; k++ {
		if k > 0 {
			fact *= float64(k)
		}
		out[k] = 1.0 / fact
	}
	return out
}

// TaylorSin returns the first n+1 Taylor coefficients of sin about 0.
func TaylorSin(n int) []float64 {
	out := make([]float64, n+1)
	fact := 1.0
	for k := 0; k <= n; k++ {
		if k > 0 {
			fact *= float64(k)
		}
		switch k % 4 {
		case 1:
			out[k] = 1.0 / fact
		case 3:
			out[k] = -1.0 / fact
		}
	}
	return out
}

// TaylorCos returns the first n+1 Taylor coefficients of cos about 0.
func TaylorCos(n int) []float64 {
	out := make([]float64, n+1)
	fact := 1.0
	for k := 0; k <= n; k++ {
		if k > 0 {
			fact *= float64(k)
		}
		switch k % 4 {
		case 0:
			out[k] = 1.0 / fact
		case 2:
			out[k] = -1.0 / fact
		}
	}
	return out
}

// TaylorLog1p returns the first n+1 Taylor coefficients of log(1+x) about 0.
func TaylorLog1p(n int) []float64 {
	out := make([]float64, n+1)
	for k := 1; k <= n; k++ {
		sign := 1.0
		if k%2 == 0 {
			sign = -1.0
		}
		out[k] = sign / float64(k)
	}
	return out
}

// PadeExp returns the [m/n] Pade approximant of exp about 0, a convenience
// wrapper around PadeApprox with the exponential Taylor series.
func PadeExp(m, n int) *PadeResult {
	res, _ := PadeApprox(TaylorExp(m+n), m, n)
	return res
}
