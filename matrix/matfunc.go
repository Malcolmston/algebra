package matrix

import (
	"fmt"
	"math"
)

// matrixPadeDegree is the fixed diagonal Padé degree used by the
// scaling-and-squaring matrix exponential. Combined with the scaling threshold
// matrixPadeScaleThreshold it yields a deterministic approximation that is
// accurate to double precision: after scaling, ‖A/2ˢ‖∞ < 1/2, and the
// truncation error of the [6/6] diagonal Padé approximant for an argument of
// that size is far below machine epsilon.
const matrixPadeDegree = 6

// matrixPadeScaleThreshold is the target infinity-norm bound for the scaled
// matrix A/2ˢ. The scaling power s is chosen as the smallest non-negative
// integer with ‖A‖∞/2ˢ < matrixPadeScaleThreshold.
const matrixPadeScaleThreshold = 0.5

// matrixPade6Coef holds the coefficients c₀…c₆ of the numerator polynomial
// p(x) = Σ cₖ xᵏ of the diagonal [6/6] Padé approximant to eˣ. The denominator
// is p(-x), so the same coefficients (with alternating signs on the odd powers)
// build both the numerator and denominator matrix polynomials, allowing the
// matrix powers Aᵏ to be computed once and reused across both series.
var matrixPade6Coef = [matrixPadeDegree + 1]float64{
	1.0,
	1.0 / 2.0,
	5.0 / 44.0,
	1.0 / 66.0,
	1.0 / 792.0,
	1.0 / 15840.0,
	1.0 / 665280.0,
}

// Exp computes the matrix exponential e^A of a real square matrix by the
// scaling-and-squaring method with a diagonal Padé approximant.
//
// The scaling power s is chosen as the smallest non-negative integer such that
// ‖A/2ˢ‖∞ < 1/2; the [6/6] diagonal Padé rational approximation to e^{A/2ˢ} is
// then evaluated by forming its numerator and denominator matrix polynomials
// (which share the powers of A/2ˢ) and solving the resulting linear system with
// Gaussian elimination and partial pivoting; finally the approximation is
// squared s times to recover e^A. The whole computation runs on flat row-major
// []float64 buffers: the matrix powers are built once and reused for both the
// numerator and denominator, and the s squarings are done in place with a single
// reused product buffer.
//
// The algorithm is deterministic (fixed Padé degree and scaling threshold) and
// returns the exact identity for the zero matrix. It returns [ErrNotSquare] for
// a non-square matrix and [ErrUnsupported] (wrapping the underlying evaluation
// error) when any entry is symbolic or otherwise not numeric. The result is
// returned as a matrix of inexact float literals via [FromFloats].
func (m *Matrix) Exp() (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	f, err := m.Floats()
	if err != nil {
		return nil, fmt.Errorf("%w: %v", ErrUnsupported, err)
	}
	n := m.rows
	a := make([]float64, n*n)
	for i := 0; i < n; i++ {
		row := f[i]
		for j := 0; j < n; j++ {
			a[i*n+j] = row[j]
		}
	}
	res := matrixExpFlat(a, n)
	return matrixFlatToMatrix(res, n, n), nil
}

// ExpScaled returns e^{tA} for the scalar t, the flow map of the linear ODE
// x' = Ax evaluated at time t. It is a convenience wrapper that scales the
// numeric entries of the matrix by t before applying the same
// scaling-and-squaring Padé algorithm as [Matrix.Exp]; in particular
// ExpScaled(1) equals Exp and ExpScaled(0) is the identity.
//
// Like [Matrix.Exp] it is deterministic, returns [ErrNotSquare] for a non-square
// matrix and [ErrUnsupported] (wrapping the underlying evaluation error) when any
// entry is symbolic or otherwise not numeric, and returns the result as a matrix
// of inexact float literals via [FromFloats].
func (m *Matrix) ExpScaled(t float64) (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	f, err := m.Floats()
	if err != nil {
		return nil, fmt.Errorf("%w: %v", ErrUnsupported, err)
	}
	n := m.rows
	a := make([]float64, n*n)
	for i := 0; i < n; i++ {
		row := f[i]
		for j := 0; j < n; j++ {
			a[i*n+j] = t * row[j]
		}
	}
	res := matrixExpFlat(a, n)
	return matrixFlatToMatrix(res, n, n), nil
}

// matrixExpFlat computes e^A for the n×n row-major buffer a and returns the
// result as a freshly allocated n×n row-major buffer. The zero matrix (and the
// degenerate 0×0 case) yields the exact identity. The algorithm is the
// scaling-and-squaring diagonal Padé method described on [Matrix.Exp].
func matrixExpFlat(a []float64, n int) []float64 {
	if n == 0 {
		return []float64{}
	}
	normInf := matrixInfNormFlat(a, n)
	if normInf == 0 {
		return matrixIdentityFlat(n)
	}

	// Choose the scaling power s so that ‖A/2ˢ‖∞ < matrixPadeScaleThreshold.
	s := 0
	for scaled := normInf; scaled >= matrixPadeScaleThreshold; scaled /= 2 {
		s++
	}
	// b = A / 2ˢ, computed with math.Ldexp to avoid overflow for large s.
	scale := math.Ldexp(1, -s)
	b := make([]float64, n*n)
	for i := range a {
		b[i] = a[i] * scale
	}

	// Powers of b: powers[k] = bᵏ, built once and reused for both the numerator
	// and denominator Padé series.
	powers := make([][]float64, matrixPadeDegree+1)
	powers[0] = matrixIdentityFlat(n)
	powers[1] = b
	for k := 2; k <= matrixPadeDegree; k++ {
		powers[k] = make([]float64, n*n)
		matrixMatMulFlat(powers[k], powers[k-1], b, n)
	}

	// Numerator num = Σ cₖ bᵏ; denominator den = Σ cₖ (-1)ᵏ bᵏ = p(-b).
	num := make([]float64, n*n)
	den := make([]float64, n*n)
	for k := 0; k <= matrixPadeDegree; k++ {
		ck := matrixPade6Coef[k]
		dk := ck
		if k%2 == 1 {
			dk = -ck
		}
		pk := powers[k]
		for idx := 0; idx < n*n; idx++ {
			num[idx] += ck * pk[idx]
			den[idx] += dk * pk[idx]
		}
	}

	// Solve den · X = num, giving X ≈ e^{A/2ˢ}.
	cur := matrixLUSolveFlat(den, num, n)

	// Square s times in place, ping-ponging between cur and a single reused
	// product buffer.
	tmp := make([]float64, n*n)
	for i := 0; i < s; i++ {
		matrixMatMulFlat(tmp, cur, cur, n)
		cur, tmp = tmp, cur
	}
	return cur
}

// matrixInfNormFlat returns the infinity norm (maximum absolute row sum) of the
// n×n row-major buffer a.
func matrixInfNormFlat(a []float64, n int) float64 {
	var norm float64
	for i := 0; i < n; i++ {
		var rowSum float64
		base := i * n
		for j := 0; j < n; j++ {
			rowSum += math.Abs(a[base+j])
		}
		if rowSum > norm {
			norm = rowSum
		}
	}
	return norm
}

// matrixIdentityFlat returns a freshly allocated n×n row-major identity buffer.
func matrixIdentityFlat(n int) []float64 {
	out := make([]float64, n*n)
	for i := 0; i < n; i++ {
		out[i*n+i] = 1
	}
	return out
}

// matrixMatMulFlat writes the product x·y of two n×n row-major buffers into dst,
// which must not alias x or y. dst is fully overwritten.
func matrixMatMulFlat(dst, x, y []float64, n int) {
	for i := 0; i < n; i++ {
		irow := i * n
		for j := 0; j < n; j++ {
			dst[irow+j] = 0
		}
		for k := 0; k < n; k++ {
			xik := x[irow+k]
			if xik == 0 {
				continue
			}
			krow := k * n
			for j := 0; j < n; j++ {
				dst[irow+j] += xik * y[krow+j]
			}
		}
	}
}

// matrixLUSolveFlat solves the linear system d·X = b for the n×n row-major
// buffers d and b, returning X as a freshly allocated n×n row-major buffer. It
// uses Gaussian elimination with partial pivoting, carrying the n right-hand-side
// columns of b through the elimination and then back-substituting. The inputs d
// and b are not modified. For the well-conditioned Padé denominator produced by
// the scaling-and-squaring method the system is always non-singular.
func matrixLUSolveFlat(d, b []float64, n int) []float64 {
	lu := make([]float64, n*n)
	copy(lu, d)
	x := make([]float64, n*n)
	copy(x, b)

	for k := 0; k < n; k++ {
		// Partial pivot: largest magnitude entry in column k at or below row k.
		p := k
		maxAbs := math.Abs(lu[k*n+k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(lu[i*n+k]); v > maxAbs {
				maxAbs = v
				p = i
			}
		}
		if p != k {
			for j := 0; j < n; j++ {
				lu[k*n+j], lu[p*n+j] = lu[p*n+j], lu[k*n+j]
				x[k*n+j], x[p*n+j] = x[p*n+j], x[k*n+j]
			}
		}
		pivot := lu[k*n+k]
		for i := k + 1; i < n; i++ {
			factor := lu[i*n+k] / pivot
			if factor == 0 {
				continue
			}
			for j := k; j < n; j++ {
				lu[i*n+j] -= factor * lu[k*n+j]
			}
			for j := 0; j < n; j++ {
				x[i*n+j] -= factor * x[k*n+j]
			}
		}
	}

	// Back-substitution: lu is now upper triangular.
	for col := 0; col < n; col++ {
		for i := n - 1; i >= 0; i-- {
			sum := x[i*n+col]
			for j := i + 1; j < n; j++ {
				sum -= lu[i*n+j] * x[j*n+col]
			}
			x[i*n+col] = sum / lu[i*n+i]
		}
	}
	return x
}
