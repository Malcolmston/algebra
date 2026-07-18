package matrix

import (
	"errors"
	"fmt"
	"math"

	"github.com/malcolmston/algebra"
)

// ErrNotPositiveDefinite reports that a matrix passed to [Matrix.Cholesky] is
// not symmetric positive definite, so no real lower-triangular Cholesky factor
// exists.
var ErrNotPositiveDefinite = errors.New("matrix: matrix is not positive definite")

// matrixEps is the double-precision machine epsilon, used to scale the
// pivot/singularity tolerances relative to the magnitude of the input.
const matrixEps = 2.220446049250313e-16

// matrixMaxAbs returns the largest absolute value in a, or 0 for an empty
// slice. It is used to size tolerances relative to the matrix scale.
func matrixMaxAbs(a []float64) float64 {
	m := 0.0
	for _, v := range a {
		if av := math.Abs(v); av > m {
			m = av
		}
	}
	return m
}

// matrixFloatsFlat evaluates every entry of m to a float64 and returns them in
// a fresh row-major flat slice together with the shape. It returns
// [ErrUnsupported] wrapping the underlying evaluation error when any entry
// contains a free symbol (or otherwise cannot be evaluated numerically), which
// is what makes the dense factorizations a numeric-only fast path.
func matrixFloatsFlat(m *Matrix) (flat []float64, rows, cols int, err error) {
	ff, e := m.Floats()
	if e != nil {
		return nil, 0, 0, fmt.Errorf("%w: %v", ErrUnsupported, e)
	}
	rows, cols = m.Rows(), m.Cols()
	flat = make([]float64, rows*cols)
	for i := 0; i < rows; i++ {
		copy(flat[i*cols:i*cols+cols], ff[i])
	}
	return flat, rows, cols, nil
}

// matrixReshape copies a row-major flat slice into a fresh [][]float64 of the
// given shape, suitable for [FromFloats].
func matrixReshape(flat []float64, rows, cols int) [][]float64 {
	out := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		out[i] = make([]float64, cols)
		copy(out[i], flat[i*cols:i*cols+cols])
	}
	return out
}

// matrixLUDecompose factors the n×n row-major matrix a in place using Doolittle
// elimination with partial pivoting, producing P*A = L*U. On return the strict
// lower triangle of a holds the unit-lower factor L's multipliers and the upper
// triangle (including the diagonal) holds U. perm is the row permutation, where
// perm[i] is the original row now occupying row i, sign is det(P) (±1), and
// singular reports whether any pivot fell at or below the scale-relative
// tolerance (in which case that column's elimination is skipped and U carries a
// vanishing diagonal entry).
//
// The permutation and a single scratch row are allocated once and reused across
// the elimination; the k-outer / i,j-inner traversal keeps every inner-loop
// access unit-stride in the flat buffer.
func matrixLUDecompose(a []float64, n int) (perm []int, sign int, singular bool) {
	perm = make([]int, n)
	for i := range perm {
		perm[i] = i
	}
	sign = 1
	scratch := make([]float64, n)

	tol := matrixMaxAbs(a) * float64(n) * matrixEps
	if tol == 0 {
		tol = matrixEps
	}

	for k := 0; k < n; k++ {
		// Deterministic pivot: largest magnitude in the column, ties broken by
		// the lowest row index (a strict '>' keeps the first maximum).
		p := k
		maxv := math.Abs(a[k*n+k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(a[i*n+k]); v > maxv {
				maxv = v
				p = i
			}
		}
		if p != k {
			copy(scratch, a[k*n:k*n+n])
			copy(a[k*n:k*n+n], a[p*n:p*n+n])
			copy(a[p*n:p*n+n], scratch)
			perm[k], perm[p] = perm[p], perm[k]
			sign = -sign
		}

		piv := a[k*n+k]
		if math.Abs(piv) <= tol {
			singular = true
			continue
		}
		for i := k + 1; i < n; i++ {
			bi := i * n
			f := a[bi+k] / piv
			a[bi+k] = f
			bk := k * n
			for j := k + 1; j < n; j++ {
				a[bi+j] -= f * a[bk+j]
			}
		}
	}
	return perm, sign, singular
}

// LU computes the partial-pivot LU factorization P*A = L*U of the (numeric)
// square matrix m. L is unit-lower-triangular, U is upper-triangular, P is a
// permutation matrix, and sign is det(P) (+1 or -1). The entries of m are
// evaluated to float64, so LU returns [ErrUnsupported] (wrapping the evaluation
// error) if any entry has a free symbol, and [ErrNotSquare] if m is not square.
// The factorization is always returned when the input is numeric; a singular
// matrix simply yields a U with a vanishing diagonal entry rather than an error.
func (m *Matrix) LU() (l, u, p *Matrix, sign int, err error) {
	if !m.IsSquare() {
		return nil, nil, nil, 0, ErrNotSquare
	}
	a, n, _, err := matrixFloatsFlat(m)
	if err != nil {
		return nil, nil, nil, 0, err
	}
	perm, sgn, _ := matrixLUDecompose(a, n)

	lf := make([][]float64, n)
	uf := make([][]float64, n)
	pf := make([][]float64, n)
	for i := 0; i < n; i++ {
		lf[i] = make([]float64, n)
		uf[i] = make([]float64, n)
		pf[i] = make([]float64, n)
		lf[i][i] = 1
		pf[i][perm[i]] = 1
		for j := 0; j < n; j++ {
			switch {
			case j < i:
				lf[i][j] = a[i*n+j]
			default:
				uf[i][j] = a[i*n+j]
			}
		}
	}
	return FromFloats(lf), FromFloats(uf), FromFloats(pf), sgn, nil
}

// DetLU returns the determinant of the (numeric) square matrix m as an
// algebra.Flt, computed as sign·∏diag(U) from the partial-pivot LU
// factorization. This is O(n^3) and is far faster than the O(n!) cofactor
// expansion of [Matrix.Det] for n ≥ 4; it applies only to numeric input. It
// returns [ErrNotSquare] for a non-square matrix and [ErrUnsupported] (wrapping
// the evaluation error) if any entry has a free symbol. A singular matrix
// yields an exact 0.
func (m *Matrix) DetLU() (algebra.Expr, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	a, n, _, err := matrixFloatsFlat(m)
	if err != nil {
		return nil, err
	}
	_, sign, singular := matrixLUDecompose(a, n)
	if singular {
		return algebra.Flt(0), nil
	}
	det := float64(sign)
	for k := 0; k < n; k++ {
		det *= a[k*n+k]
	}
	return algebra.Flt(det), nil
}

// SolveLU solves the square linear system m·x = b for the (numeric) matrix m
// using its partial-pivot LU factorization followed by forward and back
// substitution. It returns [ErrNotSquare] if m is not square, [ErrDimension] if
// len(b) differs from the matrix order, [ErrUnsupported] (wrapping the
// evaluation error) if any entry of m or b has a free symbol, and [ErrSingular]
// if a zero pivot is encountered. The result is returned as a vector of
// algebra.Flt components. Substitution reuses a single solution buffer and
// performs no per-element allocation.
func (m *Matrix) SolveLU(b *Vector) (*Vector, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	n := m.Rows()
	if b.Len() != n {
		return nil, ErrDimension
	}
	a, _, _, err := matrixFloatsFlat(m)
	if err != nil {
		return nil, err
	}
	bf := make([]float64, n)
	for i := 0; i < n; i++ {
		f, e := algebra.Evalf(b.At(i))
		if e != nil {
			return nil, fmt.Errorf("%w: %v", ErrUnsupported, e)
		}
		bf[i] = f
	}

	perm, _, singular := matrixLUDecompose(a, n)
	if singular {
		return nil, ErrSingular
	}

	// Apply the permutation: x starts as P·b.
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = bf[perm[i]]
	}
	// Forward substitution L·y = P·b (L is unit lower triangular).
	for i := 0; i < n; i++ {
		bi := i * n
		s := x[i]
		for j := 0; j < i; j++ {
			s -= a[bi+j] * x[j]
		}
		x[i] = s
	}
	// Back substitution U·x = y.
	for i := n - 1; i >= 0; i-- {
		bi := i * n
		s := x[i]
		for j := i + 1; j < n; j++ {
			s -= a[bi+j] * x[j]
		}
		x[i] = s / a[bi+i]
	}

	out := make([]algebra.Expr, n)
	for i := 0; i < n; i++ {
		out[i] = algebra.Flt(x[i])
	}
	return NewVector(out...), nil
}

// QR computes the Householder QR factorization A = Q·R of the (numeric) matrix
// m, which must have at least as many rows as columns (m ≥ n). Q is an m×m
// orthogonal matrix and R is an m×n upper-triangular matrix. It returns
// [ErrDimension] if m has fewer rows than columns and [ErrUnsupported]
// (wrapping the evaluation error) if any entry has a free symbol. The reflector
// signs are chosen deterministically to avoid cancellation.
func (m *Matrix) QR() (q, r *Matrix, err error) {
	rows, cols := m.Rows(), m.Cols()
	if rows < cols {
		return nil, nil, ErrDimension
	}
	rf, _, _, err := matrixFloatsFlat(m)
	if err != nil {
		return nil, nil, err
	}

	// Q starts as the m×m identity; R starts as a copy of A (rf, row-major).
	qf := make([]float64, rows*rows)
	for i := 0; i < rows; i++ {
		qf[i*rows+i] = 1
	}

	tol := matrixMaxAbs(rf) * float64(rows) * matrixEps
	v := make([]float64, rows) // reused reflector vector

	for k := 0; k < cols && k < rows; k++ {
		// Column norm of R[k:, k].
		norm := 0.0
		for i := k; i < rows; i++ {
			val := rf[i*cols+k]
			norm += val * val
		}
		norm = math.Sqrt(norm)
		if norm <= tol {
			continue
		}
		// Deterministic sign: reflect away from the leading entry.
		x0 := rf[k*cols+k]
		alpha := -norm
		if x0 < 0 {
			alpha = norm
		}
		// Build reflector v (nonzero only on rows k..rows-1).
		v[k] = x0 - alpha
		for i := k + 1; i < rows; i++ {
			v[i] = rf[i*cols+k]
		}
		vnorm2 := 0.0
		for i := k; i < rows; i++ {
			vnorm2 += v[i] * v[i]
		}
		if vnorm2 <= tol*tol {
			continue
		}
		// Apply H = I - 2·v·vᵀ/vnorm2 to R (columns k..cols-1).
		for j := k; j < cols; j++ {
			dot := 0.0
			for i := k; i < rows; i++ {
				dot += v[i] * rf[i*cols+j]
			}
			f := 2 * dot / vnorm2
			for i := k; i < rows; i++ {
				rf[i*cols+j] -= f * v[i]
			}
		}
		// Accumulate Q := Q·H (right-multiplication builds Q = H_0·H_1···).
		for row := 0; row < rows; row++ {
			dot := 0.0
			for i := k; i < rows; i++ {
				dot += qf[row*rows+i] * v[i]
			}
			f := 2 * dot / vnorm2
			for i := k; i < rows; i++ {
				qf[row*rows+i] -= f * v[i]
			}
		}
	}

	// Force the strict lower triangle of R to exact zeros.
	for i := 0; i < rows; i++ {
		for j := 0; j < cols && j < i; j++ {
			rf[i*cols+j] = 0
		}
	}
	return FromFloats(matrixReshape(qf, rows, rows)), FromFloats(matrixReshape(rf, rows, cols)), nil
}

// Cholesky returns the lower-triangular factor L of the (numeric) symmetric
// positive-definite matrix m, such that A = L·Lᵀ. It returns [ErrNotSquare] if
// m is not square, [ErrUnsupported] (wrapping the evaluation error) if any
// entry has a free symbol, and [ErrNotPositiveDefinite] if m is not symmetric
// or not positive definite (a diagonal pivot at or below the scale-relative
// tolerance). The strict upper triangle of the result is zero.
func (m *Matrix) Cholesky() (*Matrix, error) {
	if !m.IsSquare() {
		return nil, ErrNotSquare
	}
	a, n, _, err := matrixFloatsFlat(m)
	if err != nil {
		return nil, err
	}
	tol := matrixMaxAbs(a) * float64(n) * matrixEps

	// Reject asymmetric input.
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if math.Abs(a[i*n+j]-a[j*n+i]) > tol {
				return nil, ErrNotPositiveDefinite
			}
		}
	}

	l := make([]float64, n*n)
	for i := 0; i < n; i++ {
		bi := i * n
		for j := 0; j <= i; j++ {
			bj := j * n
			s := a[bi+j]
			for k := 0; k < j; k++ {
				s -= l[bi+k] * l[bj+k]
			}
			if i == j {
				if s <= tol {
					return nil, ErrNotPositiveDefinite
				}
				l[bi+j] = math.Sqrt(s)
			} else {
				l[bi+j] = s / l[bj+j]
			}
		}
	}
	return FromFloats(matrixReshape(l, n, n)), nil
}
