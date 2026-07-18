package matrix

import (
	"fmt"
	"math"

	"github.com/malcolmston/algebra"
)

// LeastSquares solves the m×n linear system A·x = b in the least-squares sense
// and returns the solution vector x.
//
// The routine dispatches on the shape of A:
//
//   - when A has at least as many rows as columns (m ≥ n) and full column rank,
//     it returns the x that minimizes the Euclidean residual ‖A·x − b‖₂,
//     computed by an in-place Householder QR factorization followed by
//     back-substitution against the triangular factor;
//   - when A has fewer rows than columns (m < n) and full row rank, it returns
//     the minimum-Euclidean-norm solution among the infinitely many exact
//     solutions of A·x = b, computed from the Householder QR factorization of
//     Aᵀ.
//
// The computation is numeric: every entry of A and b is evaluated to a float64.
// It returns [ErrDimension] when len(b) does not equal A.Rows(), [ErrUnsupported]
// (wrapping the underlying evaluation error) when any entry of A or b contains a
// free symbol, and [ErrSingular] when A is rank deficient so that no unique
// least-squares or minimum-norm solution is defined. The returned components are
// inexact algebra.Flt literals.
//
// For performance the factorization runs entirely inside flat row-major
// []float64 buffers with in-place Householder reflections and reused
// reflector/scratch buffers; no intermediate [Matrix] values are allocated
// between the factorization and the triangular solve.
func LeastSquares(a *Matrix, b *Vector) (*Vector, error) {
	m, n := a.rows, a.cols
	if len(b.data) != m {
		return nil, ErrDimension
	}
	af, _, _, err := matrixFloatsFlat(a)
	if err != nil {
		return nil, err
	}
	bf := make([]float64, m)
	for i := 0; i < m; i++ {
		f, e := algebra.Evalf(b.data[i])
		if e != nil {
			return nil, fmt.Errorf("%w: %v", ErrUnsupported, e)
		}
		bf[i] = f
	}
	if m >= n {
		return matrixLeastSquaresOver(af, bf, m, n)
	}
	return matrixLeastSquaresUnder(af, bf, m, n)
}

// matrixLeastSquaresOver returns the minimizer of ‖A·x − b‖₂ for an
// overdetermined (or square) system, where af holds the m×n coefficient matrix
// in row-major order and bf holds the length-m right-hand side. It reflects A
// to upper-triangular form in place with a single reused reflector buffer,
// applying each Householder reflection to bf as it goes so that bf becomes Qᵀ·b,
// then back-substitutes against the triangular factor. It returns [ErrSingular]
// if the matrix is rank deficient.
func matrixLeastSquaresOver(af, bf []float64, m, n int) (*Vector, error) {
	tol := matrixMaxAbs(af) * float64(m) * matrixEps
	v := make([]float64, m) // reused reflector / scratch buffer
	for k := 0; k < n; k++ {
		// Euclidean norm of the sub-column af[k:, k].
		norm := 0.0
		for i := k; i < m; i++ {
			val := af[i*n+k]
			norm += val * val
		}
		norm = math.Sqrt(norm)
		if norm <= tol {
			return nil, ErrSingular
		}
		// Deterministic sign: reflect away from the leading entry.
		x0 := af[k*n+k]
		alpha := -norm
		if x0 < 0 {
			alpha = norm
		}
		v[k] = x0 - alpha
		for i := k + 1; i < m; i++ {
			v[i] = af[i*n+k]
		}
		vnorm2 := 0.0
		for i := k; i < m; i++ {
			vnorm2 += v[i] * v[i]
		}
		if vnorm2 <= tol*tol {
			return nil, ErrSingular
		}
		// Apply H = I − 2·v·vᵀ/vnorm2 to the trailing columns of A.
		for j := k; j < n; j++ {
			dot := 0.0
			for i := k; i < m; i++ {
				dot += v[i] * af[i*n+j]
			}
			f := 2 * dot / vnorm2
			for i := k; i < m; i++ {
				af[i*n+j] -= f * v[i]
			}
		}
		// Apply the same reflection to b so bf accumulates Qᵀ·b.
		dot := 0.0
		for i := k; i < m; i++ {
			dot += v[i] * bf[i]
		}
		f := 2 * dot / vnorm2
		for i := k; i < m; i++ {
			bf[i] -= f * v[i]
		}
	}
	// Back-substitution R·x = (Qᵀ·b)[0:n]; R is the top n×n upper triangle of af.
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		s := bf[i]
		for j := i + 1; j < n; j++ {
			s -= af[i*n+j] * x[j]
		}
		x[i] = s / af[i*n+i]
	}
	out := make([]algebra.Expr, n)
	for i := 0; i < n; i++ {
		out[i] = algebra.Flt(x[i])
	}
	return &Vector{data: out}, nil
}

// matrixLeastSquaresUnder returns the minimum-Euclidean-norm solution of the
// underdetermined system A·x = b, where af holds the m×n coefficient matrix in
// row-major order (m < n) and bf holds the length-m right-hand side. It forms Aᵀ
// in a flat buffer, factors it with in-place Householder reflections whose
// reflector vectors are stored in a reused flat buffer, solves R₁ᵀ·z = b by
// forward-substitution, and recovers x = Q·[z; 0] by replaying the stored
// reflectors in reverse. It returns [ErrSingular] if A is not of full row rank.
func matrixLeastSquaresUnder(af, bf []float64, m, n int) (*Vector, error) {
	// Aᵀ is n×m in row-major order.
	at := make([]float64, n*m)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			at[j*m+i] = af[i*n+j]
		}
	}
	tol := matrixMaxAbs(at) * float64(n) * matrixEps
	refl := make([]float64, m*n)  // reflector k occupies refl[k*n : k*n+n]
	vnorm2s := make([]float64, m) // squared norms of the reflectors
	for k := 0; k < m; k++ {
		norm := 0.0
		for i := k; i < n; i++ {
			val := at[i*m+k]
			norm += val * val
		}
		norm = math.Sqrt(norm)
		if norm <= tol {
			return nil, ErrSingular
		}
		x0 := at[k*m+k]
		alpha := -norm
		if x0 < 0 {
			alpha = norm
		}
		base := k * n
		refl[base+k] = x0 - alpha
		for i := k + 1; i < n; i++ {
			refl[base+i] = at[i*m+k]
		}
		vnorm2 := 0.0
		for i := k; i < n; i++ {
			vnorm2 += refl[base+i] * refl[base+i]
		}
		if vnorm2 <= tol*tol {
			return nil, ErrSingular
		}
		vnorm2s[k] = vnorm2
		// Apply the reflection to the trailing columns of Aᵀ.
		for j := k; j < m; j++ {
			dot := 0.0
			for i := k; i < n; i++ {
				dot += refl[base+i] * at[i*m+j]
			}
			f := 2 * dot / vnorm2
			for i := k; i < n; i++ {
				at[i*m+j] -= f * refl[base+i]
			}
		}
	}
	// Solve R₁ᵀ·z = b by forward-substitution. R₁ is the top m×m upper triangle
	// of the factored Aᵀ, so R₁ᵀ[i][j] = at[j*m+i] for j ≤ i.
	z := make([]float64, m)
	for i := 0; i < m; i++ {
		s := bf[i]
		for j := 0; j < i; j++ {
			s -= at[j*m+i] * z[j]
		}
		z[i] = s / at[i*m+i]
	}
	// y = [z; 0] of length n, then x = Q·y = H₀·H₁···H_{m-1}·y. Replay the
	// stored reflectors from last to first, in place.
	y := make([]float64, n)
	copy(y[:m], z)
	for k := m - 1; k >= 0; k-- {
		base := k * n
		vnorm2 := vnorm2s[k]
		dot := 0.0
		for i := k; i < n; i++ {
			dot += refl[base+i] * y[i]
		}
		f := 2 * dot / vnorm2
		for i := k; i < n; i++ {
			y[i] -= f * refl[base+i]
		}
	}
	out := make([]algebra.Expr, n)
	for i := 0; i < n; i++ {
		out[i] = algebra.Flt(y[i])
	}
	return &Vector{data: out}, nil
}

// matrixPivotColumns returns, for the reduced row-echelon form rref of a matrix
// with n columns, two parallel slices giving the row and column index of every
// pivot (the first nonzero entry of each nonzero row) together with a length-n
// boolean marking the pivot columns. The pivots are reported in increasing row
// order, which is also increasing column order because rref is in echelon form.
func matrixPivotColumns(rref *Matrix, n int) (pivotRows, pivotCols []int, isPivot []bool) {
	isPivot = make([]bool, n)
	for r := 0; r < rref.rows; r++ {
		for j := 0; j < n; j++ {
			if !isZeroExpr(rref.data[r][j]) {
				pivotRows = append(pivotRows, r)
				pivotCols = append(pivotCols, j)
				isPivot[j] = true
				break
			}
		}
	}
	return pivotRows, pivotCols, isPivot
}

// NullspaceExact returns an exact basis of the kernel ker(A) = {x : A·x = 0}.
//
// The matrix is reduced to reduced row-echelon form with the exact symbolic
// [Matrix.RREF], so the result is valid even when entries contain free symbols,
// complementing the numeric SVD-based [Matrix.Nullspace]. One basis vector is
// produced for each free (non-pivot) column, taken in increasing column index:
// the free variable is set to 1, the other free variables to 0, and each pivot
// variable to the negation of the corresponding reduced entry. Every returned
// [Vector] has length A.Cols() and simplified entries. The basis is empty when A
// has full column rank. The error is always nil; it is present for signature
// symmetry with the numeric routines.
func (m *Matrix) NullspaceExact() ([]*Vector, error) {
	n := m.cols
	rref, _ := m.RREF()
	pivotRows, pivotCols, isPivot := matrixPivotColumns(rref, n)
	basis := make([]*Vector, 0, n)
	for f := 0; f < n; f++ {
		if isPivot[f] {
			continue
		}
		comps := make([]algebra.Expr, n)
		for i := range comps {
			comps[i] = zero()
		}
		comps[f] = one()
		for t := range pivotCols {
			comps[pivotCols[t]] = simp(algebra.Mul(algebra.Int(-1), rref.data[pivotRows[t]][f]))
		}
		basis = append(basis, &Vector{data: comps})
	}
	return basis, nil
}

// ColumnSpaceExact returns an exact basis of the column space (image) of A.
//
// The pivot columns are identified from the exact reduced row-echelon form and
// the corresponding original columns of A are returned unchanged, which form a
// basis of the column space. Columns are taken in increasing index and every
// returned [Vector] has length A.Rows() with simplified entries. This works with
// symbolic entries. The error is always nil; it is present for signature
// symmetry with the numeric routines.
func (m *Matrix) ColumnSpaceExact() ([]*Vector, error) {
	rref, _ := m.RREF()
	_, pivotCols, _ := matrixPivotColumns(rref, m.cols)
	basis := make([]*Vector, 0, len(pivotCols))
	for _, j := range pivotCols {
		col := make([]algebra.Expr, m.rows)
		for i := 0; i < m.rows; i++ {
			col[i] = simp(m.data[i][j])
		}
		basis = append(basis, &Vector{data: col})
	}
	return basis, nil
}

// RowSpaceExact returns an exact basis of the row space of A.
//
// The matrix is reduced to reduced row-echelon form with the exact symbolic
// [Matrix.RREF] and each nonzero row of the result is returned, in top-to-bottom
// order. These rows are linearly independent and span the same row space as A.
// Every returned [Vector] has length A.Cols() with simplified entries and works
// with symbolic entries. The error is always nil; it is present for signature
// symmetry with the numeric routines.
func (m *Matrix) RowSpaceExact() ([]*Vector, error) {
	n := m.cols
	rref, _ := m.RREF()
	basis := make([]*Vector, 0, rref.rows)
	for r := 0; r < rref.rows; r++ {
		allZero := true
		for j := 0; j < n; j++ {
			if !isZeroExpr(rref.data[r][j]) {
				allZero = false
				break
			}
		}
		if allZero {
			continue
		}
		row := make([]algebra.Expr, n)
		for j := 0; j < n; j++ {
			row[j] = simp(rref.data[r][j])
		}
		basis = append(basis, &Vector{data: row})
	}
	return basis, nil
}
