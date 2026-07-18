package matrix

import (
	"fmt"
	"math"

	"github.com/malcolmston/algebra"
)

// NormInfinity is the p selector passed to [Matrix.CondP] to request the
// infinity norm (maximum absolute row sum) instead of the 1-norm. It is a
// sentinel value distinct from any valid Lᵖ order.
const NormInfinity = -1

// NormFro returns the exact Frobenius norm of m, the square root of the sum of
// the squares of every entry. It is computed symbolically with [algebra.Sqrt],
// so an integer or rational matrix stays in simplest radical form (for example
// a matrix of ones of size 2×2 yields the literal 2, and one whose squared
// entries sum to 8 yields 2*sqrt(2)). The result is exact and works on symbolic
// entries.
func (m *Matrix) NormFro() algebra.Expr {
	terms := make([]algebra.Expr, 0, m.rows*m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			e := m.data[i][j]
			terms = append(terms, algebra.Mul(e, e))
		}
	}
	sum := simp(algebra.Add(terms...))
	return algebra.Sqrt(sum)
}

// Norm1 returns the induced 1-norm of m, the maximum absolute column sum. Every
// entry is evaluated to a float64 via [Matrix.Floats]; it returns
// [ErrUnsupported] if any entry contains a free symbol or is otherwise not
// numeric. An empty matrix has norm 0.
func (m *Matrix) Norm1() (float64, error) {
	n1, _, _, err := m.matrixNorms()
	if err != nil {
		return 0, err
	}
	return n1, nil
}

// NormInf returns the induced infinity-norm of m, the maximum absolute row sum.
// Every entry is evaluated to a float64 via [Matrix.Floats]; it returns
// [ErrUnsupported] if any entry contains a free symbol or is otherwise not
// numeric. An empty matrix has norm 0.
func (m *Matrix) NormInf() (float64, error) {
	_, nInf, _, err := m.matrixNorms()
	if err != nil {
		return 0, err
	}
	return nInf, nil
}

// NormMax returns the max-norm of m, the largest absolute value among all
// entries. Every entry is evaluated to a float64 via [Matrix.Floats]; it
// returns [ErrUnsupported] if any entry contains a free symbol or is otherwise
// not numeric. An empty matrix has norm 0.
func (m *Matrix) NormMax() (float64, error) {
	_, _, nMax, err := m.matrixNorms()
	if err != nil {
		return 0, err
	}
	return nMax, nil
}

// matrixNorms computes the induced 1-norm (maximum absolute column sum), the
// induced infinity-norm (maximum absolute row sum) and the max-norm (largest
// absolute entry) of m together in a single row-major pass over the float64
// values returned by [Matrix.Floats]. Column sums are kept in a single
// preallocated accumulator slice and the row sum and running maximum are scalar
// accumulators, so the inner loop performs no per-element allocation. It returns
// [ErrUnsupported] wrapping the underlying evaluation error when any entry is
// not numeric.
func (m *Matrix) matrixNorms() (norm1, normInf, normMax float64, err error) {
	f, ferr := m.Floats()
	if ferr != nil {
		return 0, 0, 0, fmt.Errorf("%w: %v", ErrUnsupported, ferr)
	}
	colSums := make([]float64, m.cols)
	for i := 0; i < m.rows; i++ {
		var rowSum float64
		row := f[i]
		for j := 0; j < m.cols; j++ {
			a := math.Abs(row[j])
			rowSum += a
			colSums[j] += a
			if a > normMax {
				normMax = a
			}
		}
		if rowSum > normInf {
			normInf = rowSum
		}
	}
	for _, c := range colSums {
		if c > norm1 {
			norm1 = c
		}
	}
	return norm1, normInf, normMax, nil
}

// Hadamard returns the entrywise (Schur) product m∘n, whose (i,j) entry is the
// product of the (i,j) entries of m and n, each simplified. It returns
// [ErrDimension] if the two matrices do not have identical shapes. The product
// is exact and works on symbolic entries.
func (m *Matrix) Hadamard(n *Matrix) (*Matrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrDimension
	}
	out := New(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[i][j] = simp(algebra.Mul(m.data[i][j], n.data[i][j]))
		}
	}
	return out, nil
}

// CondP returns the induced condition number of m for the norm selected by p,
// defined as ‖A‖ₚ · ‖A⁻¹‖ₚ. The only supported selectors are p == 1 for the
// 1-norm and p == [NormInfinity] for the infinity-norm; any other value yields
// [ErrUnsupported]. The inverse is obtained from [Matrix.Inverse] and both norms
// are evaluated numerically, so a non-numeric entry yields [ErrUnsupported]. It
// returns [ErrNotSquare] for a non-square matrix, and for a singular matrix it
// returns math.Inf(1) together with [ErrSingular].
func (m *Matrix) CondP(p int) (float64, error) {
	if !m.IsSquare() {
		return 0, ErrNotSquare
	}
	if p != 1 && p != NormInfinity {
		return 0, ErrUnsupported
	}
	inv, err := m.Inverse()
	if err != nil {
		if err == ErrSingular {
			return math.Inf(1), ErrSingular
		}
		return 0, err
	}
	var na, ni float64
	if p == 1 {
		if na, err = m.Norm1(); err != nil {
			return 0, err
		}
		if ni, err = inv.Norm1(); err != nil {
			return 0, err
		}
	} else {
		if na, err = m.NormInf(); err != nil {
			return 0, err
		}
		if ni, err = inv.NormInf(); err != nil {
			return 0, err
		}
	}
	return na * ni, nil
}
