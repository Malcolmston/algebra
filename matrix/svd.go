package matrix

import (
	"math"
	"sort"
)

// SVD returns the thin (economy) singular value decomposition A = U·diag(s)·Vᵀ
// of a numeric matrix. The singular values s are returned in non-increasing
// order; U is m×k and V is n×k with orthonormal columns, where k = min(m, n).
//
// The decomposition is computed with a deterministic one-sided Jacobi iteration
// using only the standard library. The sweep order and convergence threshold are
// fixed, so the result is reproducible run to run. The sign of each singular
// triplet is pinned by making the largest-magnitude entry of every column of U
// non-negative.
//
// SVD returns [ErrUnsupported] if any entry is symbolic (cannot be evaluated to
// a float64). For an exactly rank-deficient input, the columns of U that pair
// with a zero singular value are zero rather than an arbitrary orthonormal
// completion; the singular values, V, and every SVD-derived quantity are still
// correct.
func (m *Matrix) SVD() (u *Matrix, s []float64, v *Matrix, err error) {
	uf, sv, vf, mm, nn, k, err := m.matrixSVDFloats()
	if err != nil {
		return nil, nil, nil, err
	}
	u = matrixFlatToMatrix(uf, mm, k)
	v = matrixFlatToMatrix(vf, nn, k)
	return u, sv, v, nil
}

// SingularValues returns the singular values of a numeric matrix in
// non-increasing order. There are min(m, n) of them. It returns [ErrUnsupported]
// if any entry is symbolic.
func (m *Matrix) SingularValues() ([]float64, error) {
	_, s, _, _, _, _, err := m.matrixSVDFloats()
	if err != nil {
		return nil, err
	}
	return s, nil
}

// RankNumeric returns the numerical rank of a numeric matrix: the number of
// singular values strictly greater than tol. A non-positive tol selects the
// default threshold max(m, n)·eps·s[0], where s[0] is the largest singular value
// and eps is the machine epsilon. It returns [ErrUnsupported] if any entry is
// symbolic.
func (m *Matrix) RankNumeric(tol float64) (int, error) {
	_, s, _, _, _, k, err := m.matrixSVDFloats()
	if err != nil {
		return 0, err
	}
	return matrixNumericRank(s, k, tol, m.rows, m.cols), nil
}

// Nullspace returns an orthonormal basis of the kernel of a numeric matrix as
// the columns of the returned matrix. The basis vectors are the right singular
// vectors that pair with singular values at or below the rank threshold, i.e.
// the columns of V past the numeric rank. A non-positive tol selects the same
// default threshold as [Matrix.RankNumeric]. The result is an n×0 matrix when
// the matrix has full numeric rank. It returns [ErrUnsupported] if any entry is
// symbolic.
func (m *Matrix) Nullspace(tol float64) (*Matrix, error) {
	_, s, vf, _, nn, k, err := m.matrixSVDFloats()
	if err != nil {
		return nil, err
	}
	rank := matrixNumericRank(s, k, tol, m.rows, m.cols)
	ncols := k - rank
	vals := make([][]float64, nn)
	for i := 0; i < nn; i++ {
		vals[i] = make([]float64, ncols)
		for c := 0; c < ncols; c++ {
			vals[i][c] = vf[i*k+(rank+c)]
		}
	}
	return FromFloats(vals), nil
}

// Pinv returns the Moore-Penrose pseudoinverse V·diag(1/s)·Uᵀ of a numeric
// matrix, an n×m matrix. Only singular values strictly greater than tol
// contribute; smaller ones are treated as zero and dropped. A non-positive tol
// selects the same default threshold as [Matrix.RankNumeric]. It returns
// [ErrUnsupported] if any entry is symbolic.
func (m *Matrix) Pinv(tol float64) (*Matrix, error) {
	uf, s, vf, mm, nn, k, err := m.matrixSVDFloats()
	if err != nil {
		return nil, err
	}
	if tol <= 0 {
		tol = matrixDefaultTol(s, k, m.rows, m.cols)
	}
	res := make([][]float64, nn)
	for i := 0; i < nn; i++ {
		res[i] = make([]float64, mm)
	}
	for c := 0; c < k; c++ {
		if s[c] <= tol {
			continue
		}
		inv := 1.0 / s[c]
		for i := 0; i < nn; i++ {
			vic := vf[i*k+c] * inv
			for j := 0; j < mm; j++ {
				res[i][j] += vic * uf[j*k+c]
			}
		}
	}
	return FromFloats(res), nil
}

// Cond2 returns the spectral (2-norm) condition number of a numeric matrix, the
// ratio s[0]/s[k-1] of its largest to smallest singular value with k = min(m, n).
// It returns math.Inf(1) when the smallest singular value is numerically zero
// (relative to the largest), which includes the zero and empty matrices. It
// returns [ErrUnsupported] if any entry is symbolic.
func (m *Matrix) Cond2() (float64, error) {
	_, s, _, _, _, k, err := m.matrixSVDFloats()
	if err != nil {
		return 0, err
	}
	if k == 0 || s[0] == 0 {
		return math.Inf(1), nil
	}
	smin := s[k-1]
	if smin <= matrixEps*s[0] {
		return math.Inf(1), nil
	}
	return s[0] / smin, nil
}

// matrixSVDFloats computes the thin SVD in flat row-major float buffers. It
// returns U (mm×k), the singular values (length k), V (nn×k), the dimensions,
// and k = min(mm, nn). Both U and V use k as their row stride. It returns
// [ErrUnsupported] if any entry is symbolic.
func (m *Matrix) matrixSVDFloats() (uf, s, vf []float64, mm, nn, k int, err error) {
	a, ferr := m.Floats()
	if ferr != nil {
		return nil, nil, nil, 0, 0, 0, ErrUnsupported
	}
	mm, nn = m.rows, m.cols
	k = min(mm, nn)
	if k == 0 {
		return []float64{}, []float64{}, []float64{}, mm, nn, 0, nil
	}

	// One-sided Jacobi needs at least as many rows as columns. When the matrix
	// is wide (mm < nn) we factor Aᵀ instead and swap the roles of U and V.
	transposed := mm < nn
	var p, q int
	var af []float64
	if !transposed {
		p, q = mm, nn
		af = make([]float64, p*q)
		for i := 0; i < mm; i++ {
			for j := 0; j < nn; j++ {
				af[i*q+j] = a[i][j]
			}
		}
	} else {
		p, q = nn, mm
		af = make([]float64, p*q)
		for i := 0; i < mm; i++ {
			for j := 0; j < nn; j++ {
				af[j*q+i] = a[i][j]
			}
		}
	}

	uJ, sJ, vJ := matrixJacobiSVD(af, p, q)
	// After the swap, both uf and vf have row stride k (== q).
	if !transposed {
		uf, vf = uJ, vJ
	} else {
		uf, vf = vJ, uJ
	}
	s = sJ

	// Sign convention: force the largest-magnitude entry of each U column to be
	// non-negative, flipping the paired V column to preserve A = U·diag(s)·Vᵀ.
	for c := 0; c < k; c++ {
		maxAbs := 0.0
		maxRow := 0
		for r := 0; r < mm; r++ {
			if av := math.Abs(uf[r*k+c]); av > maxAbs {
				maxAbs = av
				maxRow = r
			}
		}
		if uf[maxRow*k+c] < 0 {
			for r := 0; r < mm; r++ {
				uf[r*k+c] = -uf[r*k+c]
			}
			for r := 0; r < nn; r++ {
				vf[r*k+c] = -vf[r*k+c]
			}
		}
	}
	return uf, s, vf, mm, nn, k, nil
}

// matrixJacobiSVD computes the thin SVD of the p×q (p ≥ q) row-major buffer a by
// deterministic one-sided Jacobi rotation. It orthogonalizes the columns of a in
// place, accumulating the right factor in v; on return the column norms are the
// singular values and the normalized columns are U. The rotations touch only the
// two columns of each (i, j) pair and reuse a single pair of scratch slices for
// the whole computation, so no memory is allocated inside the sweep loop.
//
// The returned buffers are freshly ordered so the singular values are
// non-increasing: uf is p×q (U), s has length q, and vf is q×q (V).
func matrixJacobiSVD(a []float64, p, q int) (uf, s, vf []float64) {
	v := make([]float64, q*q)
	for i := 0; i < q; i++ {
		v[i*q+i] = 1
	}
	// Scratch columns reused across every rotation of every sweep. Length p
	// covers both a (p rows) and v (q ≤ p rows).
	ci := make([]float64, p)
	cj := make([]float64, p)

	for sweep := 0; sweep < matrixJacobiMaxSweeps; sweep++ {
		rotated := false
		for i := 0; i < q; i++ {
			for j := i + 1; j < q; j++ {
				var alpha, beta, gamma float64
				for r := 0; r < p; r++ {
					x := a[r*q+i]
					y := a[r*q+j]
					ci[r] = x
					cj[r] = y
					alpha += x * x
					beta += y * y
					gamma += x * y
				}
				if gamma == 0 {
					continue
				}
				// Skip pairs already orthogonal to working precision.
				if math.Abs(gamma) <= matrixEps*math.Sqrt(alpha*beta) {
					continue
				}
				rotated = true

				// Jacobi rotation that diagonalizes the 2×2 Gram block.
				zeta := (beta - alpha) / (2 * gamma)
				var t float64
				if zeta >= 0 {
					t = 1.0 / (zeta + math.Sqrt(1+zeta*zeta))
				} else {
					t = -1.0 / (-zeta + math.Sqrt(1+zeta*zeta))
				}
				c := 1.0 / math.Sqrt(1+t*t)
				sn := c * t

				for r := 0; r < p; r++ {
					x := ci[r]
					y := cj[r]
					a[r*q+i] = c*x - sn*y
					a[r*q+j] = sn*x + c*y
				}
				for r := 0; r < q; r++ {
					x := v[r*q+i]
					y := v[r*q+j]
					v[r*q+i] = c*x - sn*y
					v[r*q+j] = sn*x + c*y
				}
			}
		}
		if !rotated {
			break
		}
	}

	// Column norms are the singular values; normalize to obtain U in place.
	s = make([]float64, q)
	for j := 0; j < q; j++ {
		var nrm float64
		for r := 0; r < p; r++ {
			x := a[r*q+j]
			nrm += x * x
		}
		nrm = math.Sqrt(nrm)
		s[j] = nrm
		if nrm > 0 {
			for r := 0; r < p; r++ {
				a[r*q+j] /= nrm
			}
		}
	}

	// Order the triplets by non-increasing singular value.
	idx := make([]int, q)
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(x, y int) bool { return s[idx[x]] > s[idx[y]] })

	uf = make([]float64, p*q)
	vf = make([]float64, q*q)
	ns := make([]float64, q)
	for newc, oldc := range idx {
		ns[newc] = s[oldc]
		for r := 0; r < p; r++ {
			uf[r*q+newc] = a[r*q+oldc]
		}
		for r := 0; r < q; r++ {
			vf[r*q+newc] = v[r*q+oldc]
		}
	}
	return uf, ns, vf
}

// matrixDefaultTol returns the default singular-value threshold
// max(rows, cols)·eps·s[0], or 0 when there are no singular values.
func matrixDefaultTol(s []float64, k, rows, cols int) float64 {
	if k == 0 {
		return 0
	}
	return float64(max(rows, cols)) * matrixEps * s[0]
}

// matrixNumericRank counts the singular values strictly above tol, applying the
// default threshold from matrixDefaultTol when tol is non-positive.
func matrixNumericRank(s []float64, k int, tol float64, rows, cols int) int {
	if k == 0 {
		return 0
	}
	if tol <= 0 {
		tol = matrixDefaultTol(s, k, rows, cols)
	}
	r := 0
	for i := 0; i < k; i++ {
		if s[i] > tol {
			r++
		}
	}
	return r
}

// matrixFlatToMatrix wraps a flat row-major buffer of rows×cols floats into a
// numeric *Matrix.
func matrixFlatToMatrix(buf []float64, rows, cols int) *Matrix {
	vals := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		vals[i] = make([]float64, cols)
		for j := 0; j < cols; j++ {
			vals[i][j] = buf[i*cols+j]
		}
	}
	return FromFloats(vals)
}
