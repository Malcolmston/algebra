package probability

import "math"

// probabilityCopyMatrix returns a deep copy of a rectangular matrix.
func probabilityCopyMatrix(a [][]float64) [][]float64 {
	out := make([][]float64, len(a))
	for i := range a {
		out[i] = make([]float64, len(a[i]))
		copy(out[i], a[i])
	}
	return out
}

// probabilityIdentity returns the n-by-n identity matrix.
func probabilityIdentity(n int) [][]float64 {
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
		m[i][i] = 1
	}
	return m
}

// probabilityMatMul multiplies an m-by-k matrix a by a k-by-n matrix b and
// returns the m-by-n product. It panics if the inner dimensions disagree, which
// only occurs on a programming error inside this package.
func probabilityMatMul(a, b [][]float64) [][]float64 {
	m := len(a)
	k := len(b)
	n := 0
	if k > 0 {
		n = len(b[0])
	}
	out := make([][]float64, m)
	for i := 0; i < m; i++ {
		out[i] = make([]float64, n)
		for t := 0; t < k; t++ {
			aij := a[i][t]
			if aij == 0 {
				continue
			}
			for j := 0; j < n; j++ {
				out[i][j] += aij * b[t][j]
			}
		}
	}
	return out
}

// probabilityMatVec multiplies matrix a by column vector v.
func probabilityMatVec(a [][]float64, v []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		s := 0.0
		for j := range v {
			s += a[i][j] * v[j]
		}
		out[i] = s
	}
	return out
}

// probabilityVecMat multiplies row vector v by matrix a (v has length equal to
// the number of rows of a) and returns the resulting row vector.
func probabilityVecMat(v []float64, a [][]float64) []float64 {
	n := 0
	if len(a) > 0 {
		n = len(a[0])
	}
	out := make([]float64, n)
	for i := range v {
		vi := v[i]
		if vi == 0 {
			continue
		}
		for j := 0; j < n; j++ {
			out[j] += vi * a[i][j]
		}
	}
	return out
}

// probabilityMatPow returns a raised to the non-negative integer power p using
// exponentiation by squaring. probabilityMatPow(a, 0) is the identity.
func probabilityMatPow(a [][]float64, p int) [][]float64 {
	n := len(a)
	result := probabilityIdentity(n)
	base := probabilityCopyMatrix(a)
	for p > 0 {
		if p&1 == 1 {
			result = probabilityMatMul(result, base)
		}
		p >>= 1
		if p > 0 {
			base = probabilityMatMul(base, base)
		}
	}
	return result
}

// probabilitySolve solves the linear system A x = b for a square matrix A using
// Gaussian elimination with partial pivoting. A and b are not modified. It
// returns an error if A is singular.
func probabilitySolve(a [][]float64, b []float64) ([]float64, error) {
	n := len(a)
	m := make([][]float64, n)
	for i := range a {
		m[i] = make([]float64, n+1)
		copy(m[i], a[i])
		m[i][n] = b[i]
	}
	for col := 0; col < n; col++ {
		// Partial pivot: find the row with the largest magnitude in this column.
		piv := col
		best := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if best < 1e-300 {
			return nil, probabilityErrorf("probabilitySolve: singular matrix")
		}
		m[col], m[piv] = m[piv], m[col]
		// Eliminate below.
		for r := col + 1; r < n; r++ {
			f := m[r][col] / m[col][col]
			if f == 0 {
				continue
			}
			for c := col; c <= n; c++ {
				m[r][c] -= f * m[col][c]
			}
		}
	}
	// Back-substitution.
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		s := m[i][n]
		for j := i + 1; j < n; j++ {
			s -= m[i][j] * x[j]
		}
		x[i] = s / m[i][i]
	}
	return x, nil
}

// probabilityInverse returns the inverse of a square matrix using Gauss-Jordan
// elimination with partial pivoting. It returns an error if the matrix is
// singular.
func probabilityInverse(a [][]float64) ([][]float64, error) {
	n := len(a)
	// Augment [A | I].
	m := make([][]float64, n)
	for i := range a {
		m[i] = make([]float64, 2*n)
		copy(m[i], a[i])
		m[i][n+i] = 1
	}
	for col := 0; col < n; col++ {
		piv := col
		best := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if best < 1e-300 {
			return nil, probabilityErrorf("probabilityInverse: singular matrix")
		}
		m[col], m[piv] = m[piv], m[col]
		// Normalize pivot row.
		pv := m[col][col]
		for c := 0; c < 2*n; c++ {
			m[col][c] /= pv
		}
		// Eliminate all other rows.
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := m[r][col]
			if f == 0 {
				continue
			}
			for c := 0; c < 2*n; c++ {
				m[r][c] -= f * m[col][c]
			}
		}
	}
	inv := make([][]float64, n)
	for i := 0; i < n; i++ {
		inv[i] = make([]float64, n)
		copy(inv[i], m[i][n:])
	}
	return inv, nil
}
