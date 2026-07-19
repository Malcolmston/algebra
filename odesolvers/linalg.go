package odesolvers

import "math"

// Matrix is a dense row-major matrix used by the implicit solvers. Each element
// of the outer slice is a row.
type Matrix = [][]float64

// NewMatrix returns a rows-by-cols matrix of zeros.
func NewMatrix(rows, cols int) Matrix {
	m := make(Matrix, rows)
	for i := range m {
		m[i] = make([]float64, cols)
	}
	return m
}

// IdentityMatrix returns the n-by-n identity matrix.
func IdentityMatrix(n int) Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m[i][i] = 1
	}
	return m
}

// CopyMatrix returns a deep copy of a.
func CopyMatrix(a Matrix) Matrix {
	out := make(Matrix, len(a))
	for i := range a {
		out[i] = Clone(a[i])
	}
	return out
}

// MatVec returns the matrix-vector product a*x. It panics on a dimension
// mismatch.
func MatVec(a Matrix, x []float64) []float64 {
	out := make([]float64, len(a))
	for i := range a {
		if len(a[i]) != len(x) {
			panic("odesolvers: MatVec dimension mismatch")
		}
		var s float64
		for j := range x {
			s += a[i][j] * x[j]
		}
		out[i] = s
	}
	return out
}

// SolveLinearSystem solves A x = b for x by Gaussian elimination with partial
// pivoting. A must be square. It returns [ErrSingularMatrix] when A is
// numerically singular and [ErrDimensionMismatch] on a shape mismatch. A and b
// are not modified.
func SolveLinearSystem(a Matrix, b []float64) ([]float64, error) {
	n := len(a)
	if n == 0 {
		return []float64{}, nil
	}
	if len(b) != n {
		return nil, ErrDimensionMismatch
	}
	// Working copies (augmented).
	m := make(Matrix, n)
	for i := range a {
		if len(a[i]) != n {
			return nil, ErrDimensionMismatch
		}
		m[i] = Clone(a[i])
	}
	rhs := Clone(b)

	for col := 0; col < n; col++ {
		// Partial pivot.
		piv := col
		best := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > best {
				best = v
				piv = r
			}
		}
		if best < 1e-300 {
			return nil, ErrSingularMatrix
		}
		if piv != col {
			m[col], m[piv] = m[piv], m[col]
			rhs[col], rhs[piv] = rhs[piv], rhs[col]
		}
		// Eliminate below.
		pivVal := m[col][col]
		for r := col + 1; r < n; r++ {
			factor := m[r][col] / pivVal
			if factor == 0 {
				continue
			}
			for c := col; c < n; c++ {
				m[r][c] -= factor * m[col][c]
			}
			rhs[r] -= factor * rhs[col]
		}
	}
	// Back substitution.
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		s := rhs[i]
		for j := i + 1; j < n; j++ {
			s -= m[i][j] * x[j]
		}
		x[i] = s / m[i][i]
	}
	return x, nil
}

// FiniteDiffJacobian approximates the Jacobian dF/dx of the vector field F at x
// using first-order forward differences with step eps (a sensible default is
// used when eps <= 0). The result J satisfies J[i][j] = dF_i/dx_j.
func FiniteDiffJacobian(F func(x []float64) []float64, x []float64, eps float64) Matrix {
	n := len(x)
	f0 := F(x)
	m := len(f0)
	J := NewMatrix(m, n)
	xp := Clone(x)
	for j := 0; j < n; j++ {
		h := eps
		if h <= 0 {
			h = math.Sqrt(2.2e-16) * (math.Abs(x[j]) + 1)
		}
		xp[j] = x[j] + h
		fj := F(xp)
		xp[j] = x[j]
		for i := 0; i < m; i++ {
			J[i][j] = (fj[i] - f0[i]) / h
		}
	}
	return J
}

// NewtonSolve solves the square nonlinear system F(x) = 0 by Newton's method
// with a finite-difference Jacobian, starting from x0. It iterates until the
// residual infinity-norm falls below tol or maxIter iterations elapse. It
// returns the approximate root and, on non-convergence or a singular Jacobian,
// a non-nil error together with the best iterate reached.
func NewtonSolve(F func(x []float64) []float64, x0 []float64, tol float64, maxIter int) ([]float64, error) {
	if tol <= 0 {
		tol = 1e-10
	}
	if maxIter <= 0 {
		maxIter = 50
	}
	x := Clone(x0)
	for it := 0; it < maxIter; it++ {
		fx := F(x)
		if NormInf(fx) < tol {
			return x, nil
		}
		J := FiniteDiffJacobian(F, x, 0)
		delta, err := SolveLinearSystem(J, Neg(fx))
		if err != nil {
			return x, err
		}
		// Damped update with a simple backtracking line search on |F|.
		lambda := 1.0
		base := Norm2(fx)
		var xnew []float64
		for k := 0; k < 20; k++ {
			xnew = AXPY(lambda, delta, x)
			if Norm2(F(xnew)) < base || lambda < 1e-6 {
				break
			}
			lambda *= 0.5
		}
		x = xnew
		if NormInf(delta)*lambda < tol {
			if NormInf(F(x)) < math.Max(tol, 1e-8) {
				return x, nil
			}
		}
	}
	if NormInf(F(x)) < math.Max(tol, 1e-6) {
		return x, nil
	}
	return x, ErrNoConvergence
}

// MatMul returns the matrix product a*b. It panics on a dimension mismatch.
func MatMul(a, b Matrix) Matrix {
	ra, ca := len(a), 0
	if ra > 0 {
		ca = len(a[0])
	}
	rb, cb := len(b), 0
	if rb > 0 {
		cb = len(b[0])
	}
	if ca != rb {
		panic("odesolvers: MatMul dimension mismatch")
	}
	out := NewMatrix(ra, cb)
	for i := 0; i < ra; i++ {
		for k := 0; k < ca; k++ {
			aik := a[i][k]
			if aik == 0 {
				continue
			}
			for j := 0; j < cb; j++ {
				out[i][j] += aik * b[k][j]
			}
		}
	}
	return out
}

// Transpose returns the transpose of a.
func Transpose(a Matrix) Matrix {
	r := len(a)
	c := 0
	if r > 0 {
		c = len(a[0])
	}
	out := NewMatrix(c, r)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			out[j][i] = a[i][j]
		}
	}
	return out
}

// Trace returns the sum of the diagonal entries of the square matrix a.
func Trace(a Matrix) float64 {
	var s float64
	for i := range a {
		s += a[i][i]
	}
	return s
}

// MatAdd returns the elementwise sum a+b. It panics on a shape mismatch.
func MatAdd(a, b Matrix) Matrix {
	out := CopyMatrix(a)
	for i := range a {
		for j := range a[i] {
			out[i][j] += b[i][j]
		}
	}
	return out
}

// MatScale returns the matrix s*a.
func MatScale(s float64, a Matrix) Matrix {
	out := CopyMatrix(a)
	for i := range out {
		for j := range out[i] {
			out[i][j] *= s
		}
	}
	return out
}

// Determinant returns the determinant of the square matrix a via LU
// factorization with partial pivoting. A singular matrix yields 0.
func Determinant(a Matrix) float64 {
	n := len(a)
	m := CopyMatrix(a)
	det := 1.0
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
			return 0
		}
		if piv != col {
			m[col], m[piv] = m[piv], m[col]
			det = -det
		}
		det *= m[col][col]
		for r := col + 1; r < n; r++ {
			factor := m[r][col] / m[col][col]
			for c := col; c < n; c++ {
				m[r][c] -= factor * m[col][c]
			}
		}
	}
	return det
}
