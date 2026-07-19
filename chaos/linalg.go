package chaos

import (
	"math"
	"math/cmplx"
)

// Vec is a real vector represented as a slice of float64 coordinates.
type Vec []float64

// Mat is a dense real matrix stored as a slice of rows; every row must have
// the same length, which is the number of columns.
type Mat []Vec

// NewVec returns a zero vector of length n.
func NewVec(n int) Vec {
	return make(Vec, n)
}

// VecOf returns a Vec containing the given components.
func VecOf(xs ...float64) Vec {
	v := make(Vec, len(xs))
	copy(v, xs)
	return v
}

// Clone returns an independent copy of v.
func (v Vec) Clone() Vec {
	c := make(Vec, len(v))
	copy(c, v)
	return c
}

// Dim returns the number of components of v.
func (v Vec) Dim() int { return len(v) }

// Add returns the elementwise sum v+w. The shorter length governs the result.
func (v Vec) Add(w Vec) Vec {
	n := min(len(v), len(w))
	r := make(Vec, n)
	for i := 0; i < n; i++ {
		r[i] = v[i] + w[i]
	}
	return r
}

// Sub returns the elementwise difference v-w.
func (v Vec) Sub(w Vec) Vec {
	n := min(len(v), len(w))
	r := make(Vec, n)
	for i := 0; i < n; i++ {
		r[i] = v[i] - w[i]
	}
	return r
}

// Scale returns the vector s*v.
func (v Vec) Scale(s float64) Vec {
	r := make(Vec, len(v))
	for i, x := range v {
		r[i] = s * x
	}
	return r
}

// AddScaled returns v + s*w, the fused multiply-add of two vectors.
func (v Vec) AddScaled(s float64, w Vec) Vec {
	n := min(len(v), len(w))
	r := make(Vec, n)
	for i := 0; i < n; i++ {
		r[i] = v[i] + s*w[i]
	}
	return r
}

// Dot returns the Euclidean inner product of v and w.
func (v Vec) Dot(w Vec) float64 {
	n := min(len(v), len(w))
	var s float64
	for i := 0; i < n; i++ {
		s += v[i] * w[i]
	}
	return s
}

// Norm returns the Euclidean (L2) norm of v.
func (v Vec) Norm() float64 {
	return math.Sqrt(v.Dot(v))
}

// Norm1 returns the L1 (taxicab) norm of v.
func (v Vec) Norm1() float64 {
	var s float64
	for _, x := range v {
		s += math.Abs(x)
	}
	return s
}

// NormInf returns the maximum-absolute-value (Chebyshev) norm of v.
func (v Vec) NormInf() float64 {
	var m float64
	for _, x := range v {
		if a := math.Abs(x); a > m {
			m = a
		}
	}
	return m
}

// Distance returns the Euclidean distance between v and w.
func (v Vec) Distance(w Vec) float64 {
	return v.Sub(w).Norm()
}

// Normalize returns v scaled to unit Euclidean norm; a zero vector is returned
// unchanged.
func (v Vec) Normalize() Vec {
	n := v.Norm()
	if n == 0 {
		return v.Clone()
	}
	return v.Scale(1 / n)
}

// Max returns the largest component of v, or NaN if v is empty.
func (v Vec) Max() float64 {
	if len(v) == 0 {
		return math.NaN()
	}
	m := v[0]
	for _, x := range v[1:] {
		if x > m {
			m = x
		}
	}
	return m
}

// Min returns the smallest component of v, or NaN if v is empty.
func (v Vec) Min() float64 {
	if len(v) == 0 {
		return math.NaN()
	}
	m := v[0]
	for _, x := range v[1:] {
		if x < m {
			m = x
		}
	}
	return m
}

// Sum returns the sum of the components of v.
func (v Vec) Sum() float64 {
	var s float64
	for _, x := range v {
		s += x
	}
	return s
}

// Mean returns the arithmetic mean of the components of v, or NaN if empty.
func (v Vec) Mean() float64 {
	if len(v) == 0 {
		return math.NaN()
	}
	return v.Sum() / float64(len(v))
}

// NewMat returns an r-by-c matrix of zeros.
func NewMat(r, c int) Mat {
	m := make(Mat, r)
	for i := range m {
		m[i] = make(Vec, c)
	}
	return m
}

// Eye returns the n-by-n identity matrix.
func Eye(n int) Mat {
	m := NewMat(n, n)
	for i := 0; i < n; i++ {
		m[i][i] = 1
	}
	return m
}

// MatFromRows builds a matrix from the given rows, copying each one.
func MatFromRows(rows ...Vec) Mat {
	m := make(Mat, len(rows))
	for i, r := range rows {
		m[i] = r.Clone()
	}
	return m
}

// Rows returns the number of rows of A.
func (A Mat) Rows() int { return len(A) }

// Cols returns the number of columns of A (zero for an empty matrix).
func (A Mat) Cols() int {
	if len(A) == 0 {
		return 0
	}
	return len(A[0])
}

// Clone returns an independent copy of A.
func (A Mat) Clone() Mat {
	c := make(Mat, len(A))
	for i, r := range A {
		c[i] = r.Clone()
	}
	return c
}

// IsSquare reports whether A has equal numbers of rows and columns.
func (A Mat) IsSquare() bool { return A.Rows() == A.Cols() }

// Transpose returns the transpose of A.
func (A Mat) Transpose() Mat {
	r, c := A.Rows(), A.Cols()
	t := NewMat(c, r)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			t[j][i] = A[i][j]
		}
	}
	return t
}

// MulVec returns the matrix-vector product A*x.
func (A Mat) MulVec(x Vec) Vec {
	r := make(Vec, A.Rows())
	for i, row := range A {
		r[i] = row.Dot(x)
	}
	return r
}

// Mul returns the matrix product A*B.
func (A Mat) Mul(B Mat) Mat {
	n, m, p := A.Rows(), A.Cols(), B.Cols()
	C := NewMat(n, p)
	for i := 0; i < n; i++ {
		for k := 0; k < m; k++ {
			a := A[i][k]
			if a == 0 {
				continue
			}
			for j := 0; j < p; j++ {
				C[i][j] += a * B[k][j]
			}
		}
	}
	return C
}

// Add returns the elementwise sum A+B.
func (A Mat) Add(B Mat) Mat {
	r, c := A.Rows(), A.Cols()
	C := NewMat(r, c)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			C[i][j] = A[i][j] + B[i][j]
		}
	}
	return C
}

// Sub returns the elementwise difference A-B.
func (A Mat) Sub(B Mat) Mat {
	r, c := A.Rows(), A.Cols()
	C := NewMat(r, c)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			C[i][j] = A[i][j] - B[i][j]
		}
	}
	return C
}

// Scale returns the matrix s*A.
func (A Mat) Scale(s float64) Mat {
	C := NewMat(A.Rows(), A.Cols())
	for i := range A {
		for j := range A[i] {
			C[i][j] = s * A[i][j]
		}
	}
	return C
}

// Trace returns the sum of the diagonal entries of a square matrix.
func (A Mat) Trace() float64 {
	var s float64
	n := min(A.Rows(), A.Cols())
	for i := 0; i < n; i++ {
		s += A[i][i]
	}
	return s
}

// Col returns column j of A as a new vector.
func (A Mat) Col(j int) Vec {
	c := make(Vec, A.Rows())
	for i := range A {
		c[i] = A[i][j]
	}
	return c
}

// FrobeniusNorm returns the Frobenius norm of A.
func (A Mat) FrobeniusNorm() float64 {
	var s float64
	for _, row := range A {
		for _, x := range row {
			s += x * x
		}
	}
	return math.Sqrt(s)
}

// Det2 returns the determinant of the 2-by-2 matrix A.
func Det2(A Mat) float64 {
	return A[0][0]*A[1][1] - A[0][1]*A[1][0]
}

// Det3 returns the determinant of the 3-by-3 matrix A.
func Det3(A Mat) float64 {
	return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) -
		A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) +
		A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
}

// Det returns the determinant of a square matrix by Gaussian elimination with
// partial pivoting.
func Det(A Mat) float64 {
	n := A.Rows()
	if n == 0 || !A.IsSquare() {
		return math.NaN()
	}
	m := A.Clone()
	det := 1.0
	for k := 0; k < n; k++ {
		p := k
		max := math.Abs(m[k][k])
		for i := k + 1; i < n; i++ {
			if a := math.Abs(m[i][k]); a > max {
				max, p = a, i
			}
		}
		if max == 0 {
			return 0
		}
		if p != k {
			m[k], m[p] = m[p], m[k]
			det = -det
		}
		det *= m[k][k]
		for i := k + 1; i < n; i++ {
			f := m[i][k] / m[k][k]
			for j := k; j < n; j++ {
				m[i][j] -= f * m[k][j]
			}
		}
	}
	return det
}

// SolveLinear solves A*x = b for a square matrix A by Gaussian elimination
// with partial pivoting. It returns ErrSingular if A is singular.
func SolveLinear(A Mat, b Vec) (Vec, error) {
	n := A.Rows()
	if !A.IsSquare() {
		return nil, ErrNonSquare
	}
	if len(b) != n {
		return nil, ErrDimensionMismatch
	}
	m := A.Clone()
	x := b.Clone()
	for k := 0; k < n; k++ {
		p := k
		max := math.Abs(m[k][k])
		for i := k + 1; i < n; i++ {
			if a := math.Abs(m[i][k]); a > max {
				max, p = a, i
			}
		}
		if max < 1e-300 {
			return nil, ErrSingular
		}
		if p != k {
			m[k], m[p] = m[p], m[k]
			x[k], x[p] = x[p], x[k]
		}
		for i := k + 1; i < n; i++ {
			f := m[i][k] / m[k][k]
			for j := k; j < n; j++ {
				m[i][j] -= f * m[k][j]
			}
			x[i] -= f * x[k]
		}
	}
	for i := n - 1; i >= 0; i-- {
		s := x[i]
		for j := i + 1; j < n; j++ {
			s -= m[i][j] * x[j]
		}
		x[i] = s / m[i][i]
	}
	return x, nil
}

// Inverse returns the inverse of a square matrix, or ErrSingular if A is not
// invertible to working precision.
func Inverse(A Mat) (Mat, error) {
	n := A.Rows()
	if !A.IsSquare() {
		return nil, ErrNonSquare
	}
	inv := NewMat(n, n)
	for j := 0; j < n; j++ {
		e := make(Vec, n)
		e[j] = 1
		col, err := SolveLinear(A, e)
		if err != nil {
			return nil, err
		}
		for i := 0; i < n; i++ {
			inv[i][j] = col[i]
		}
	}
	return inv, nil
}

// QR computes the reduced QR factorisation of A by modified Gram-Schmidt,
// returning an orthonormal Q with the same shape as A and an upper-triangular
// R such that A = Q*R.
func QR(A Mat) (Q, R Mat) {
	n, m := A.Rows(), A.Cols()
	Q = NewMat(n, m)
	R = NewMat(m, m)
	// Work with columns.
	cols := make([]Vec, m)
	for j := 0; j < m; j++ {
		cols[j] = A.Col(j)
	}
	for j := 0; j < m; j++ {
		v := cols[j].Clone()
		for i := 0; i < j; i++ {
			qi := Q.Col(i)
			r := qi.Dot(v)
			R[i][j] = r
			v = v.AddScaled(-r, qi)
		}
		nrm := v.Norm()
		R[j][j] = nrm
		if nrm > 0 {
			v = v.Scale(1 / nrm)
		}
		for i := 0; i < n; i++ {
			Q[i][j] = v[i]
		}
	}
	return Q, R
}

// GramSchmidt orthonormalises the given set of vectors in place-independent
// fashion, returning a new orthonormal basis for their span (with zero vectors
// dropped to length zero when linearly dependent).
func GramSchmidt(vs []Vec) []Vec {
	out := make([]Vec, len(vs))
	for j := range vs {
		v := vs[j].Clone()
		for i := 0; i < j; i++ {
			if len(out[i]) == 0 {
				continue
			}
			v = v.AddScaled(-out[i].Dot(v), out[i])
		}
		n := v.Norm()
		if n > 1e-300 {
			out[j] = v.Scale(1 / n)
		} else {
			out[j] = make(Vec, len(v))
		}
	}
	return out
}

// Eigenvalues2 returns the two (possibly complex) eigenvalues of a 2-by-2
// matrix A.
func Eigenvalues2(A Mat) [2]complex128 {
	tr := A[0][0] + A[1][1]
	det := Det2(A)
	disc := complex(tr*tr-4*det, 0)
	s := cmplx.Sqrt(disc)
	half := complex(0.5, 0)
	l1 := (complex(tr, 0) + s) * half
	l2 := (complex(tr, 0) - s) * half
	return [2]complex128{l1, l2}
}

// Eigenvalues3 returns the three (possibly complex) eigenvalues of a 3-by-3
// matrix A, obtained by solving its characteristic cubic.
func Eigenvalues3(A Mat) [3]complex128 {
	// Characteristic polynomial: lambda^3 - c2 lambda^2 + c1 lambda - c0.
	c2 := A.Trace()
	c1 := A[0][0]*A[1][1] - A[0][1]*A[1][0] +
		A[0][0]*A[2][2] - A[0][2]*A[2][0] +
		A[1][1]*A[2][2] - A[1][2]*A[2][1]
	c0 := Det3(A)
	roots := solveCubic(1, -c2, c1, -c0)
	return roots
}

// SpectralRadius returns the largest magnitude among the eigenvalues of a
// small (2-by-2 or 3-by-3) matrix.
func SpectralRadius(A Mat) float64 {
	switch A.Rows() {
	case 1:
		return math.Abs(A[0][0])
	case 2:
		e := Eigenvalues2(A)
		return math.Max(cmplx.Abs(e[0]), cmplx.Abs(e[1]))
	case 3:
		e := Eigenvalues3(A)
		m := 0.0
		for _, l := range e {
			if a := cmplx.Abs(l); a > m {
				m = a
			}
		}
		return m
	default:
		return math.NaN()
	}
}

// solveCubic returns the three complex roots of a*x^3 + b*x^2 + c*x + d.
func solveCubic(a, b, c, d float64) [3]complex128 {
	if a == 0 {
		q := solveQuadratic(b, c, d)
		return [3]complex128{q[0], q[1], 0}
	}
	// Normalise and depress.
	B := b / a
	C := c / a
	D := d / a
	// x = t - B/3
	p := C - B*B/3
	q := 2*B*B*B/27 - B*C/3 + D
	shift := complex(-B/3, 0)
	// Roots of t^3 + p t + q = 0 via Cardano with complex arithmetic.
	pp := complex(p, 0)
	qq := complex(q, 0)
	disc := cmplx.Sqrt(qq*qq/4 + pp*pp*pp/27)
	u3 := -qq/2 + disc
	u := cmplx.Pow(u3, complex(1.0/3.0, 0))
	var roots [3]complex128
	if cmplx.Abs(u) < 1e-300 {
		// u ~ 0, use the other branch.
		v3 := -qq/2 - disc
		v := cmplx.Pow(v3, complex(1.0/3.0, 0))
		u = v
	}
	omega := complex(-0.5, math.Sqrt(3)/2)
	for k := 0; k < 3; k++ {
		uk := u
		for j := 0; j < k; j++ {
			uk *= omega
		}
		var vk complex128
		if cmplx.Abs(uk) > 1e-300 {
			vk = -pp / (3 * uk)
		}
		roots[k] = uk + vk + shift
	}
	return roots
}

// solveQuadratic returns the two complex roots of a*x^2 + b*x + c.
func solveQuadratic(a, b, c float64) [2]complex128 {
	if a == 0 {
		if b == 0 {
			return [2]complex128{0, 0}
		}
		return [2]complex128{complex(-c/b, 0), complex(-c/b, 0)}
	}
	disc := cmplx.Sqrt(complex(b*b-4*a*c, 0))
	da := complex(2*a, 0)
	return [2]complex128{(-complex(b, 0) + disc) / da, (-complex(b, 0) - disc) / da}
}
