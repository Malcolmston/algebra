package projectivegeom

import (
	"fmt"
	"math"
)

// Mat3 is a 3x3 real matrix in row-major order. It represents planar
// homographies, conics (as symmetric matrices) and general linear maps of
// homogeneous 3-vectors.
type Mat3 [3][3]float64

// NewMat3 builds a matrix from its nine entries given in row-major order.
func NewMat3(a, b, c, d, e, f, g, h, i float64) Mat3 {
	return Mat3{{a, b, c}, {d, e, f}, {g, h, i}}
}

// Identity3 returns the 3x3 identity matrix.
func Identity3() Mat3 { return Mat3{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} }

// Zero3 returns the 3x3 zero matrix.
func Zero3() Mat3 { return Mat3{} }

// Diag3 returns the diagonal matrix with entries a, b, c.
func Diag3(a, b, c float64) Mat3 { return Mat3{{a, 0, 0}, {0, b, 0}, {0, 0, c}} }

// At returns the entry in row i, column j.
func (m Mat3) At(i, j int) float64 { return m[i][j] }

// Row returns row i as a Vec3.
func (m Mat3) Row(i int) Vec3 { return Vec3{m[i][0], m[i][1], m[i][2]} }

// Col returns column j as a Vec3.
func (m Mat3) Col(j int) Vec3 { return Vec3{m[0][j], m[1][j], m[2][j]} }

// Add returns the entry-wise sum m+n.
func (m Mat3) Add(n Mat3) Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r[i][j] = m[i][j] + n[i][j]
		}
	}
	return r
}

// Sub returns the entry-wise difference m-n.
func (m Mat3) Sub(n Mat3) Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r[i][j] = m[i][j] - n[i][j]
		}
	}
	return r
}

// Scale returns m multiplied by the scalar s.
func (m Mat3) Scale(s float64) Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r[i][j] = m[i][j] * s
		}
	}
	return r
}

// Mul returns the matrix product m*n.
func (m Mat3) Mul(n Mat3) Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			var s float64
			for k := 0; k < 3; k++ {
				s += m[i][k] * n[k][j]
			}
			r[i][j] = s
		}
	}
	return r
}

// MulVec returns the matrix-vector product m*v.
func (m Mat3) MulVec(v Vec3) Vec3 {
	return Vec3{
		m[0][0]*v.X + m[0][1]*v.Y + m[0][2]*v.Z,
		m[1][0]*v.X + m[1][1]*v.Y + m[1][2]*v.Z,
		m[2][0]*v.X + m[2][1]*v.Y + m[2][2]*v.Z,
	}
}

// Transpose returns the transpose of m.
func (m Mat3) Transpose() Mat3 {
	var r Mat3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r[i][j] = m[j][i]
		}
	}
	return r
}

// Trace returns the sum of the diagonal entries.
func (m Mat3) Trace() float64 { return m[0][0] + m[1][1] + m[2][2] }

// Det returns the determinant of m.
func (m Mat3) Det() float64 {
	return det3(
		m[0][0], m[0][1], m[0][2],
		m[1][0], m[1][1], m[1][2],
		m[2][0], m[2][1], m[2][2],
	)
}

// Adjugate returns the classical adjugate (transpose of the cofactor matrix),
// so that m.Mul(m.Adjugate()) equals Det*Identity.
func (m Mat3) Adjugate() Mat3 {
	c := func(r0, r1, c0, c1 int) float64 {
		return m[r0][c0]*m[r1][c1] - m[r0][c1]*m[r1][c0]
	}
	var a Mat3
	a[0][0] = +c(1, 2, 1, 2)
	a[0][1] = -c(0, 2, 1, 2)
	a[0][2] = +c(0, 1, 1, 2)
	a[1][0] = -c(1, 2, 0, 2)
	a[1][1] = +c(0, 2, 0, 2)
	a[1][2] = -c(0, 1, 0, 2)
	a[2][0] = +c(1, 2, 0, 1)
	a[2][1] = -c(0, 2, 0, 1)
	a[2][2] = +c(0, 1, 0, 1)
	return a
}

// Inverse returns the inverse of m and true, or the zero matrix and false when
// m is (numerically) singular.
func (m Mat3) Inverse() (Mat3, bool) {
	d := m.Det()
	if math.Abs(d) < Eps*Eps {
		return Mat3{}, false
	}
	return m.Adjugate().Scale(1 / d), true
}

// IsSymmetric reports whether m equals its transpose within tol.
func (m Mat3) IsSymmetric(tol float64) bool {
	return math.Abs(m[0][1]-m[1][0]) <= tol &&
		math.Abs(m[0][2]-m[2][0]) <= tol &&
		math.Abs(m[1][2]-m[2][1]) <= tol
}

// Symmetrize returns (m + m^T)/2, the symmetric part of m.
func (m Mat3) Symmetrize() Mat3 { return m.Add(m.Transpose()).Scale(0.5) }

// ApproxEqual reports whether m and n agree entry-wise within tol.
func (m Mat3) ApproxEqual(n Mat3, tol float64) bool {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if math.Abs(m[i][j]-n[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

// MaxAbs returns the largest absolute entry of m.
func (m Mat3) MaxAbs() float64 {
	var mx float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if a := math.Abs(m[i][j]); a > mx {
				mx = a
			}
		}
	}
	return mx
}

// Quad returns the quadratic form u^T m v.
func (m Mat3) Quad(u, v Vec3) float64 { return u.Dot(m.MulVec(v)) }

// String renders m on three bracketed rows.
func (m Mat3) String() string {
	return fmt.Sprintf("[%g %g %g; %g %g %g; %g %g %g]",
		m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2])
}

// Mat4 is a 4x4 real matrix in row-major order. It represents collineations of
// RP^3 and quadric matrices.
type Mat4 [4][4]float64

// Identity4 returns the 4x4 identity matrix.
func Identity4() Mat4 {
	var m Mat4
	for i := 0; i < 4; i++ {
		m[i][i] = 1
	}
	return m
}

// At returns the entry in row i, column j.
func (m Mat4) At(i, j int) float64 { return m[i][j] }

// Mul returns the matrix product m*n.
func (m Mat4) Mul(n Mat4) Mat4 {
	var r Mat4
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			var s float64
			for k := 0; k < 4; k++ {
				s += m[i][k] * n[k][j]
			}
			r[i][j] = s
		}
	}
	return r
}

// MulVec returns the matrix-vector product m*v.
func (m Mat4) MulVec(v Vec4) Vec4 {
	c := [4]float64{v.X, v.Y, v.Z, v.W}
	var o [4]float64
	for i := 0; i < 4; i++ {
		var s float64
		for k := 0; k < 4; k++ {
			s += m[i][k] * c[k]
		}
		o[i] = s
	}
	return Vec4{o[0], o[1], o[2], o[3]}
}

// Transpose returns the transpose of m.
func (m Mat4) Transpose() Mat4 {
	var r Mat4
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			r[i][j] = m[j][i]
		}
	}
	return r
}

// Trace returns the sum of the diagonal entries.
func (m Mat4) Trace() float64 { return m[0][0] + m[1][1] + m[2][2] + m[3][3] }

// Det returns the determinant of m.
func (m Mat4) Det() float64 {
	d, _ := detN(m.dense())
	return d
}

// Inverse returns the inverse of m and true, or the zero matrix and false when
// m is (numerically) singular.
func (m Mat4) Inverse() (Mat4, bool) {
	inv, ok := invertDense(m.dense())
	if !ok {
		return Mat4{}, false
	}
	var r Mat4
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			r[i][j] = inv[i][j]
		}
	}
	return r, true
}

// Quad returns the quadratic form u^T m v.
func (m Mat4) Quad(u, v Vec4) float64 { return u.Dot(m.MulVec(v)) }

// ApproxEqual reports whether m and n agree entry-wise within tol.
func (m Mat4) ApproxEqual(n Mat4, tol float64) bool {
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if math.Abs(m[i][j]-n[i][j]) > tol {
				return false
			}
		}
	}
	return true
}

func (m Mat4) dense() [][]float64 {
	d := make([][]float64, 4)
	for i := 0; i < 4; i++ {
		d[i] = []float64{m[i][0], m[i][1], m[i][2], m[i][3]}
	}
	return d
}

// detN returns the determinant of a square dense matrix via partial-pivot
// Gaussian elimination, together with false when the matrix is empty or
// non-square.
func detN(a [][]float64) (float64, bool) {
	n := len(a)
	if n == 0 {
		return 0, false
	}
	m := make([][]float64, n)
	for i := range a {
		if len(a[i]) != n {
			return 0, false
		}
		m[i] = append([]float64(nil), a[i]...)
	}
	det := 1.0
	for col := 0; col < n; col++ {
		piv := col
		best := math.Abs(m[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(m[r][col]); v > best {
				best, piv = v, r
			}
		}
		if best < 1e-300 {
			return 0, true
		}
		if piv != col {
			m[col], m[piv] = m[piv], m[col]
			det = -det
		}
		det *= m[col][col]
		inv := 1 / m[col][col]
		for r := col + 1; r < n; r++ {
			f := m[r][col] * inv
			if f == 0 {
				continue
			}
			for c := col; c < n; c++ {
				m[r][c] -= f * m[col][c]
			}
		}
	}
	return det, true
}

// invertDense returns the inverse of a square dense matrix via Gauss-Jordan
// elimination with partial pivoting, or false when the matrix is singular.
func invertDense(a [][]float64) ([][]float64, bool) {
	n := len(a)
	aug := make([][]float64, n)
	for i := 0; i < n; i++ {
		aug[i] = make([]float64, 2*n)
		copy(aug[i], a[i])
		aug[i][n+i] = 1
	}
	for col := 0; col < n; col++ {
		piv := col
		best := math.Abs(aug[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(aug[r][col]); v > best {
				best, piv = v, r
			}
		}
		if best < 1e-12 {
			return nil, false
		}
		aug[col], aug[piv] = aug[piv], aug[col]
		inv := 1 / aug[col][col]
		for c := 0; c < 2*n; c++ {
			aug[col][c] *= inv
		}
		for r := 0; r < n; r++ {
			if r == col {
				continue
			}
			f := aug[r][col]
			if f == 0 {
				continue
			}
			for c := 0; c < 2*n; c++ {
				aug[r][c] -= f * aug[col][c]
			}
		}
	}
	out := make([][]float64, n)
	for i := 0; i < n; i++ {
		out[i] = append([]float64(nil), aug[i][n:]...)
	}
	return out, true
}

// nullVector returns a non-zero solution x of the homogeneous system a*x=0 for
// a matrix with rows count one less than the number of columns (so the null
// space is generically one-dimensional). It uses Gaussian elimination and back
// substitution over the pivot structure, returning false only if the input is
// malformed.
func nullVector(a [][]float64) ([]float64, bool) {
	rows := len(a)
	if rows == 0 {
		return nil, false
	}
	cols := len(a[0])
	m := make([][]float64, rows)
	for i := range a {
		if len(a[i]) != cols {
			return nil, false
		}
		m[i] = append([]float64(nil), a[i]...)
	}
	pivotCol := make([]int, 0, rows)
	r := 0
	for c := 0; c < cols && r < rows; c++ {
		piv := -1
		best := 1e-12
		for i := r; i < rows; i++ {
			if v := math.Abs(m[i][c]); v > best {
				best, piv = v, i
			}
		}
		if piv < 0 {
			continue
		}
		m[r], m[piv] = m[piv], m[r]
		inv := 1 / m[r][c]
		for j := 0; j < cols; j++ {
			m[r][j] *= inv
		}
		for i := 0; i < rows; i++ {
			if i == r {
				continue
			}
			f := m[i][c]
			if f == 0 {
				continue
			}
			for j := 0; j < cols; j++ {
				m[i][j] -= f * m[r][j]
			}
		}
		pivotCol = append(pivotCol, c)
		r++
	}
	// Find a free column (one without a pivot); set it to 1 and back out the
	// pivot variables.
	isPivot := make([]bool, cols)
	for _, c := range pivotCol {
		isPivot[c] = true
	}
	free := -1
	for c := 0; c < cols; c++ {
		if !isPivot[c] {
			free = c
			break
		}
	}
	x := make([]float64, cols)
	if free < 0 {
		// Fully determined up to scale only when rank == cols-1 is unmet; fall
		// back to the last column direction.
		free = cols - 1
	}
	x[free] = 1
	for i, c := range pivotCol {
		// Row i in reduced form: x[c] + sum over free cols m[i][f]*x[f] = 0.
		x[c] = -m[i][free]
	}
	// Guard against an all-zero result.
	allZero := true
	for _, v := range x {
		if math.Abs(v) > 1e-300 {
			allZero = false
			break
		}
	}
	if allZero {
		x[cols-1] = 1
	}
	return x, true
}
