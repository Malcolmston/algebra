package interval

import "errors"

// ErrDimension reports that the shapes of the interval matrices or vectors in
// an operation are incompatible.
var ErrDimension = errors.New("interval: incompatible dimensions")

// Matrix is a dense matrix whose entries are intervals. It is used for verified
// linear algebra: each entry encloses the corresponding real value, so matrix
// operations propagate rigorous bounds. The zero value is not usable; construct
// matrices with [NewMatrix], [MatrixFromFloat] or [Identity].
type Matrix struct {
	rows, cols int
	data       []Interval // row-major, length rows*cols
}

// NewMatrix returns a rows-by-cols interval matrix with every entry set to the
// degenerate interval [0, 0]. It panics if rows or cols is negative.
func NewMatrix(rows, cols int) *Matrix {
	if rows < 0 || cols < 0 {
		panic("interval: negative matrix dimension")
	}
	return &Matrix{rows: rows, cols: cols, data: make([]Interval, rows*cols)}
}

// MatrixFromFloat builds a rows-by-cols interval matrix whose entries are the
// degenerate (point) intervals of the given row-major values. It returns
// [ErrDimension] if len(values) != rows*cols.
func MatrixFromFloat(rows, cols int, values []float64) (*Matrix, error) {
	if len(values) != rows*cols {
		return nil, ErrDimension
	}
	m := NewMatrix(rows, cols)
	for i, v := range values {
		m.data[i] = Point(v)
	}
	return m, nil
}

// Identity returns the n-by-n interval identity matrix, with degenerate 1s on
// the diagonal and 0s elsewhere.
func Identity(n int) *Matrix {
	m := NewMatrix(n, n)
	for i := 0; i < n; i++ {
		m.data[i*n+i] = Point(1)
	}
	return m
}

// Rows returns the number of rows in the matrix.
func (m *Matrix) Rows() int { return m.rows }

// Cols returns the number of columns in the matrix.
func (m *Matrix) Cols() int { return m.cols }

// At returns the interval stored at row i, column j (both 0-indexed). It panics
// if the indices are out of range.
func (m *Matrix) At(i, j int) Interval {
	if i < 0 || i >= m.rows || j < 0 || j >= m.cols {
		panic("interval: matrix index out of range")
	}
	return m.data[i*m.cols+j]
}

// Set stores the interval v at row i, column j (both 0-indexed). It panics if
// the indices are out of range.
func (m *Matrix) Set(i, j int, v Interval) {
	if i < 0 || i >= m.rows || j < 0 || j >= m.cols {
		panic("interval: matrix index out of range")
	}
	m.data[i*m.cols+j] = v
}

// Clone returns a deep copy of the matrix that shares no storage with the
// receiver.
func (m *Matrix) Clone() *Matrix {
	out := NewMatrix(m.rows, m.cols)
	copy(out.data, m.data)
	return out
}

// Equal reports whether m and n have the same shape and every corresponding
// entry denotes the same interval.
func (m *Matrix) Equal(n *Matrix) bool {
	if m.rows != n.rows || m.cols != n.cols {
		return false
	}
	for i := range m.data {
		if !m.data[i].Equal(n.data[i]) {
			return false
		}
	}
	return true
}

// Add returns the entrywise sum m + n as a new matrix, each entry an outward
// rounded interval sum. It returns [ErrDimension] if the shapes differ.
func (m *Matrix) Add(n *Matrix) (*Matrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrDimension
	}
	out := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		out.data[i] = m.data[i].Add(n.data[i])
	}
	return out, nil
}

// Sub returns the entrywise difference m - n as a new matrix. It returns
// [ErrDimension] if the shapes differ.
func (m *Matrix) Sub(n *Matrix) (*Matrix, error) {
	if m.rows != n.rows || m.cols != n.cols {
		return nil, ErrDimension
	}
	out := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		out.data[i] = m.data[i].Sub(n.data[i])
	}
	return out, nil
}

// Scale returns a new matrix with every entry multiplied by the interval s.
func (m *Matrix) Scale(s Interval) *Matrix {
	out := NewMatrix(m.rows, m.cols)
	for i := range m.data {
		out.data[i] = m.data[i].Mul(s)
	}
	return out
}

// Mul returns the matrix product m * n as a new matrix. Each entry is an
// outward rounded interval dot product, so the result encloses the true product
// of any real matrices contained in m and n. It returns [ErrDimension] if the
// inner dimensions do not match.
func (m *Matrix) Mul(n *Matrix) (*Matrix, error) {
	if m.cols != n.rows {
		return nil, ErrDimension
	}
	out := NewMatrix(m.rows, n.cols)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < n.cols; j++ {
			sum := Point(0)
			for k := 0; k < m.cols; k++ {
				sum = sum.Add(m.data[i*m.cols+k].Mul(n.data[k*n.cols+j]))
			}
			out.data[i*n.cols+j] = sum
		}
	}
	return out, nil
}

// MulVec returns the matrix-vector product m * x as a new interval vector. It
// returns [ErrDimension] if len(x) != m.Cols().
func (m *Matrix) MulVec(x []Interval) ([]Interval, error) {
	if len(x) != m.cols {
		return nil, ErrDimension
	}
	out := make([]Interval, m.rows)
	for i := 0; i < m.rows; i++ {
		sum := Point(0)
		for k := 0; k < m.cols; k++ {
			sum = sum.Add(m.data[i*m.cols+k].Mul(x[k]))
		}
		out[i] = sum
	}
	return out, nil
}

// Transpose returns a new matrix that is the transpose of m.
func (m *Matrix) Transpose() *Matrix {
	out := NewMatrix(m.cols, m.rows)
	for i := 0; i < m.rows; i++ {
		for j := 0; j < m.cols; j++ {
			out.data[j*m.rows+i] = m.data[i*m.cols+j]
		}
	}
	return out
}

// MaxWidth returns the largest entry width in the matrix, a scalar measure of
// how tight the enclosure is. An empty matrix returns 0.
func (m *Matrix) MaxWidth() float64 {
	w := 0.0
	for i := range m.data {
		if d := m.data[i].Width(); d > w {
			w = d
		}
	}
	return w
}

// Midpoint returns the real matrix of entrywise midpoints as a row-major slice
// of length Rows*Cols.
func (m *Matrix) Midpoint() []float64 {
	out := make([]float64, len(m.data))
	for i := range m.data {
		out[i] = m.data[i].Midpoint()
	}
	return out
}

// DotVec returns an outward rounded enclosure of the dot product of interval
// vectors x and y. It returns [ErrDimension] if their lengths differ.
func DotVec(x, y []Interval) (Interval, error) {
	if len(x) != len(y) {
		return Empty(), ErrDimension
	}
	sum := Point(0)
	for i := range x {
		sum = sum.Add(x[i].Mul(y[i]))
	}
	return sum, nil
}
