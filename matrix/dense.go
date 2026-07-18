package matrix

import "fmt"

// matrixDenseBlock is the cache-blocking tile edge used by the tiled multiply
// and transpose kernels. Sub-blocks of this size are chosen so that the working
// set of the inner loops stays resident in a typical L1/L2 cache while keeping
// the innermost loop long enough to amortize loop overhead. 64 is a good
// general-purpose default for float64 data.
const matrixDenseBlock = 64

// Dense is a performance-focused, flat row-major matrix of float64 values. It
// deliberately avoids the per-element algebra.Expr allocation and the [][]Expr
// row-of-slices indirection of the symbolic [Matrix] core: every element lives
// in a single contiguous backing slice with a fixed stride equal to the column
// count, so the numeric hot paths (multiply, AXPY, scale, transpose) run with
// unit-stride memory access and allocate nothing once their buffers are sized.
//
// The zero value is not usable; construct a Dense with [NewDense],
// [Matrix.ToDense] or the buffers produced by the in-place operations. Element
// (i, j) is stored at data[i*cols+j]. All operations are deterministic.
type Dense struct {
	rows, cols int
	data       []float64
}

// NewDense returns a rows×cols dense matrix with every element initialized to
// zero. It panics if rows or cols is negative. The backing slice is a single
// contiguous row-major buffer of length rows*cols with stride cols.
func NewDense(rows, cols int) *Dense {
	if rows < 0 || cols < 0 {
		panic("matrix: negative dimension")
	}
	return &Dense{rows: rows, cols: cols, data: make([]float64, rows*cols)}
}

// Rows returns the number of rows.
func (d *Dense) Rows() int { return d.rows }

// Cols returns the number of columns.
func (d *Dense) Cols() int { return d.cols }

// At returns the element at row i, column j (both 0-based). It panics if the
// indices are out of range.
func (d *Dense) At(i, j int) float64 {
	if i < 0 || i >= d.rows || j < 0 || j >= d.cols {
		panic("matrix: index out of range")
	}
	return d.data[i*d.cols+j]
}

// Set stores v at row i, column j (both 0-based). It panics if the indices are
// out of range.
func (d *Dense) Set(i, j int, v float64) {
	if i < 0 || i >= d.rows || j < 0 || j >= d.cols {
		panic("matrix: index out of range")
	}
	d.data[i*d.cols+j] = v
}

// ToDense converts the symbolic matrix to a numeric [Dense] using the numeric
// fast-path of [Matrix.Floats]. It returns [ErrUnsupported] wrapping the
// underlying evaluation error when any entry contains a free symbol or
// otherwise cannot be evaluated to a float64. The returned Dense owns a fresh
// contiguous buffer and shares no storage with the receiver.
func (m *Matrix) ToDense() (*Dense, error) {
	ff, err := m.Floats()
	if err != nil {
		return nil, fmt.Errorf("%w: %v", ErrUnsupported, err)
	}
	d := NewDense(m.rows, m.cols)
	for i := 0; i < m.rows; i++ {
		copy(d.data[i*d.cols:i*d.cols+d.cols], ff[i])
	}
	return d, nil
}

// ToMatrix converts the dense matrix back into a symbolic [Matrix] whose entries
// are inexact algebra.Flt literals, via [FromFloats]. The result shares no
// storage with the receiver.
func (d *Dense) ToMatrix() *Matrix {
	vals := make([][]float64, d.rows)
	for i := 0; i < d.rows; i++ {
		vals[i] = make([]float64, d.cols)
		copy(vals[i], d.data[i*d.cols:i*d.cols+d.cols])
	}
	return FromFloats(vals)
}

// matrixDenseEnsure resizes d to the given shape, reusing d's backing slice when
// it already has enough capacity and allocating a fresh, zeroed buffer only when
// it does not. It is the mechanism that lets the in-place kernels allocate
// nothing in steady state. The (re)used buffer is zeroed so callers may
// accumulate into it.
func matrixDenseEnsure(d *Dense, rows, cols int) {
	n := rows * cols
	if cap(d.data) >= n {
		d.data = d.data[:n]
		for i := range d.data {
			d.data[i] = 0
		}
	} else {
		d.data = make([]float64, n)
	}
	d.rows = rows
	d.cols = cols
}

// MulInto computes the matrix product a·b and writes it into the receiver,
// returning the receiver for convenience. It returns [ErrDimension] if the inner
// dimensions disagree (a.Cols() != b.Rows()). The receiver's backing slice is
// reused in place, grown only when its capacity is insufficient, so repeated
// calls with compatibly shaped operands allocate nothing.
//
// The kernel uses an i-k-j, cache-blocked (tiled) loop order with a tile edge of
// matrixDenseBlock. Iterating k in the middle and j innermost makes the inner
// loop a unit-stride AXPY over contiguous rows of a, b and the result buffer,
// which maximizes sequential memory throughput and lets the compiler keep the
// scalar a[i,k] in a register across the row. The receiver may alias neither a
// nor b; pass a distinct destination.
func (d *Dense) MulInto(a, b *Dense) (*Dense, error) {
	if a.cols != b.rows {
		return nil, ErrDimension
	}
	m, k, n := a.rows, a.cols, b.cols
	if d == a || d == b {
		panic("matrix: MulInto destination must not alias an operand")
	}
	matrixDenseEnsure(d, m, n)
	bl := matrixDenseBlock
	for ii := 0; ii < m; ii += bl {
		iMax := ii + bl
		if iMax > m {
			iMax = m
		}
		for kk := 0; kk < k; kk += bl {
			kMax := kk + bl
			if kMax > k {
				kMax = k
			}
			for jj := 0; jj < n; jj += bl {
				jMax := jj + bl
				if jMax > n {
					jMax = n
				}
				for i := ii; i < iMax; i++ {
					crow := d.data[i*n : i*n+n]
					arow := a.data[i*k : i*k+k]
					for kp := kk; kp < kMax; kp++ {
						av := arow[kp]
						if av == 0 {
							continue
						}
						brow := b.data[kp*n : kp*n+n]
						for j := jj; j < jMax; j++ {
							crow[j] += av * brow[j]
						}
					}
				}
			}
		}
	}
	return d, nil
}

// Scale multiplies every element of the receiver by s in place and returns the
// receiver. It allocates nothing.
func (d *Dense) Scale(s float64) *Dense {
	for i := range d.data {
		d.data[i] *= s
	}
	return d
}

// AddScaled performs the in-place AXPY update d += s*o and returns the receiver.
// It returns [ErrDimension] if o does not have the same shape as the receiver.
// The update walks the two contiguous backing slices with unit stride and
// allocates nothing.
func (d *Dense) AddScaled(o *Dense, s float64) (*Dense, error) {
	if d.rows != o.rows || d.cols != o.cols {
		return nil, ErrDimension
	}
	dd := d.data
	od := o.data
	for i := range dd {
		dd[i] += s * od[i]
	}
	return d, nil
}

// TransposeInto writes the transpose of the receiver into dst and returns dst.
// The destination's backing slice is reused in place, grown only when its
// capacity is insufficient, so repeated calls allocate nothing. dst must not
// alias the receiver.
//
// The copy is cache-blocked with a tile edge of matrixDenseBlock: transposing
// one small tile at a time keeps both the source and destination tiles resident
// in cache, avoiding the cache-thrashing strided access pattern of a naive
// element-by-element transpose on large matrices.
func (d *Dense) TransposeInto(dst *Dense) *Dense {
	if dst == d {
		panic("matrix: TransposeInto destination must not alias the source")
	}
	m, n := d.rows, d.cols
	matrixDenseEnsure(dst, n, m)
	bl := matrixDenseBlock
	for ii := 0; ii < m; ii += bl {
		iMax := ii + bl
		if iMax > m {
			iMax = m
		}
		for jj := 0; jj < n; jj += bl {
			jMax := jj + bl
			if jMax > n {
				jMax = n
			}
			for i := ii; i < iMax; i++ {
				srow := d.data[i*n : i*n+n]
				for j := jj; j < jMax; j++ {
					dst.data[j*m+i] = srow[j]
				}
			}
		}
	}
	return dst
}
