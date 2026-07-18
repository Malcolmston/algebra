package matrix

import (
	"errors"
	"math"
	"testing"

	"github.com/malcolmston/algebra"
)

// denseFromRows builds a Dense from a rectangular literal for test setup.
func denseFromRows(rows [][]float64) *Dense {
	r := len(rows)
	c := 0
	if r > 0 {
		c = len(rows[0])
	}
	d := NewDense(r, c)
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			d.Set(i, j, rows[i][j])
		}
	}
	return d
}

// denseEqual reports whether two Dense matrices match within tol.
func denseEqual(a, b *Dense, tol float64) bool {
	if a.Rows() != b.Rows() || a.Cols() != b.Cols() {
		return false
	}
	for i := 0; i < a.Rows(); i++ {
		for j := 0; j < a.Cols(); j++ {
			if math.Abs(a.At(i, j)-b.At(i, j)) > tol {
				return false
			}
		}
	}
	return true
}

func TestNewDenseZeroAndAccess(t *testing.T) {
	d := NewDense(2, 3)
	if d.Rows() != 2 || d.Cols() != 3 {
		t.Fatalf("shape = %dx%d, want 2x3", d.Rows(), d.Cols())
	}
	for i := 0; i < 2; i++ {
		for j := 0; j < 3; j++ {
			if d.At(i, j) != 0 {
				t.Fatalf("At(%d,%d) = %v, want 0", i, j, d.At(i, j))
			}
		}
	}
	d.Set(1, 2, 7.5)
	if d.At(1, 2) != 7.5 {
		t.Fatalf("At(1,2) = %v, want 7.5", d.At(1, 2))
	}
	// Confirm contiguous row-major layout: element (1,2) is at index 1*3+2.
	if d.data[1*3+2] != 7.5 {
		t.Fatalf("backing layout mismatch: data[5] = %v, want 7.5", d.data[5])
	}
}

func TestNewDenseNegativePanics(t *testing.T) {
	defer func() {
		if recover() == nil {
			t.Fatalf("expected panic on negative dimension")
		}
	}()
	NewDense(-1, 2)
}

func TestDenseRoundTripMatrix(t *testing.T) {
	src := FromFloats([][]float64{{1, 2, 3}, {4, 5, 6}})
	d, err := src.ToDense()
	if err != nil {
		t.Fatalf("ToDense: %v", err)
	}
	if d.At(0, 0) != 1 || d.At(1, 2) != 6 {
		t.Fatalf("ToDense values wrong: %v %v", d.At(0, 0), d.At(1, 2))
	}
	back := d.ToMatrix()
	got, err := back.Floats()
	if err != nil {
		t.Fatalf("Floats: %v", err)
	}
	want := [][]float64{{1, 2, 3}, {4, 5, 6}}
	for i := range want {
		for j := range want[i] {
			if got[i][j] != want[i][j] {
				t.Fatalf("round trip (%d,%d) = %v, want %v", i, j, got[i][j], want[i][j])
			}
		}
	}
}

func TestToDenseSymbolicUnsupported(t *testing.T) {
	m := New(1, 1)
	m.Set(0, 0, algebra.Sym("x"))
	_, err := m.ToDense()
	if !errors.Is(err, ErrUnsupported) {
		t.Fatalf("err = %v, want ErrUnsupported", err)
	}
}

func TestMulIntoKnownAnswers(t *testing.T) {
	tests := []struct {
		name string
		a, b [][]float64
		want [][]float64
	}{
		{
			name: "2x2 identity",
			a:    [][]float64{{1, 0}, {0, 1}},
			b:    [][]float64{{3, 4}, {5, 6}},
			want: [][]float64{{3, 4}, {5, 6}},
		},
		{
			name: "2x3 times 3x2",
			a:    [][]float64{{1, 2, 3}, {4, 5, 6}},
			b:    [][]float64{{7, 8}, {9, 10}, {11, 12}},
			want: [][]float64{{58, 64}, {139, 154}},
		},
		{
			name: "1x1",
			a:    [][]float64{{2}},
			b:    [][]float64{{3}},
			want: [][]float64{{6}},
		},
		{
			name: "negatives",
			a:    [][]float64{{-1, 2}, {3, -4}},
			b:    [][]float64{{5, -6}, {-7, 8}},
			want: [][]float64{{-19, 22}, {43, -50}},
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			a := denseFromRows(tc.a)
			b := denseFromRows(tc.b)
			want := denseFromRows(tc.want)
			out := NewDense(0, 0)
			got, err := out.MulInto(a, b)
			if err != nil {
				t.Fatalf("MulInto: %v", err)
			}
			if !denseEqual(got, want, 1e-12) {
				t.Fatalf("MulInto result mismatch\n got %v\nwant %v", got.data, want.data)
			}
		})
	}
}

func TestMulIntoLargeAgainstNaive(t *testing.T) {
	// Deterministic pseudo-random-ish fill via a fixed formula, then compare the
	// blocked kernel against a straightforward triple loop.
	const m, k, n = 40, 33, 27
	a := NewDense(m, k)
	b := NewDense(k, n)
	for i := 0; i < m; i++ {
		for j := 0; j < k; j++ {
			a.Set(i, j, float64((i*7+j*3)%11)-5)
		}
	}
	for i := 0; i < k; i++ {
		for j := 0; j < n; j++ {
			b.Set(i, j, float64((i*5+j*2)%13)-6)
		}
	}
	naive := NewDense(m, n)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			var s float64
			for kp := 0; kp < k; kp++ {
				s += a.At(i, kp) * b.At(kp, j)
			}
			naive.Set(i, j, s)
		}
	}
	out := NewDense(0, 0)
	got, err := out.MulInto(a, b)
	if err != nil {
		t.Fatalf("MulInto: %v", err)
	}
	if !denseEqual(got, naive, 1e-9) {
		t.Fatalf("blocked multiply disagrees with naive")
	}
}

func TestMulIntoDimensionError(t *testing.T) {
	a := denseFromRows([][]float64{{1, 2, 3}})
	b := denseFromRows([][]float64{{1, 2}})
	out := NewDense(0, 0)
	if _, err := out.MulInto(a, b); !errors.Is(err, ErrDimension) {
		t.Fatalf("err = %v, want ErrDimension", err)
	}
}

func TestMulIntoReusesBuffer(t *testing.T) {
	a := denseFromRows([][]float64{{1, 2}, {3, 4}})
	b := denseFromRows([][]float64{{5, 6}, {7, 8}})
	out := NewDense(2, 2)
	if _, err := out.MulInto(a, b); err != nil {
		t.Fatalf("MulInto: %v", err)
	}
	first := &out.data[0]
	firstCap := cap(out.data)
	// Second call with the same shapes must not reallocate.
	if _, err := out.MulInto(a, b); err != nil {
		t.Fatalf("MulInto: %v", err)
	}
	if &out.data[0] != first || cap(out.data) != firstCap {
		t.Fatalf("MulInto reallocated backing slice on reuse")
	}
	want := denseFromRows([][]float64{{19, 22}, {43, 50}})
	if !denseEqual(out, want, 1e-12) {
		t.Fatalf("reused result wrong: %v", out.data)
	}
}

func TestScaleInPlace(t *testing.T) {
	d := denseFromRows([][]float64{{1, -2}, {3, 4}})
	before := &d.data[0]
	got := d.Scale(2)
	if got != d {
		t.Fatalf("Scale did not return receiver")
	}
	if &d.data[0] != before {
		t.Fatalf("Scale reallocated")
	}
	want := denseFromRows([][]float64{{2, -4}, {6, 8}})
	if !denseEqual(d, want, 1e-12) {
		t.Fatalf("Scale wrong: %v", d.data)
	}
}

func TestAddScaledAXPY(t *testing.T) {
	d := denseFromRows([][]float64{{1, 2}, {3, 4}})
	o := denseFromRows([][]float64{{10, 20}, {30, 40}})
	before := &d.data[0]
	got, err := d.AddScaled(o, 0.5)
	if err != nil {
		t.Fatalf("AddScaled: %v", err)
	}
	if got != d || &d.data[0] != before {
		t.Fatalf("AddScaled must update in place and return receiver")
	}
	want := denseFromRows([][]float64{{6, 12}, {18, 24}})
	if !denseEqual(d, want, 1e-12) {
		t.Fatalf("AddScaled wrong: %v", d.data)
	}
}

func TestAddScaledDimensionError(t *testing.T) {
	d := denseFromRows([][]float64{{1, 2}})
	o := denseFromRows([][]float64{{1, 2, 3}})
	if _, err := d.AddScaled(o, 1); !errors.Is(err, ErrDimension) {
		t.Fatalf("err = %v, want ErrDimension", err)
	}
}

func TestTransposeInto(t *testing.T) {
	tests := [][][]float64{
		{{1, 2, 3}, {4, 5, 6}},
		{{1}},
		{{1, 2}, {3, 4}, {5, 6}},
	}
	for _, rows := range tests {
		src := denseFromRows(rows)
		dst := NewDense(0, 0)
		got := src.TransposeInto(dst)
		if got != dst {
			t.Fatalf("TransposeInto did not return dst")
		}
		if got.Rows() != src.Cols() || got.Cols() != src.Rows() {
			t.Fatalf("transpose shape %dx%d, want %dx%d", got.Rows(), got.Cols(), src.Cols(), src.Rows())
		}
		for i := 0; i < src.Rows(); i++ {
			for j := 0; j < src.Cols(); j++ {
				if got.At(j, i) != src.At(i, j) {
					t.Fatalf("transpose (%d,%d) mismatch", i, j)
				}
			}
		}
	}
}

func TestTransposeIntoLargeBlocked(t *testing.T) {
	const m, n = 70, 45
	src := NewDense(m, n)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			src.Set(i, j, float64(i*100+j))
		}
	}
	dst := NewDense(0, 0)
	src.TransposeInto(dst)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			if dst.At(j, i) != src.At(i, j) {
				t.Fatalf("blocked transpose (%d,%d) mismatch", i, j)
			}
		}
	}
}

func TestTransposeIntoReusesBuffer(t *testing.T) {
	src := denseFromRows([][]float64{{1, 2, 3}, {4, 5, 6}})
	dst := NewDense(3, 2)
	src.TransposeInto(dst)
	first := &dst.data[0]
	firstCap := cap(dst.data)
	src.TransposeInto(dst)
	if &dst.data[0] != first || cap(dst.data) != firstCap {
		t.Fatalf("TransposeInto reallocated on reuse")
	}
}

// denseBench builds a deterministically filled n×n Dense for benchmarks.
func denseBench(n int) *Dense {
	d := NewDense(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			d.Set(i, j, float64((i*31+j*17)%97)-48)
		}
	}
	return d
}

func BenchmarkMulInto(b *testing.B) {
	const n = 128
	a := denseBench(n)
	c := denseBench(n)
	out := NewDense(n, n)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := out.MulInto(a, c); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkAddScaled(b *testing.B) {
	const n = 256
	d := denseBench(n)
	o := denseBench(n)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := d.AddScaled(o, 1.0000001); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkScale(b *testing.B) {
	const n = 256
	d := denseBench(n)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		d.Scale(1.0000001)
	}
}

func BenchmarkTransposeInto(b *testing.B) {
	const n = 256
	src := denseBench(n)
	dst := NewDense(n, n)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		src.TransposeInto(dst)
	}
}
