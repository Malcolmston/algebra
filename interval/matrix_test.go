package interval

import (
	"math"
	"testing"
)

func TestMatrixBasics(t *testing.T) {
	m, err := MatrixFromFloat(2, 2, []float64{1, 2, 3, 4})
	if err != nil {
		t.Fatal(err)
	}
	if m.Rows() != 2 || m.Cols() != 2 {
		t.Errorf("shape = %dx%d, want 2x2", m.Rows(), m.Cols())
	}
	if !m.At(0, 1).Equal(Point(2)) {
		t.Errorf("At(0,1) = %v, want [2,2]", m.At(0, 1))
	}
	m.Set(0, 0, New(0, 1))
	if !m.At(0, 0).Equal(New(0, 1)) {
		t.Error("Set/At mismatch")
	}
	if !m.Clone().Equal(m) {
		t.Error("Clone should equal original")
	}
	if _, err := MatrixFromFloat(2, 2, []float64{1}); err == nil {
		t.Error("MatrixFromFloat should reject wrong length")
	}
}

func TestMatrixAddMul(t *testing.T) {
	a, _ := MatrixFromFloat(2, 2, []float64{1, 2, 3, 4})
	b, _ := MatrixFromFloat(2, 2, []float64{5, 6, 7, 8})
	sum, err := a.Add(b)
	if err != nil {
		t.Fatal(err)
	}
	if !sum.At(0, 0).Contains(6) || !sum.At(1, 1).Contains(12) {
		t.Errorf("Add wrong: %v %v", sum.At(0, 0), sum.At(1, 1))
	}
	// Product [[1,2],[3,4]] * [[5,6],[7,8]] = [[19,22],[43,50]].
	prod, err := a.Mul(b)
	if err != nil {
		t.Fatal(err)
	}
	want := []float64{19, 22, 43, 50}
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			if !prod.At(i, j).Contains(want[i*2+j]) {
				t.Errorf("Mul[%d][%d]=%v want enclose %g", i, j, prod.At(i, j), want[i*2+j])
			}
		}
	}
	// Dimension mismatch.
	c := NewMatrix(3, 2)
	if _, err := a.Mul(c); err == nil {
		t.Error("Mul should reject inner-dim mismatch")
	}
}

func TestMatrixVecAndTranspose(t *testing.T) {
	m, _ := MatrixFromFloat(2, 3, []float64{1, 2, 3, 4, 5, 6})
	x := []Interval{Point(1), Point(0), Point(-1)}
	// [1,2,3]*[1,0,-1] = -2 ; [4,5,6]*[1,0,-1] = -2.
	y, err := m.MulVec(x)
	if err != nil {
		t.Fatal(err)
	}
	if !y[0].Contains(-2) || !y[1].Contains(-2) {
		t.Errorf("MulVec = %v %v, want enclose -2,-2", y[0], y[1])
	}
	tr := m.Transpose()
	if tr.Rows() != 3 || tr.Cols() != 2 || !tr.At(2, 1).Contains(6) {
		t.Errorf("Transpose wrong: %dx%d At(2,1)=%v", tr.Rows(), tr.Cols(), tr.At(2, 1))
	}
}

func TestIdentityScaleWidth(t *testing.T) {
	id := Identity(3)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			exp := 0.0
			if i == j {
				exp = 1
			}
			if !id.At(i, j).Equal(Point(exp)) {
				t.Errorf("Identity[%d][%d]=%v want %g", i, j, id.At(i, j), exp)
			}
		}
	}
	// Identity * A = A (enclosure).
	a, _ := MatrixFromFloat(3, 1, []float64{7, 8, 9})
	got, _ := id.Mul(a)
	for i := 0; i < 3; i++ {
		if !got.At(i, 0).Contains(a.At(i, 0).Midpoint()) {
			t.Errorf("I*A row %d = %v", i, got.At(i, 0))
		}
	}
	s := NewMatrix(1, 1)
	s.Set(0, 0, New(1, 3))
	if w := s.MaxWidth(); math.Abs(w-2) > 1e-9 {
		t.Errorf("MaxWidth = %g, want 2", w)
	}
	if md := s.Midpoint(); math.Abs(md[0]-2) > 1e-12 {
		t.Errorf("Midpoint = %g, want 2", md[0])
	}
}

func TestDotVec(t *testing.T) {
	x := []Interval{Point(1), Point(2), Point(3)}
	y := []Interval{Point(4), Point(5), Point(6)}
	got, err := DotVec(x, y)
	if err != nil {
		t.Fatal(err)
	}
	encloses(t, "DotVec", got, 32, 1e-9) // 4+10+18
	if _, err := DotVec(x, y[:2]); err == nil {
		t.Error("DotVec should reject length mismatch")
	}
}

// BenchmarkMatrixMul exercises the heaviest routine, dense interval
// matrix-matrix multiplication, which is cubic in the dimension and performs an
// outward-rounded interval multiply-add at each of the O(n^3) inner steps.
func BenchmarkMatrixMul(b *testing.B) {
	const n = 32
	vals := make([]float64, n*n)
	for i := range vals {
		vals[i] = float64(i%7) - 3
	}
	m, _ := MatrixFromFloat(n, n, vals)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := m.Mul(m); err != nil {
			b.Fatal(err)
		}
	}
}
