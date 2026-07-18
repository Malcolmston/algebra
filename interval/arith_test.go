package interval

import (
	"math"
	"testing"
)

func TestAddSubEnclosure(t *testing.T) {
	tests := []struct {
		a, b       Interval
		wantLo, hi float64
	}{
		{New(1, 2), New(3, 4), 4, 6},
		{New(-1, 1), New(-2, 3), -3, 4},
		{New(0.1, 0.1), New(0.2, 0.2), 0.3, 0.3},
	}
	for _, tc := range tests {
		got := tc.a.Add(tc.b)
		if !got.Contains(tc.wantLo) || !got.Contains(tc.hi) {
			t.Errorf("Add(%v,%v) = %v, does not enclose [%g,%g]", tc.a, tc.b, got, tc.wantLo, tc.hi)
		}
		// Width should stay tight (a couple ULPs of slack only).
		if got.Width() > (tc.hi-tc.wantLo)+1e-9 {
			t.Errorf("Add(%v,%v) width %g too large", tc.a, tc.b, got.Width())
		}
	}
	// Sub: [1,2]-[3,4] = [-3,-1].
	got := New(1, 2).Sub(New(3, 4))
	if !got.Contains(-3) || !got.Contains(-1) {
		t.Errorf("Sub = %v, want enclose [-3,-1]", got)
	}
	// Empty propagation.
	if !New(1, 2).Add(Empty()).IsEmpty() {
		t.Error("Add with empty should be empty")
	}
}

func TestMulEnclosure(t *testing.T) {
	tests := []struct {
		a, b       Interval
		wantLo, hi float64
	}{
		{New(2, 3), New(4, 5), 8, 15},
		{New(-2, 3), New(-4, 5), -12, 15}, // min(-2*5,3*-4)=-12 ... check: products {8,-10,-12,15}
		{New(-3, -1), New(-2, -1), 1, 6},
		{New(-1, 1), New(-1, 1), -1, 1},
	}
	for _, tc := range tests {
		got := tc.a.Mul(tc.b)
		if !got.Contains(tc.wantLo) || !got.Contains(tc.hi) {
			t.Errorf("Mul(%v,%v) = %v, want enclose [%g,%g]", tc.a, tc.b, got, tc.wantLo, tc.hi)
		}
		if got.Width() > (tc.hi-tc.wantLo)+1e-9 {
			t.Errorf("Mul(%v,%v) width %g too large", tc.a, tc.b, got.Width())
		}
	}
}

func TestDivAndRecip(t *testing.T) {
	// [10,20]/[2,4] = [2.5,10].
	got := New(10, 20).Div(New(2, 4))
	if !got.Contains(2.5) || !got.Contains(10) {
		t.Errorf("Div = %v, want enclose [2.5,10]", got)
	}
	// point/point stays tight.
	encloses(t, "Div point", New(6, 6).Div(New(2, 2)), 3, 1e-12)
	// Reciprocal of [2,4] = [0.25,0.5].
	r := New(2, 4).Recip()
	if !r.Contains(0.25) || !r.Contains(0.5) {
		t.Errorf("Recip = %v, want enclose [0.25,0.5]", r)
	}
	// Division by interval containing 0 gives entire.
	if !New(1, 2).Div(New(-1, 1)).IsEntire() {
		t.Error("Div by zero-straddling interval should be entire")
	}
	// Reciprocal of [0,0] is empty.
	if !Point(0).Recip().IsEmpty() {
		t.Error("Recip of [0,0] should be empty")
	}
	// One-sided: 1/[0,2] = [0.5, +Inf).
	half := New(0, 2).Recip()
	if !half.Contains(0.5) || !math.IsInf(half.Hi, 1) {
		t.Errorf("Recip([0,2]) = %v, want [0.5,+Inf)", half)
	}
}

func TestNegAbsSquare(t *testing.T) {
	if got := New(2, 5).Neg(); !got.Equal(New(-5, -2)) {
		t.Errorf("Neg = %v, want [-5,-2]", got)
	}
	if got := New(-3, 2).Abs(); !got.Equal(New(0, 3)) {
		t.Errorf("Abs = %v, want [0,3]", got)
	}
	// Square of [-3,2] = [0,9] (tighter than Mul).
	sq := New(-3, 2).Square()
	if !sq.Contains(0) || !sq.Contains(9) || sq.Lo < 0 {
		t.Errorf("Square = %v, want enclose [0,9]", sq)
	}
	// Square of [2,3] = [4,9].
	sq2 := New(2, 3).Square()
	if !sq2.Contains(4) || !sq2.Contains(9) {
		t.Errorf("Square([2,3]) = %v, want enclose [4,9]", sq2)
	}
	// Square of a point stays tight.
	encloses(t, "Square point", Point(3).Square(), 9, 1e-9)
}

func TestMinMaxScalarOps(t *testing.T) {
	if got := New(1, 4).Min(New(2, 3)); !got.Equal(New(1, 3)) {
		t.Errorf("Min = %v, want [1,3]", got)
	}
	if got := New(1, 4).Max(New(2, 3)); !got.Equal(New(2, 4)) {
		t.Errorf("Max = %v, want [2,4]", got)
	}
	encloses(t, "AddFloat", Point(1).AddFloat(0.5), 1.5, 1e-9)
	encloses(t, "MulFloat", Point(2).MulFloat(3), 6, 1e-9)
	encloses(t, "DivFloat", Point(6).DivFloat(2), 3, 1e-9)
}

func TestFreeFunctionForms(t *testing.T) {
	a, b := New(1, 2), New(3, 4)
	if !Add(a, b).Equal(a.Add(b)) || !Sub(a, b).Equal(a.Sub(b)) ||
		!Mul(a, b).Equal(a.Mul(b)) || !Div(a, b).Equal(a.Div(b)) ||
		!Hull(a, b).Equal(a.Hull(b)) || !Intersect(a, b).Equal(a.Intersect(b)) {
		t.Error("free function forms disagree with methods")
	}
}
