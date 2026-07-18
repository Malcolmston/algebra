package interval

import (
	"math"
	"testing"
)

func TestSqrtExpLog(t *testing.T) {
	encloses(t, "Sqrt(4)", Point(4).Sqrt(), 2, 1e-12)
	encloses(t, "Sqrt(2)", Point(2).Sqrt(), math.Sqrt2, 1e-12)
	// sqrt over an interval.
	s := New(4, 9).Sqrt()
	if !s.Contains(2) || !s.Contains(3) {
		t.Errorf("Sqrt([4,9]) = %v, want enclose [2,3]", s)
	}
	// negative part dropped.
	if !New(-5, -1).Sqrt().IsEmpty() {
		t.Error("Sqrt of negative interval should be empty")
	}
	encloses(t, "Exp(0)", Point(0).Exp(), 1, 1e-12)
	encloses(t, "Exp(1)", Point(1).Exp(), math.E, 1e-12)
	encloses(t, "Log(e)", Point(math.E).Log(), 1, 1e-12)
	encloses(t, "Log(1)", Point(1).Log(), 0, 1e-12)
	// exp(log(x)) encloses x.
	x := New(2, 5)
	round := x.Log().Exp()
	if !round.Contains(2) || !round.Contains(5) {
		t.Errorf("exp(log([2,5])) = %v, want enclose [2,5]", round)
	}
}

func TestHyperbolicAndAtan(t *testing.T) {
	encloses(t, "Sinh(1)", Point(1).Sinh(), math.Sinh(1), 1e-12)
	encloses(t, "Cosh(1)", Point(1).Cosh(), math.Cosh(1), 1e-12)
	encloses(t, "Tanh(1)", Point(1).Tanh(), math.Tanh(1), 1e-12)
	encloses(t, "Atan(1)", Point(1).Atan(), math.Pi/4, 1e-12)
	// cosh minimum at 0.
	c := New(-1, 2).Cosh()
	if !c.Contains(1) || c.Lo > 1+1e-12 {
		t.Errorf("Cosh([-1,2]) = %v, want lower bound near 1", c)
	}
}

func TestSinCos(t *testing.T) {
	tests := []struct {
		name string
		in   Interval
		fn   func(Interval) Interval
		ref  func(float64) float64
		pts  []float64
	}{
		{"sin(0)", Point(0), Interval.Sin, math.Sin, []float64{0}},
		{"sin(pi/2)", Point(math.Pi / 2), Interval.Sin, math.Sin, []float64{math.Pi / 2}},
		{"cos(0)", Point(0), Interval.Cos, math.Cos, []float64{0}},
		{"cos(pi)", Point(math.Pi), Interval.Cos, math.Cos, []float64{math.Pi}},
		{"sin range", New(0, 1), Interval.Sin, math.Sin, []float64{0, 0.5, 1}},
		{"cos range", New(0, 1), Interval.Cos, math.Cos, []float64{0, 0.5, 1}},
	}
	for _, tc := range tests {
		got := tc.fn(tc.in)
		for _, p := range tc.pts {
			if !got.Contains(tc.ref(p)) {
				t.Errorf("%s: %v does not enclose f(%g)=%g", tc.name, got, p, tc.ref(p))
			}
		}
	}
	// sin over interval containing pi/2 must reach up to 1.
	s := New(1, 2).Sin() // pi/2 ~= 1.5708 inside
	if s.Hi < 1-1e-12 {
		t.Errorf("Sin([1,2]) = %v, expected upper bound 1", s)
	}
	// cos over interval containing pi reaches -1.
	c := New(3, 3.3).Cos()
	if c.Lo > -1+1e-9 {
		t.Errorf("Cos([3,3.3]) = %v, expected lower bound near -1", c)
	}
	// wide interval => full range.
	full := New(0, 10).Sin()
	if full.Lo > -1 || full.Hi < 1 {
		t.Errorf("Sin([0,10]) = %v, want [-1,1]", full)
	}
}

func TestTan(t *testing.T) {
	encloses(t, "Tan(0)", Point(0).Tan(), 0, 1e-9)
	encloses(t, "Tan(pi/4)", Point(math.Pi/4).Tan(), 1, 1e-9)
	// interval straddling pi/2 => cos encloses 0 => entire.
	if !New(1, 2).Tan().IsEntire() {
		t.Errorf("Tan([1,2]) should be entire (crosses pi/2)")
	}
}

func TestIntPow(t *testing.T) {
	if !Point(2).IntPow(0).Equal(Point(1)) {
		t.Error("x^0 should be [1,1]")
	}
	encloses(t, "2^10", Point(2).IntPow(10), 1024, 1e-6)
	// even power of straddling interval: [-3,2]^2 = [0,9].
	p := New(-3, 2).IntPow(2)
	if !p.Contains(0) || !p.Contains(9) || p.Lo < 0 {
		t.Errorf("[-3,2]^2 = %v, want enclose [0,9] with Lo>=0", p)
	}
	// odd power monotone: [-2,3]^3 = [-8,27].
	q := New(-2, 3).IntPow(3)
	if !q.Contains(-8) || !q.Contains(27) {
		t.Errorf("[-2,3]^3 = %v, want enclose [-8,27]", q)
	}
	// negative exponent: 2^-1 = 0.5.
	encloses(t, "2^-1", Point(2).IntPow(-1), 0.5, 1e-12)
}

func TestPowReal(t *testing.T) {
	encloses(t, "2^0.5", Point(2).Pow(0.5), math.Sqrt2, 1e-9)
	encloses(t, "4^1.5", Point(4).Pow(1.5), 8, 1e-6)
	// non-positive base is out of domain.
	if !New(-1, 2).Pow(0.5).IsEmpty() {
		t.Error("Pow of non-positive base should be empty")
	}
}
