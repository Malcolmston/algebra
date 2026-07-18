package interval

import (
	"math"
	"testing"
)

func TestNewtonSqrt2(t *testing.T) {
	// f(x) = x^2 - 2, root sqrt(2). f'(x) = 2x.
	f := func(x Interval) Interval { return x.Square().SubFloat(2) }
	df := func(x Interval) Interval { return x.MulFloat(2) }
	got, verified := Newton(f, df, New(1, 2), 1e-12, 100)
	if got.IsEmpty() {
		t.Fatal("Newton returned empty for a bracketed root")
	}
	if !got.Contains(math.Sqrt2) {
		t.Errorf("Newton = %v, does not enclose sqrt(2)=%g", got, math.Sqrt2)
	}
	if !verified {
		t.Error("Newton should verify existence and uniqueness on [1,2]")
	}
	if got.Width() > 1e-10 {
		t.Errorf("Newton width %g too wide", got.Width())
	}
}

func TestNewtonCosRoot(t *testing.T) {
	// f(x) = cos(x), root pi/2 in [1,2]. f'(x) = -sin(x).
	f := func(x Interval) Interval { return x.Cos() }
	df := func(x Interval) Interval { return x.Sin().Neg() }
	got, verified := Newton(f, df, New(1, 2), 1e-10, 100)
	if !got.Contains(math.Pi / 2) {
		t.Errorf("Newton = %v, does not enclose pi/2", got)
	}
	if !verified {
		t.Error("expected verified root for cos on [1,2]")
	}
}

func TestNewtonNoRoot(t *testing.T) {
	// f(x) = x^2 + 1 has no real root.
	f := func(x Interval) Interval { return x.Square().AddFloat(1) }
	df := func(x Interval) Interval { return x.MulFloat(2) }
	got, verified := Newton(f, df, New(1, 2), 1e-12, 100)
	if !got.IsEmpty() {
		t.Errorf("Newton should prove no root, got %v", got)
	}
	if verified {
		t.Error("no-root case should not report verified")
	}
}

func TestBisect(t *testing.T) {
	l, r := Bisect(New(0, 4))
	if !l.Equal(New(0, 2)) || !r.Equal(New(2, 4)) {
		t.Errorf("Bisect([0,4]) = %v,%v want [0,2],[2,4]", l, r)
	}
}
