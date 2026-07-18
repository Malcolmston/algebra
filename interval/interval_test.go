package interval

import (
	"math"
	"testing"
)

// encloses is a test helper asserting that the interval a contains the true
// value v and that its width does not exceed tol.
func encloses(t *testing.T, name string, a Interval, v, tol float64) {
	t.Helper()
	if a.IsEmpty() {
		t.Fatalf("%s: got empty interval, want enclosure of %v", name, v)
	}
	if !a.Contains(v) {
		t.Errorf("%s: %v does not contain true value %v", name, a, v)
	}
	if w := a.Width(); w > tol {
		t.Errorf("%s: width %g exceeds tolerance %g (interval %v)", name, w, tol, a)
	}
}

func TestConstructorsAndPredicates(t *testing.T) {
	if !Empty().IsEmpty() {
		t.Error("Empty should be empty")
	}
	if Empty().Contains(0) {
		t.Error("Empty contains nothing")
	}
	if !Entire().IsEntire() {
		t.Error("Entire should be entire")
	}
	if !Point(3).IsPoint() {
		t.Error("Point(3) should be a point")
	}
	// New swaps reversed bounds.
	got := New(5, 1)
	if got.Lo != 1 || got.Hi != 5 {
		t.Errorf("New(5,1) = %v, want [1,5]", got)
	}
	if !New(math.NaN(), 1).IsEmpty() {
		t.Error("New with NaN should be empty")
	}
	if !Point(2).IsBounded() || Entire().IsBounded() {
		t.Error("boundedness classification wrong")
	}
}

func TestMeasures(t *testing.T) {
	a := New(2, 6)
	if w := a.Width(); w < 4 || w > 4+1e-12 {
		t.Errorf("Width = %g, want ~4", w)
	}
	if r := a.Radius(); r < 2 || r > 2+1e-12 {
		t.Errorf("Radius = %g, want ~2", r)
	}
	if m := a.Midpoint(); math.Abs(m-4) > 1e-12 {
		t.Errorf("Midpoint = %g, want 4", m)
	}
	if math.Abs(a.Mag()-6) > 0 {
		t.Errorf("Mag = %g, want 6", a.Mag())
	}
	if math.Abs(a.Mig()-2) > 0 {
		t.Errorf("Mig = %g, want 2", a.Mig())
	}
	b := New(-3, 5)
	if b.Mig() != 0 {
		t.Errorf("Mig of straddling interval = %g, want 0", b.Mig())
	}
	if b.Mag() != 5 {
		t.Errorf("Mag = %g, want 5", b.Mag())
	}
}

func TestContainsAndSetOps(t *testing.T) {
	a := New(1, 4)
	b := New(3, 6)
	if !a.Contains(2) || a.Contains(5) {
		t.Error("Contains wrong")
	}
	if !a.ContainsInterval(New(2, 3)) || a.ContainsInterval(b) {
		t.Error("ContainsInterval wrong")
	}
	if !a.Overlaps(b) || a.Overlaps(New(10, 11)) {
		t.Error("Overlaps wrong")
	}
	if got := a.Intersect(b); !got.Equal(New(3, 4)) {
		t.Errorf("Intersect = %v, want [3,4]", got)
	}
	if got := a.Intersect(New(10, 11)); !got.IsEmpty() {
		t.Errorf("disjoint Intersect should be empty, got %v", got)
	}
	if got := a.Hull(b); !got.Equal(New(1, 6)) {
		t.Errorf("Hull = %v, want [1,6]", got)
	}
	if got := a.Hull(Empty()); !got.Equal(a) {
		t.Errorf("Hull with empty = %v, want %v", got, a)
	}
}

func TestString(t *testing.T) {
	if s := New(1, 2).String(); s != "[1, 2]" {
		t.Errorf("String = %q, want [1, 2]", s)
	}
	if s := Empty().String(); s != "[empty]" {
		t.Errorf("empty String = %q", s)
	}
}
