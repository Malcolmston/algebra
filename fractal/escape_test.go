package fractal

import (
	"math"
	"testing"
)

func approx(t *testing.T, got, want, tol float64, name string) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s: got %v, want %v (tol %v)", name, got, want, tol)
	}
}

func TestOrbit(t *testing.T) {
	// z0=0, c=1: 0 -> 1 -> 2 -> 5.
	got := Orbit(0, complex(1, 0), 3)
	want := []complex128{0, 1, 2, 5}
	if len(got) != len(want) {
		t.Fatalf("length: got %d want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("orbit[%d]: got %v want %v", i, got[i], want[i])
		}
	}
	if Orbit(0, 0, -1) != nil {
		t.Errorf("negative n should return nil")
	}
}

func TestMandelbrotEscape(t *testing.T) {
	// c=2: 0 -> 2 -> 6, escapes at n=2 with bailout 2.
	r := MandelbrotEscape(complex(2, 0), 100, 2)
	if !r.Escaped || r.Iterations != 2 {
		t.Errorf("c=2: got escaped=%v iters=%d, want true,2", r.Escaped, r.Iterations)
	}
	// c=0: never escapes.
	r = MandelbrotEscape(0, 50, 2)
	if r.Escaped || r.Iterations != 50 {
		t.Errorf("c=0: got escaped=%v iters=%d, want false,50", r.Escaped, r.Iterations)
	}
}

func TestInMandelbrotSet(t *testing.T) {
	cases := []struct {
		c    complex128
		want bool
	}{
		{0, true},
		{complex(-1, 0), true},
		{complex(-0.5, 0.5), true},
		{complex(1, 0), false},
		{complex(2, 0), false},
		{complex(0.4, 0.4), false},
	}
	for _, c := range cases {
		if got := InMandelbrotSet(c.c, 200); got != c.want {
			t.Errorf("InMandelbrotSet(%v): got %v want %v", c.c, got, c.want)
		}
	}
}

func TestInteriorTests(t *testing.T) {
	// Main cardioid.
	if !InMainCardioid(0) {
		t.Errorf("0 should be in main cardioid")
	}
	if !InMainCardioid(complex(0.25, 0)) { // cusp, on boundary
		t.Errorf("0.25 (cusp) should be in closed cardioid")
	}
	if InMainCardioid(complex(-1, 0)) {
		t.Errorf("-1 should not be in main cardioid")
	}
	// Period-2 bulb: disk radius 1/4 at -1.
	if !InPeriod2Bulb(complex(-1, 0)) {
		t.Errorf("-1 should be in period-2 bulb")
	}
	if !InPeriod2Bulb(complex(-0.75, 0)) { // boundary point c=-3/4
		t.Errorf("-0.75 should be on period-2 bulb boundary")
	}
	if InPeriod2Bulb(0) {
		t.Errorf("0 should not be in period-2 bulb")
	}
	// Points in either closed-form region are in the set.
	for _, c := range []complex128{0, complex(-0.2, 0.3), complex(-1, 0), complex(-0.9, 0.1)} {
		if (InMainCardioid(c) || InPeriod2Bulb(c)) && !InMandelbrotSet(c, 500) {
			t.Errorf("closed-form interior point %v not classified in set", c)
		}
	}
}

func TestJuliaEscape(t *testing.T) {
	// c=0: Julia set is the unit circle; |z|<1 stays bounded, |z|>1 escapes.
	if !InJuliaSet(complex(0.5, 0), 0, 200) {
		t.Errorf("|z|=0.5 should be in Julia set of c=0")
	}
	if InJuliaSet(complex(1.5, 0), 0, 200) {
		t.Errorf("|z|=1.5 should escape Julia set of c=0")
	}
}

func TestSmooth(t *testing.T) {
	// Non-escaped result returns the integer iteration count.
	r := EscapeResult{Escaped: false, Iterations: 100}
	if got := r.Smooth(256); got != 100 {
		t.Errorf("non-escaped Smooth: got %v want 100", got)
	}
	// Escaped smooth value stays close to the integer count and is finite.
	er := MandelbrotEscape(complex(0.4, 0.3), 1000, 256)
	if !er.Escaped {
		t.Fatalf("expected escape for c=0.4+0.3i")
	}
	s := er.Smooth(256)
	if math.IsNaN(s) || math.IsInf(s, 0) {
		t.Fatalf("smooth value not finite: %v", s)
	}
	if s < float64(er.Iterations)-1 || s > float64(er.Iterations)+1 {
		t.Errorf("smooth %v not within +-1 of iters %d", s, er.Iterations)
	}
}

func TestViewportPixelToComplex(t *testing.T) {
	v := Viewport{-1, 1, -1, 1}
	// Top-left pixel.
	if got := v.PixelToComplex(0, 0, 2, 2); got != complex(-1, 1) {
		t.Errorf("top-left: got %v want (-1+1i)", got)
	}
	// Bottom-right pixel.
	if got := v.PixelToComplex(1, 1, 2, 2); got != complex(1, -1) {
		t.Errorf("bottom-right: got %v want (1-1i)", got)
	}
	if v.SpanX() != 2 || v.SpanY() != 2 || v.Center() != 0 {
		t.Errorf("span/center wrong: %v %v %v", v.SpanX(), v.SpanY(), v.Center())
	}
	z := NewViewport(0, 0, 2).Zoom(0.5)
	approx(t, z.SpanX(), 2, 1e-12, "zoom spanX")
}

func TestMandelbrotGrid(t *testing.T) {
	v := NewViewport(-0.5, 0, 1.5)
	g := MandelbrotGrid(v, 40, 40, 100, 256)
	if g.Width != 40 || g.Height != 40 {
		t.Fatalf("grid dims wrong")
	}
	// The center region contains set members (value == maxIter).
	if g.Max() < 100 {
		t.Errorf("expected some interior cells at maxIter=100, max=%v", g.Max())
	}
	// Corner (2, 2i approx) escapes quickly -> small value.
	if g.At(39, 0) > 20 {
		t.Errorf("corner should escape quickly, got %v", g.At(39, 0))
	}
	// JuliaGrid smoke test.
	jg := JuliaGrid(complex(-0.8, 0.156), v, 20, 20, 100, 256)
	if jg.CountFinite() != 400 {
		t.Errorf("julia grid has non-finite cells: %d finite", jg.CountFinite())
	}
}

func BenchmarkMandelbrotGrid(b *testing.B) {
	v := NewViewport(-0.75, 0, 1.5)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		MandelbrotGrid(v, 128, 128, 256, 256)
	}
}
