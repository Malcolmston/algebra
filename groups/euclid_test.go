package groups

import "testing"

func TestGcd(t *testing.T) {
	cases := []struct {
		a, b, want int
	}{
		{12, 8, 4},
		{54, 24, 6},
		{17, 5, 1},
		{0, 0, 0},
		{0, 9, 9},
		{-12, 8, 4},
		{-12, -8, 4},
		{1071, 462, 21},
	}
	for _, c := range cases {
		if got := Gcd(c.a, c.b); got != c.want {
			t.Errorf("Gcd(%d,%d)=%d want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestExtendedGcd(t *testing.T) {
	cases := [][2]int{{240, 46}, {17, 5}, {-12, 8}, {1071, 462}, {0, 7}}
	for _, c := range cases {
		g, x, y := ExtendedGcd(c[0], c[1])
		if g != Gcd(c[0], c[1]) {
			t.Errorf("ExtendedGcd(%d,%d) g=%d want %d", c[0], c[1], g, Gcd(c[0], c[1]))
		}
		if c[0]*x+c[1]*y != g {
			t.Errorf("Bezout failed: %d*%d + %d*%d != %d", c[0], x, c[1], y, g)
		}
	}
}

func TestLcm(t *testing.T) {
	cases := []struct{ a, b, want int }{
		{4, 6, 12},
		{21, 6, 42},
		{0, 5, 0},
		{7, 1, 7},
		{-4, 6, 12},
	}
	for _, c := range cases {
		if got := Lcm(c.a, c.b); got != c.want {
			t.Errorf("Lcm(%d,%d)=%d want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestGcdLcmMany(t *testing.T) {
	if g := GcdMany(24, 36, 48); g != 12 {
		t.Errorf("GcdMany=%d want 12", g)
	}
	if l := LcmMany(4, 6, 10); l != 60 {
		t.Errorf("LcmMany=%d want 60", l)
	}
	if LcmMany() != 1 {
		t.Errorf("empty LcmMany want 1")
	}
	if GcdMany() != 0 {
		t.Errorf("empty GcdMany want 0")
	}
}

func TestCoprime(t *testing.T) {
	if !Coprime(9, 16) {
		t.Errorf("9,16 should be coprime")
	}
	if Coprime(9, 15) {
		t.Errorf("9,15 not coprime")
	}
}
