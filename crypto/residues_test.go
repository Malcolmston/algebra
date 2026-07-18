package crypto

import (
	"math/big"
	"testing"
)

func TestJacobi(t *testing.T) {
	// Reference values from standard tables.
	cases := []struct {
		a, n, want int64
	}{
		{1, 3, 1},
		{2, 3, -1},
		{2, 7, 1}, // 3^2=9≡2 mod7
		{3, 7, -1},
		{1001, 9907, -1},
		{5, 21, 1}, // (5/3)(5/7)=(-1)(-1)=1
		{6, 7, -1},
		{0, 3, 0},
	}
	for _, c := range cases {
		if got := Jacobi(bi(c.a), bi(c.n)); int64(got) != c.want {
			t.Errorf("Jacobi(%d,%d)=%d want %d", c.a, c.n, got, c.want)
		}
	}
}

func TestLegendre(t *testing.T) {
	// Squares mod 11 are {1,3,4,5,9}.
	residues := map[int64]bool{1: true, 3: true, 4: true, 5: true, 9: true}
	for a := int64(1); a < 11; a++ {
		want := -1
		if residues[a] {
			want = 1
		}
		if got := Legendre(bi(a), bi(11)); got != want {
			t.Errorf("Legendre(%d,11)=%d want %d", a, got, want)
		}
	}
	if Legendre(bi(11), bi(11)) != 0 {
		t.Error("Legendre(11,11) want 0")
	}
	// Jacobi must agree with Legendre for prime modulus.
	for a := int64(1); a < 11; a++ {
		if Jacobi(bi(a), bi(11)) != Legendre(bi(a), bi(11)) {
			t.Errorf("Jacobi/Legendre disagree at a=%d", a)
		}
	}
}

func TestModSqrt(t *testing.T) {
	primes := []int64{7, 11, 13, 17, 23, 31, 101} // mix of p≡3 and p≡1 mod4
	for _, pp := range primes {
		p := bi(pp)
		for a := int64(0); a < pp; a++ {
			r, err := ModSqrt(bi(a), p)
			if Legendre(bi(a), p) == -1 {
				if err == nil {
					t.Errorf("ModSqrt(%d,%d) expected non-residue error", a, pp)
				}
				continue
			}
			if err != nil {
				t.Fatalf("ModSqrt(%d,%d) err %v", a, pp, err)
			}
			chk := new(big.Int).Mul(r, r)
			chk.Mod(chk, p)
			if chk.Int64() != a%pp {
				t.Errorf("ModSqrt(%d,%d)=%d, r^2=%d", a, pp, r.Int64(), chk.Int64())
			}
		}
	}
}

func TestMultiplicativeOrder(t *testing.T) {
	cases := []struct{ a, n, want int64 }{
		{2, 7, 3}, // 2,4,1
		{3, 7, 6}, // primitive root
		{2, 9, 6},
		{4, 7, 3},
		{10, 3, 1},
	}
	for _, c := range cases {
		got, err := MultiplicativeOrder(bi(c.a), bi(c.n))
		if err != nil {
			t.Fatalf("MultiplicativeOrder(%d,%d) err %v", c.a, c.n, err)
		}
		if got.Int64() != c.want {
			t.Errorf("MultiplicativeOrder(%d,%d)=%d want %d", c.a, c.n, got.Int64(), c.want)
		}
	}
}

func TestPrimitiveRoot(t *testing.T) {
	cases := []struct{ p, want int64 }{
		{7, 3},
		{11, 2},
		{13, 2},
		{23, 5},
	}
	for _, c := range cases {
		g := PrimitiveRoot(bi(c.p))
		if g.Int64() != c.want {
			t.Errorf("PrimitiveRoot(%d)=%d want %d", c.p, g.Int64(), c.want)
		}
		if !IsPrimitiveRoot(g, bi(c.p)) {
			t.Errorf("IsPrimitiveRoot(%d,%d) want true", g.Int64(), c.p)
		}
		// Its order must equal p-1.
		ord, _ := MultiplicativeOrder(g, bi(c.p))
		if ord.Int64() != c.p-1 {
			t.Errorf("order of primitive root %d mod %d = %d want %d", g.Int64(), c.p, ord.Int64(), c.p-1)
		}
	}
}
