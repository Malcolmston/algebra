package seq

import "testing"

func TestPolygonalFamilies(t *testing.T) {
	tri := []uint64{0, 1, 3, 6, 10, 15, 21, 28, 36, 45}
	sq := []uint64{0, 1, 4, 9, 16, 25, 36, 49, 64, 81}
	pent := []uint64{0, 1, 5, 12, 22, 35, 51, 70, 92, 117}
	hex := []uint64{0, 1, 6, 15, 28, 45, 66, 91, 120, 153}
	hept := []uint64{0, 1, 7, 18, 34, 55, 81, 112, 148, 189}
	oct := []uint64{0, 1, 8, 21, 40, 65, 96, 133, 176, 225}
	for n := 0; n < 10; n++ {
		if Triangular(n) != tri[n] {
			t.Errorf("Triangular(%d) = %d, want %d", n, Triangular(n), tri[n])
		}
		if Square(n) != sq[n] {
			t.Errorf("Square(%d) = %d, want %d", n, Square(n), sq[n])
		}
		if Pentagonal(n) != pent[n] {
			t.Errorf("Pentagonal(%d) = %d, want %d", n, Pentagonal(n), pent[n])
		}
		if Hexagonal(n) != hex[n] {
			t.Errorf("Hexagonal(%d) = %d, want %d", n, Hexagonal(n), hex[n])
		}
		if Heptagonal(n) != hept[n] {
			t.Errorf("Heptagonal(%d) = %d, want %d", n, Heptagonal(n), hept[n])
		}
		if Octagonal(n) != oct[n] {
			t.Errorf("Octagonal(%d) = %d, want %d", n, Octagonal(n), oct[n])
		}
		// Polygonal must agree with the specialised formulas.
		if Polygonal(3, n) != tri[n] || Polygonal(4, n) != sq[n] ||
			Polygonal(5, n) != pent[n] || Polygonal(6, n) != hex[n] ||
			Polygonal(7, n) != hept[n] || Polygonal(8, n) != oct[n] {
			t.Errorf("Polygonal disagrees with specialised formula at n=%d", n)
		}
	}
}

func TestGeneralizedPentagonal(t *testing.T) {
	// Signed index sequence 0,1,-1,2,-2,3,-3 -> 0,1,2,5,7,12,15.
	idx := []int{0, 1, -1, 2, -2, 3, -3}
	want := []int64{0, 1, 2, 5, 7, 12, 15}
	for i, k := range idx {
		if got := GeneralizedPentagonal(k); got != want[i] {
			t.Errorf("GeneralizedPentagonal(%d) = %d, want %d", k, got, want[i])
		}
	}
}

func TestPredicates(t *testing.T) {
	triSet := map[uint64]bool{}
	pentSet := map[uint64]bool{}
	hexSet := map[uint64]bool{}
	pronSet := map[uint64]bool{}
	for n := 0; n < 200; n++ {
		triSet[Triangular(n)] = true
		pentSet[Pentagonal(n)] = true
		hexSet[Hexagonal(n)] = true
		pronSet[Pronic(n)] = true
	}
	for x := uint64(0); x <= 500; x++ {
		if IsTriangular(x) != triSet[x] {
			t.Errorf("IsTriangular(%d) wrong", x)
		}
		if IsSquare(x) != IsPerfectSquare(x) {
			t.Errorf("IsSquare(%d) wrong", x)
		}
		if IsPentagonal(x) != pentSet[x] {
			t.Errorf("IsPentagonal(%d) wrong", x)
		}
		if IsHexagonal(x) != hexSet[x] {
			t.Errorf("IsHexagonal(%d) wrong", x)
		}
		if IsPronic(x) != pronSet[x] {
			t.Errorf("IsPronic(%d) wrong", x)
		}
	}
}

func TestCentered(t *testing.T) {
	if got := []uint64{CenteredTriangular(0), CenteredTriangular(1), CenteredTriangular(2), CenteredTriangular(3)}; got[0] != 1 || got[1] != 4 || got[2] != 10 || got[3] != 19 {
		t.Errorf("CenteredTriangular mismatch: %v", got)
	}
	if got := []uint64{CenteredSquare(0), CenteredSquare(1), CenteredSquare(2), CenteredSquare(3)}; got[0] != 1 || got[1] != 5 || got[2] != 13 || got[3] != 25 {
		t.Errorf("CenteredSquare mismatch: %v", got)
	}
	if got := []uint64{CenteredHexagonal(0), CenteredHexagonal(1), CenteredHexagonal(2), CenteredHexagonal(3)}; got[0] != 1 || got[1] != 7 || got[2] != 19 || got[3] != 37 {
		t.Errorf("CenteredHexagonal mismatch: %v", got)
	}
}

func TestSpatialFigurate(t *testing.T) {
	tet := []uint64{0, 1, 4, 10, 20, 35, 56, 84, 120}
	pyr := []uint64{0, 1, 5, 14, 30, 55, 91, 140, 204}
	ptp := []uint64{0, 1, 5, 15, 35, 70, 126, 210, 330}
	pron := []uint64{0, 2, 6, 12, 20, 30, 42, 56, 72}
	cube := []uint64{0, 1, 8, 27, 64, 125, 216, 343, 512}
	for n := 0; n < 9; n++ {
		if Tetrahedral(n) != tet[n] {
			t.Errorf("Tetrahedral(%d) = %d, want %d", n, Tetrahedral(n), tet[n])
		}
		if SquarePyramidal(n) != pyr[n] {
			t.Errorf("SquarePyramidal(%d) = %d, want %d", n, SquarePyramidal(n), pyr[n])
		}
		if Pentatope(n) != ptp[n] {
			t.Errorf("Pentatope(%d) = %d, want %d", n, Pentatope(n), ptp[n])
		}
		if Pronic(n) != pron[n] {
			t.Errorf("Pronic(%d) = %d, want %d", n, Pronic(n), pron[n])
		}
		if Cube(n) != cube[n] {
			t.Errorf("Cube(%d) = %d, want %d", n, Cube(n), cube[n])
		}
	}
}

func TestStarGnomonic(t *testing.T) {
	star := []uint64{1, 13, 37, 73, 121, 181}
	for i, w := range star {
		if got := StarNumber(i + 1); got != w {
			t.Errorf("StarNumber(%d) = %d, want %d", i+1, got, w)
		}
	}
	for n := 1; n <= 10; n++ {
		if Gnomonic(n) != uint64(2*n-1) {
			t.Errorf("Gnomonic(%d) wrong", n)
		}
	}
}

func TestIntSqrt(t *testing.T) {
	for v := uint64(0); v <= 1000; v++ {
		r := IntSqrt(v)
		if r*r > v || (r+1)*(r+1) <= v {
			t.Errorf("IntSqrt(%d) = %d incorrect", v, r)
		}
	}
	// A large exact square near the top of the range.
	big := uint64(4000000000)
	if IntSqrt(big*big) != big {
		t.Errorf("IntSqrt of large square wrong")
	}
}
