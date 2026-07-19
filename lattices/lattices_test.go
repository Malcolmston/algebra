package lattices

import (
	"fmt"
	"math"
	"math/big"
	"testing"
)

const tol = 1e-9

func almost(a, b float64) bool { return math.Abs(a-b) <= 1e-9 }

// ExampleBasis_LLL demonstrates LLL reduction shrinking a skewed basis while
// preserving the lattice determinant.
func ExampleBasis_LLL() {
	b := NewBasis(
		NewVec(1, 1, 1),
		NewVec(-1, 0, 2),
		NewVec(3, 5, 6),
	)
	red, _ := b.LLL(0.75)
	fmt.Println(red.IsLLLReduced(0.75, 1e-9))
	fmt.Printf("%.0f\n", red.Determinant())
	// Output:
	// true
	// 3
}

// ExampleBasis_ShortestVector finds the shortest nonzero vector of the D2
// checkerboard lattice, whose length is sqrt(2).
func ExampleBasis_ShortestVector() {
	b := NewBasis(NewVec(1, 1), NewVec(2, 0))
	res, _ := b.ShortestVector()
	fmt.Printf("%.6f\n", res.Norm())
	// Output:
	// 1.414214
}

func TestVecOps(t *testing.T) {
	v := NewVec(3, 4)
	w := NewVec(1, 2)
	if !almost(v.Norm(), 5) {
		t.Errorf("norm = %v want 5", v.Norm())
	}
	if !almost(v.Norm2(), 25) {
		t.Errorf("norm2 = %v", v.Norm2())
	}
	if !almost(v.Dot(w), 11) {
		t.Errorf("dot = %v want 11", v.Dot(w))
	}
	if !v.Add(w).Equal(NewVec(4, 6)) {
		t.Errorf("add = %v", v.Add(w))
	}
	if !v.Sub(w).Equal(NewVec(2, 2)) {
		t.Errorf("sub = %v", v.Sub(w))
	}
	if !v.Scale(2).Equal(NewVec(6, 8)) {
		t.Errorf("scale = %v", v.Scale(2))
	}
	if !almost(v.L1(), 7) {
		t.Errorf("l1 = %v", v.L1())
	}
	if !almost(v.LInf(), 4) {
		t.Errorf("linf = %v", v.LInf())
	}
	if !almost(NewVec(1, 0).Angle(NewVec(0, 1)), math.Pi/2) {
		t.Errorf("angle wrong")
	}
	if !almost(NewVec(3, 0).Dist(NewVec(0, 4)), 5) {
		t.Errorf("dist wrong")
	}
}

func TestCombine(t *testing.T) {
	got := Combine([]float64{2, -1}, []Vec{NewVec(1, 0), NewVec(0, 1)})
	if !got.Equal(NewVec(2, -1)) {
		t.Errorf("combine = %v", got)
	}
}

func TestGramAndDeterminant(t *testing.T) {
	tests := []struct {
		name string
		b    Basis
		det  float64
	}{
		{"identity2", IdentityBasis(2), 1},
		{"identity3", IdentityBasis(3), 1},
		{"scaled", NewBasis(NewVec(2, 0), NewVec(0, 3)), 6},
		{"checker", NewBasis(NewVec(1, 1), NewVec(2, 0)), 2},
		{"skew3", NewBasis(NewVec(1, 1, 1), NewVec(-1, 0, 2), NewVec(3, 5, 6)), 3},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if d := tt.b.Determinant(); !almost(d, tt.det) {
				t.Errorf("det = %v want %v", d, tt.det)
			}
			g := tt.b.Gram()
			if !g.IsSymmetric(tol) {
				t.Errorf("gram not symmetric")
			}
		})
	}
}

func TestGramDeterminantExact(t *testing.T) {
	b := NewBasis(NewVec(2, 0), NewVec(0, 3))
	d := b.GramDeterminantRat()
	if d.Cmp(big.NewRat(36, 1)) != 0 {
		t.Errorf("exact gram det = %v want 36", d)
	}
}

func TestGramSchmidt(t *testing.T) {
	b := NewBasis(NewVec(1, 1, 0), NewVec(1, 0, 1), NewVec(0, 1, 1))
	gs := b.Orthogonalize()
	// b*_0 == b_0
	if !gs.Star[0].ApproxEqual(NewVec(1, 1, 0), tol) {
		t.Errorf("star0 = %v", gs.Star[0])
	}
	// orthogonality of the star vectors
	for i := 0; i < 3; i++ {
		for j := i + 1; j < 3; j++ {
			if !almost(gs.Star[i].Dot(gs.Star[j]), 0) {
				t.Errorf("star %d,%d not orthogonal: %v", i, j, gs.Star[i].Dot(gs.Star[j]))
			}
		}
	}
	// product of star norms == covolume
	if !almost(b.ProductOfStarNorms(), b.Determinant()) {
		t.Errorf("prod star norms %v != det %v", b.ProductOfStarNorms(), b.Determinant())
	}
}

func TestLLLReducesAndPreservesLattice(t *testing.T) {
	tests := []Basis{
		NewBasis(NewVec(1, 1, 1), NewVec(-1, 0, 2), NewVec(3, 5, 6)),
		NewBasis(NewVec(201, 37), NewVec(1537, 279)),
		NewBasis(NewVec(1, 0, 0, 1), NewVec(0, 1, 0, 1), NewVec(0, 0, 1, 1), NewVec(1, 1, 1, 2)),
	}
	for i, b := range tests {
		t.Run(fmt.Sprintf("case%d", i), func(t *testing.T) {
			red, err := b.LLL(0.75)
			if err != nil {
				t.Fatalf("LLL error: %v", err)
			}
			if !red.IsLLLReduced(0.75, 1e-7) {
				t.Errorf("not LLL reduced")
			}
			if rel := math.Abs(red.Determinant()-b.Determinant()) / (1 + math.Abs(b.Determinant())); rel > 1e-9 {
				t.Errorf("determinant changed: %v vs %v", red.Determinant(), b.Determinant())
			}
			// reduction never increases the shortest basis vector length below lambda_1
			same, err := SameLattice(b, red)
			if err != nil {
				t.Fatalf("SameLattice error: %v", err)
			}
			if !same {
				t.Errorf("LLL changed the lattice")
			}
		})
	}
}

func TestLLLBadParameter(t *testing.T) {
	b := IdentityBasis(2)
	if _, err := b.LLL(2); err != ErrBadParameter {
		t.Errorf("expected ErrBadParameter, got %v", err)
	}
}

func TestSizeReduced(t *testing.T) {
	b := NewBasis(NewVec(1, 0), NewVec(10, 1))
	sr := b.SizeReduced()
	if !sr.IsSizeReduced(1e-9) {
		t.Errorf("not size reduced: %v", sr)
	}
}

func TestShortestVector(t *testing.T) {
	tests := []struct {
		name   string
		b      Basis
		lambda float64
	}{
		{"Z2", IdentityBasis(2), 1},
		{"scaled", NewBasis(NewVec(2, 0), NewVec(0, 3)), 2},
		{"checker", NewBasis(NewVec(1, 1), NewVec(2, 0)), math.Sqrt2},
		{"skew", NewBasis(NewVec(12, 2), NewVec(13, 4)), math.Sqrt(5)},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			res, err := tt.b.ShortestVector()
			if err != nil {
				t.Fatalf("error: %v", err)
			}
			if !almost(res.Norm(), tt.lambda) {
				t.Errorf("lambda_1 = %v want %v", res.Norm(), tt.lambda)
			}
			// the returned vector must actually be a lattice point of that norm
			if !almost(res.Vector.Norm(), tt.lambda) {
				t.Errorf("vector norm mismatch")
			}
		})
	}
}

func TestFirstMinimumBelowBounds(t *testing.T) {
	b := NewBasis(NewVec(1, 1), NewVec(2, 0))
	l1, err := b.FirstMinimum()
	if err != nil {
		t.Fatal(err)
	}
	if l1 > b.MinkowskiBoundBasis()+1e-9 {
		t.Errorf("lambda_1 %v exceeds Minkowski bound %v", l1, b.MinkowskiBoundBasis())
	}
	if l1 > b.HermiteBoundBasis()+1e-9 {
		t.Errorf("lambda_1 %v exceeds Hermite bound %v", l1, b.HermiteBoundBasis())
	}
}

func TestClosestVectorAndBabai(t *testing.T) {
	b := NewBasis(NewVec(2, 0), NewVec(0, 2))
	target := NewVec(1.1, 2.9)
	res, err := b.ClosestVector(target)
	if err != nil {
		t.Fatal(err)
	}
	if !res.Vector.ApproxEqual(NewVec(2, 2), tol) {
		t.Errorf("closest = %v want [2 2]", res.Vector)
	}
	// Babai methods return lattice points
	if br := b.BabaiRound(target); !br.ApproxEqual(NewVec(2, 2), tol) {
		t.Errorf("babai round = %v", br)
	}
	if bnp := b.BabaiNearestPlane(target); !bnp.ApproxEqual(NewVec(2, 2), tol) {
		t.Errorf("babai nearest plane = %v", bnp)
	}
}

func TestClosestVectorIsAtLeastAsGoodAsBabai(t *testing.T) {
	b := NewBasis(NewVec(2, 3), NewVec(3, 5))
	target := NewVec(4.4, 6.1)
	res, err := b.ClosestVector(target)
	if err != nil {
		t.Fatal(err)
	}
	babai := b.BabaiNearestPlane(target)
	if res.Norm2 > babai.Dist2(target)+1e-9 {
		t.Errorf("enumeration worse than Babai: %v vs %v", res.Norm2, babai.Dist2(target))
	}
}

func TestDual(t *testing.T) {
	b := NewBasis(NewVec(2, 0), NewVec(0, 4))
	dual, err := b.Dual()
	if err != nil {
		t.Fatal(err)
	}
	if !dual[0].ApproxEqual(NewVec(0.5, 0), tol) || !dual[1].ApproxEqual(NewVec(0, 0.25), tol) {
		t.Errorf("dual = %v", dual)
	}
	if !b.IsDualTo(dual, 1e-9) {
		t.Errorf("IsDualTo failed")
	}
	if !almost(b.DualDeterminant(), 1.0/b.Determinant()) {
		t.Errorf("dual determinant wrong")
	}
	bi, err := b.Biorthogonality()
	if err != nil {
		t.Fatal(err)
	}
	if !bi.ApproxEqual(IdentityMatrix(2), 1e-9) {
		t.Errorf("biorthogonality not identity: %v", bi)
	}
}

func TestHermiteNormalForm(t *testing.T) {
	tests := []struct {
		name string
		in   [][]int64
		want [][]int64
	}{
		{"2x2", [][]int64{{2, 3}, {4, 5}}, [][]int64{{2, 0}, {0, 1}}},
		{"identityish", [][]int64{{1, 0}, {0, 1}}, [][]int64{{1, 0}, {0, 1}}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			h := NewIntMatrix(tt.in).HermiteNormalForm()
			if !h.Equal(NewIntMatrix(tt.want)) {
				t.Errorf("HNF = \n%v\nwant\n%v", h, NewIntMatrix(tt.want))
			}
			if !h.IsHermiteNormalForm() {
				t.Errorf("output not recognized as HNF")
			}
		})
	}
}

func TestHNFTransform(t *testing.T) {
	m := NewIntMatrix([][]int64{{2, 3, 1}, {4, 1, 5}, {7, 2, 3}})
	h, u := m.HermiteNormalFormWithTransform()
	if !u.Mul(m).Equal(h) {
		t.Errorf("U*M != H")
	}
	// determinant of unimodular U is +-1
	d, _ := u.Det()
	if d.CmpAbs(big.NewInt(1)) != 0 {
		t.Errorf("U not unimodular, det = %v", d)
	}
}

func TestIntDeterminant(t *testing.T) {
	tests := []struct {
		in   [][]int64
		want int64
	}{
		{[][]int64{{2, 3}, {4, 5}}, -2},
		{[][]int64{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 1},
		{[][]int64{{6, 1, 1}, {4, -2, 5}, {2, 8, 7}}, -306},
	}
	for _, tt := range tests {
		d, err := NewIntMatrix(tt.in).Det()
		if err != nil {
			t.Fatal(err)
		}
		if d.Cmp(big.NewInt(tt.want)) != 0 {
			t.Errorf("det = %v want %v", d, tt.want)
		}
	}
}

func TestSameLattice(t *testing.T) {
	a := IdentityBasis(2)
	b := NewBasis(NewVec(1, 1), NewVec(0, 1)) // unimodular transform of Z^2
	same, err := SameLattice(a, b)
	if err != nil {
		t.Fatal(err)
	}
	if !same {
		t.Errorf("expected same lattice")
	}
	c := NewBasis(NewVec(2, 0), NewVec(0, 1)) // index-2 sublattice
	same, _ = SameLattice(a, c)
	if same {
		t.Errorf("expected different lattices")
	}
}

func TestSuccessiveMinima(t *testing.T) {
	tests := []struct {
		name string
		b    Basis
		want []float64
	}{
		{"Z2", IdentityBasis(2), []float64{1, 1}},
		{"scaled", NewBasis(NewVec(1, 0), NewVec(0, 3)), []float64{1, 3}},
		{"checker", NewBasis(NewVec(1, 1), NewVec(2, 0)), []float64{math.Sqrt2, math.Sqrt2}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			sm, err := tt.b.SuccessiveMinima()
			if err != nil {
				t.Fatal(err)
			}
			if len(sm) != len(tt.want) {
				t.Fatalf("len = %d want %d", len(sm), len(tt.want))
			}
			for i := range sm {
				if !almost(sm[i], tt.want[i]) {
					t.Errorf("lambda_%d = %v want %v", i+1, sm[i], tt.want[i])
				}
			}
		})
	}
}

func TestKissingNumber(t *testing.T) {
	tests := []struct {
		b    Basis
		want int
	}{
		{IdentityBasis(2), 4},
		{NewBasis(NewVec(1, 1), NewVec(2, 0)), 4}, // D2 ~ Z2 rotated
	}
	for i, tt := range tests {
		kn, err := tt.b.KissingNumber()
		if err != nil {
			t.Fatal(err)
		}
		if kn != tt.want {
			t.Errorf("case %d: kissing = %d want %d", i, kn, tt.want)
		}
	}
}

func TestBallVolume(t *testing.T) {
	if !almost(UnitBallVolume(2), math.Pi) {
		t.Errorf("V2 = %v want pi", UnitBallVolume(2))
	}
	if !almost(UnitBallVolume(3), 4.0/3*math.Pi) {
		t.Errorf("V3 = %v", UnitBallVolume(3))
	}
	if !almost(BallVolume(2, 2), 4*math.Pi) {
		t.Errorf("BallVolume(2,2) = %v", BallVolume(2, 2))
	}
}

func TestHermiteConstant(t *testing.T) {
	if !almost(HermiteConstant(1), 1) {
		t.Errorf("gamma_1 wrong")
	}
	if !almost(HermiteConstant(2), 2/math.Sqrt(3)) {
		t.Errorf("gamma_2 wrong")
	}
	if !almost(HermiteConstant(8), 2) {
		t.Errorf("gamma_8 wrong")
	}
}

func TestMatrixInverseAndSolve(t *testing.T) {
	m := NewMatrix([][]float64{{4, 7}, {2, 6}})
	inv, err := m.Inverse()
	if err != nil {
		t.Fatal(err)
	}
	prod := m.Mul(inv)
	if !prod.ApproxEqual(IdentityMatrix(2), 1e-9) {
		t.Errorf("M*Minv != I: %v", prod)
	}
	x, err := m.Solve(NewVec(1, 0))
	if err != nil {
		t.Fatal(err)
	}
	if !m.MulVec(x).ApproxEqual(NewVec(1, 0), 1e-9) {
		t.Errorf("solve wrong")
	}
	d, _ := m.Det()
	if !almost(d, 10) {
		t.Errorf("det = %v want 10", d)
	}
}

func TestRatMatrixExact(t *testing.T) {
	m := NewRatMatrix([][]int64{{1, 2}, {3, 4}})
	d, err := m.Det()
	if err != nil {
		t.Fatal(err)
	}
	if d.Cmp(big.NewRat(-2, 1)) != 0 {
		t.Errorf("exact det = %v want -2", d)
	}
	inv, err := m.Inverse()
	if err != nil {
		t.Fatal(err)
	}
	if !inv.Mul(m).Equal(IdentityRatMatrix(2)) {
		t.Errorf("rat inverse wrong")
	}
}

func TestOrthogonalityDefect(t *testing.T) {
	// orthogonal basis => defect 1
	b := NewBasis(NewVec(3, 0), NewVec(0, 5))
	if !almost(b.OrthogonalityDefect(), 1) {
		t.Errorf("defect = %v want 1", b.OrthogonalityDefect())
	}
	// skewed basis => defect > 1
	sk := NewBasis(NewVec(1, 0), NewVec(10, 1))
	if sk.OrthogonalityDefect() <= 1 {
		t.Errorf("skewed defect should exceed 1")
	}
}

func TestErrorsPaths(t *testing.T) {
	var empty Basis
	if _, err := empty.ShortestVector(); err != ErrEmpty {
		t.Errorf("expected ErrEmpty, got %v", err)
	}
	if _, err := empty.LLL(0.75); err != ErrEmpty {
		t.Errorf("expected ErrEmpty, got %v", err)
	}
	m := NewMatrix([][]float64{{1, 2, 3}})
	if _, err := m.Det(); err != ErrNotSquare {
		t.Errorf("expected ErrNotSquare, got %v", err)
	}
}

func TestPointAndCoordinates(t *testing.T) {
	b := NewBasis(NewVec(1, 1), NewVec(2, 0))
	p := b.Point([]int64{2, -1})
	if !p.Equal(NewVec(0, 2)) {
		t.Errorf("point = %v want [0 2]", p)
	}
	coords, err := b.Coordinates(NewVec(0, 2))
	if err != nil {
		t.Fatal(err)
	}
	if !coords.ApproxEqual(NewVec(2, -1), 1e-9) {
		t.Errorf("coords = %v want [2 -1]", coords)
	}
}
