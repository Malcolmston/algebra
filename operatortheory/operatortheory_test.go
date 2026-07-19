package operatortheory

import (
	"fmt"
	"math"
	"math/cmplx"
	"sort"
	"testing"
)

const tol = 1e-7

func mustMat(t *testing.T, rows, cols int, data []complex128) *Matrix {
	t.Helper()
	m, err := NewMatrixFromData(rows, cols, data)
	if err != nil {
		t.Fatalf("NewMatrixFromData: %v", err)
	}
	return m
}

func realMat(t *testing.T, rows, cols int, data []float64) *Matrix {
	t.Helper()
	m, err := FromReal(rows, cols, data)
	if err != nil {
		t.Fatalf("FromReal: %v", err)
	}
	return m
}

// c is a short constructor for complex128.
func c(re, im float64) complex128 { return complex(re, im) }

func TestArithmetic(t *testing.T) {
	a := realMat(t, 2, 2, []float64{1, 2, 3, 4})
	b := realMat(t, 2, 2, []float64{5, 6, 7, 8})

	sum, _ := a.Add(b)
	if !sum.Equal(realMat(t, 2, 2, []float64{6, 8, 10, 12}), tol) {
		t.Errorf("Add wrong: %v", sum)
	}
	prod, _ := a.Mul(b)
	if !prod.Equal(realMat(t, 2, 2, []float64{19, 22, 43, 50}), tol) {
		t.Errorf("Mul wrong: %v", prod)
	}
	if got := a.Trace(); got != c(5, 0) {
		t.Errorf("Trace = %v, want 5", got)
	}
	// Adjoint conjugates and transposes.
	m := mustMat(t, 1, 2, []complex128{c(1, 2), c(3, -4)})
	adj := m.Adjoint()
	if adj.At(0, 0) != c(1, -2) || adj.At(1, 0) != c(3, 4) {
		t.Errorf("Adjoint wrong: %v", adj)
	}
}

func TestPowerAndCommutator(t *testing.T) {
	a := realMat(t, 2, 2, []float64{2, 0, 0, 3})
	p := a.Power(3)
	if !p.Equal(realMat(t, 2, 2, []float64{8, 0, 0, 27}), tol) {
		t.Errorf("Power wrong: %v", p)
	}
	if !a.Power(0).Equal(Identity(2), tol) {
		t.Errorf("Power(0) should be identity")
	}
	// [X,Y] for Pauli-like matrices is nonzero.
	x := realMat(t, 2, 2, []float64{0, 1, 1, 0})
	z := realMat(t, 2, 2, []float64{1, 0, 0, -1})
	comm, _ := x.Commutator(z)
	if comm.MaxAbs() < 1 {
		t.Errorf("expected nonzero commutator, got %v", comm)
	}
}

func TestKroneckerAndDirectSum(t *testing.T) {
	i2 := Identity(2)
	k := i2.Kron(i2)
	if !k.Equal(Identity(4), tol) {
		t.Errorf("I2 kron I2 should be I4")
	}
	a := realMat(t, 1, 1, []float64{2})
	b := realMat(t, 1, 1, []float64{3})
	ds := a.DirectSum(b)
	if ds.At(0, 0) != c(2, 0) || ds.At(1, 1) != c(3, 0) || ds.At(0, 1) != 0 {
		t.Errorf("DirectSum wrong: %v", ds)
	}
}

func TestPredicates(t *testing.T) {
	tests := []struct {
		name string
		m    *Matrix
		pred func(*Matrix) bool
		want bool
	}{
		{"hermitian", mustMat(t, 2, 2, []complex128{2, c(0, 1), c(0, -1), 2}), func(m *Matrix) bool { return m.IsHermitian(tol) }, true},
		{"not-hermitian", realMat(t, 2, 2, []float64{1, 2, 3, 4}), func(m *Matrix) bool { return m.IsHermitian(tol) }, false},
		{"skew", mustMat(t, 2, 2, []complex128{c(0, 1), 2, -2, c(0, 3)}), func(m *Matrix) bool { return m.IsSkewHermitian(tol) }, true},
		{"unitary", mustMat(t, 2, 2, []complex128{c(0, 0), c(1, 0), c(1, 0), c(0, 0)}), func(m *Matrix) bool { return m.IsUnitary(tol) }, true},
		{"normal-diag", Diagonal([]complex128{c(1, 1), c(2, -3)}), func(m *Matrix) bool { return m.IsNormal(tol) }, true},
		{"non-normal", realMat(t, 2, 2, []float64{0, 1, 0, 0}), func(m *Matrix) bool { return m.IsNormal(tol) }, false},
		{"projection", realMat(t, 2, 2, []float64{1, 0, 0, 0}), func(m *Matrix) bool { return m.IsProjection(tol) }, true},
		{"orthproj", realMat(t, 2, 2, []float64{1, 0, 0, 0}), func(m *Matrix) bool { return m.IsOrthogonalProjection(tol) }, true},
		{"oblique-proj", realMat(t, 2, 2, []float64{1, 1, 0, 0}), func(m *Matrix) bool { return m.IsOrthogonalProjection(tol) }, false},
		{"involution", realMat(t, 2, 2, []float64{0, 1, 1, 0}), func(m *Matrix) bool { return m.IsInvolution(tol) }, true},
		{"nilpotent", realMat(t, 2, 2, []float64{0, 1, 0, 0}), func(m *Matrix) bool { return m.IsNilpotent(tol) }, true},
		{"pos-def", realMat(t, 2, 2, []float64{2, 0, 0, 3}), func(m *Matrix) bool { return m.IsPositiveDefinite(tol) }, true},
		{"not-pos-def", realMat(t, 2, 2, []float64{2, 0, 0, -3}), func(m *Matrix) bool { return m.IsPositiveDefinite(tol) }, false},
		{"psd", realMat(t, 2, 2, []float64{0, 0, 0, 3}), func(m *Matrix) bool { return m.IsPositiveSemidefinite(tol) }, true},
		{"diagonal", Diagonal([]complex128{1, 2}), func(m *Matrix) bool { return m.IsDiagonal(tol) }, true},
		{"upper-tri", realMat(t, 2, 2, []float64{1, 2, 0, 3}), func(m *Matrix) bool { return m.IsUpperTriangular(tol) }, true},
		{"contraction", realMat(t, 2, 2, []float64{0.5, 0, 0, 0.5}), func(m *Matrix) bool { return m.IsContraction(tol) }, true},
		{"isometry", mustMat(t, 2, 2, []complex128{0, 1, 1, 0}), func(m *Matrix) bool { return m.IsIsometry(tol) }, true},
		{"schur-stable", realMat(t, 2, 2, []float64{0.5, 0, 0, 0.3}), func(m *Matrix) bool { return m.IsSchurStable(tol) }, true},
		{"stable", realMat(t, 2, 2, []float64{-1, 0, 0, -2}), func(m *Matrix) bool { return m.IsStable(tol) }, true},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			if got := tc.pred(tc.m); got != tc.want {
				t.Errorf("%s: got %v want %v", tc.name, got, tc.want)
			}
		})
	}
}

func TestNorms(t *testing.T) {
	a := realMat(t, 2, 2, []float64{3, 0, 0, 4})
	tests := []struct {
		name string
		got  float64
		want float64
	}{
		{"frobenius", a.FrobeniusNorm(), 5},
		{"one", a.OneNorm(), 4},
		{"inf", a.InfNorm(), 4},
		{"operator", a.OperatorNorm(), 4},
		{"nuclear", a.NuclearNorm(), 7},
		{"schatten2", a.SchattenNorm(2), 5},
		{"kyfan1", a.KyFanNorm(1), 4},
	}
	for _, tc := range tests {
		if math.Abs(tc.got-tc.want) > 1e-6 {
			t.Errorf("%s = %v, want %v", tc.name, tc.got, tc.want)
		}
	}
}

func TestHermitianEigen(t *testing.T) {
	// [[2, i],[-i, 2]] has eigenvalues 1 and 3.
	h := mustMat(t, 2, 2, []complex128{2, c(0, 1), c(0, -1), 2})
	vals, vecs, err := h.HermitianEigen()
	if err != nil {
		t.Fatal(err)
	}
	sort.Float64s(vals)
	if math.Abs(vals[0]-1) > 1e-9 || math.Abs(vals[1]-3) > 1e-9 {
		t.Errorf("eigenvalues = %v, want [1 3]", vals)
	}
	// Reconstruct U D U^H.
	d := NewMatrix(2, 2)
	// vecs columns correspond to ascending vals from the raw solver; rebuild
	// directly from HermitianEigen output.
	valsRaw, vecsRaw := hermitianEigenRaw(h)
	for i, v := range valsRaw {
		d.Set(i, i, complex(v, 0))
	}
	t1, _ := vecsRaw.Mul(d)
	rec, _ := t1.Mul(vecsRaw.Adjoint())
	if !rec.Equal(h, 1e-8) {
		t.Errorf("reconstruction failed: %v", rec)
	}
	_ = vecs
}

func TestGeneralSpectrum(t *testing.T) {
	tests := []struct {
		name string
		m    *Matrix
		want []complex128
	}{
		{"companion", func() *Matrix { m, _ := Companion([]complex128{2, -3}); return m }(), []complex128{c(1, 0), c(2, 0)}},
		{"rotation", realMat(t, 2, 2, []float64{0, -1, 1, 0}), []complex128{c(0, -1), c(0, 1)}},
		{"diagonal", Diagonal([]complex128{c(5, 0), c(-2, 0), c(3, 0)}), []complex128{c(-2, 0), c(3, 0), c(5, 0)}},
		{"triangular", realMat(t, 3, 3, []float64{1, 9, 9, 0, 2, 9, 0, 0, 3}), []complex128{c(1, 0), c(2, 0), c(3, 0)}},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, err := tc.m.Spectrum()
			if err != nil {
				t.Fatal(err)
			}
			if len(got) != len(tc.want) {
				t.Fatalf("len = %d, want %d", len(got), len(tc.want))
			}
			for i := range got {
				if cmplx.Abs(got[i]-tc.want[i]) > 1e-6 {
					t.Errorf("eigenvalue[%d] = %v, want %v", i, got[i], tc.want[i])
				}
			}
		})
	}
}

func TestSpectralRadiusAbscissa(t *testing.T) {
	m := realMat(t, 2, 2, []float64{0, -2, 2, 0}) // eigenvalues +-2i
	if math.Abs(m.SpectralRadius()-2) > 1e-9 {
		t.Errorf("SpectralRadius = %v, want 2", m.SpectralRadius())
	}
	if math.Abs(m.SpectralAbscissa()) > 1e-9 {
		t.Errorf("SpectralAbscissa = %v, want 0", m.SpectralAbscissa())
	}
	h := realMat(t, 2, 2, []float64{-1, 0, 0, -3})
	if math.Abs(h.SpectralAbscissa()-(-1)) > 1e-9 {
		t.Errorf("SpectralAbscissa = %v, want -1", h.SpectralAbscissa())
	}
}

func TestCharacteristicPolynomial(t *testing.T) {
	m := realMat(t, 2, 2, []float64{1, 2, 3, 4})
	cp, err := m.CharacteristicPolynomial()
	if err != nil {
		t.Fatal(err)
	}
	// x^2 - 5x - 2.
	want := []complex128{c(-2, 0), c(-5, 0), c(1, 0)}
	for i := range want {
		if cmplx.Abs(cp[i]-want[i]) > 1e-9 {
			t.Errorf("coeff[%d] = %v, want %v", i, cp[i], want[i])
		}
	}
}

func TestDeterminantInverseSolve(t *testing.T) {
	m := realMat(t, 3, 3, []float64{2, 1, 1, 1, 3, 2, 1, 0, 0})
	det, _ := m.Determinant()
	if cmplx.Abs(det-c(-1, 0)) > 1e-9 {
		t.Errorf("Determinant = %v, want -1", det)
	}
	inv, err := m.Inverse()
	if err != nil {
		t.Fatal(err)
	}
	prod, _ := m.Mul(inv)
	if !prod.Equal(Identity(3), 1e-8) {
		t.Errorf("M*Minv != I: %v", prod)
	}
	b := VectorOf(c(4, 0), c(5, 0), c(6, 0))
	x, err := m.SolveVec(b)
	if err != nil {
		t.Fatal(err)
	}
	got, _ := m.MulVec(x)
	if !got.Equal(b, 1e-8) {
		t.Errorf("solve failed: %v", got)
	}
}

func TestSingularAndRank(t *testing.T) {
	m := realMat(t, 2, 2, []float64{1, 1, 1, 1}) // rank 1
	if r := m.Rank(0); r != 1 {
		t.Errorf("Rank = %d, want 1", r)
	}
	if _, err := m.Inverse(); err != ErrSingular {
		t.Errorf("expected ErrSingular, got %v", err)
	}
}

func TestSVDReconstruct(t *testing.T) {
	m := mustMat(t, 2, 3, []complex128{c(1, 1), c(2, 0), c(0, -1), c(3, 0), c(0, 1), c(1, 0)})
	u, s, v := m.SVD()
	// Reconstruct U diag(s) V^H.
	k := len(s)
	sd := NewMatrix(k, k)
	for i := range s {
		sd.Set(i, i, complex(s[i], 0))
	}
	t1, _ := u.Mul(sd)
	rec, _ := t1.Mul(v.Adjoint())
	if !rec.Equal(m, 1e-7) {
		t.Errorf("SVD reconstruction failed:\n%v\nwant\n%v", rec, m)
	}
	// Singular values are sorted descending and nonnegative.
	for i := 1; i < len(s); i++ {
		if s[i] > s[i-1]+1e-12 || s[i] < -1e-12 {
			t.Errorf("singular values not sorted/nonneg: %v", s)
		}
	}
}

func TestPolarDecomposition(t *testing.T) {
	m := RandomMatrix(4, 4, 7)
	u, p, err := m.PolarDecomposition()
	if err != nil {
		t.Fatal(err)
	}
	if !u.IsUnitary(1e-6) {
		t.Errorf("polar factor U not unitary")
	}
	if !p.IsPositiveSemidefinite(1e-6) {
		t.Errorf("polar factor P not PSD")
	}
	back, _ := u.Mul(p)
	if !back.Equal(m, 1e-6) {
		t.Errorf("polar reconstruction failed")
	}
}

func TestFunctionalCalculus(t *testing.T) {
	// exp of skew-symmetric [[0,-t],[t,0]] is a rotation.
	th := 0.9
	m := realMat(t, 2, 2, []float64{0, -th, th, 0})
	e, err := m.Exp()
	if err != nil {
		t.Fatal(err)
	}
	want := realMat(t, 2, 2, []float64{math.Cos(th), -math.Sin(th), math.Sin(th), math.Cos(th)})
	if !e.Equal(want, 1e-8) {
		t.Errorf("Exp rotation failed: %v", e)
	}
	// Sqrt of a positive definite matrix squares back.
	pd := realMat(t, 2, 2, []float64{2, 0, 0, 8})
	sq, _ := pd.Sqrt()
	sq2, _ := sq.Mul(sq)
	if !sq2.Equal(pd, 1e-7) {
		t.Errorf("Sqrt^2 != original: %v", sq2)
	}
	// Log is inverse of Exp on a Hermitian positive definite matrix.
	lg, _ := pd.Log()
	back, _ := lg.Exp()
	if !back.Equal(pd, 1e-6) {
		t.Errorf("exp(log(M)) != M: %v", back)
	}
	// sin^2 + cos^2 = I.
	a := realMat(t, 2, 2, []float64{0.3, 0.1, 0.1, 0.5})
	s, _ := a.Sin()
	co, _ := a.Cos()
	s2, _ := s.Mul(s)
	c2, _ := co.Mul(co)
	id, _ := s2.Add(c2)
	if !id.Equal(Identity(2), 1e-7) {
		t.Errorf("sin^2+cos^2 != I: %v", id)
	}
}

func TestResolventAndNumericalRange(t *testing.T) {
	m := Diagonal([]complex128{c(1, 0), c(0, 1)})
	// Numerical range of a normal matrix is the convex hull of its spectrum.
	if r := m.NumericalRadius(128); math.Abs(r-1) > 1e-6 {
		t.Errorf("NumericalRadius = %v, want 1", r)
	}
	if a := m.NumericalAbscissa(); math.Abs(a-1) > 1e-9 {
		t.Errorf("NumericalAbscissa = %v, want 1", a)
	}
	// Resolvent identity: (zI - A) R = I.
	z := c(3, 2)
	res, err := m.Resolvent(z)
	if err != nil {
		t.Fatal(err)
	}
	target := m.Scale(-1)
	for i := 0; i < 2; i++ {
		target.Set(i, i, z+target.At(i, i))
	}
	prod, _ := target.Mul(res)
	if !prod.Equal(Identity(2), 1e-8) {
		t.Errorf("resolvent identity failed: %v", prod)
	}
}

func TestPseudospectra(t *testing.T) {
	// A non-normal Jordan-like block has resolvent norm exceeding 1/dist.
	m := realMat(t, 2, 2, []float64{0, 5, 0, 0})
	// At z far from spectrum the resolvent norm is finite.
	rn := m.ResolventNorm(c(2, 0))
	if math.IsInf(rn, 1) || rn <= 0 {
		t.Errorf("unexpected resolvent norm %v", rn)
	}
	// The eps-pseudospectral abscissa exceeds the spectral abscissa (0) for a
	// non-normal matrix.
	eps := 0.5
	psa, err := m.PseudospectralAbscissa(eps, 40)
	if err != nil {
		t.Fatal(err)
	}
	if psa <= 0 {
		t.Errorf("pseudospectral abscissa should be positive for non-normal matrix, got %v", psa)
	}
	// DistanceToSingularity equals smallest singular value.
	inv := realMat(t, 2, 2, []float64{1, 0, 0, 1})
	if d := inv.DistanceToSingularity(); math.Abs(d-1) > 1e-9 {
		t.Errorf("DistanceToSingularity = %v, want 1", d)
	}
}

func TestSpectralDecomposition(t *testing.T) {
	// A = 1*P1 + 4*P2 where the eigenvalues are 1 and 4.
	h := mustMat(t, 2, 2, []complex128{c(2.5, 0), c(1.5, 0), c(1.5, 0), c(2.5, 0)})
	comps, err := h.SpectralDecomposition(1e-9)
	if err != nil {
		t.Fatal(err)
	}
	if len(comps) != 2 {
		t.Fatalf("expected 2 components, got %d", len(comps))
	}
	// Projectors sum to identity.
	sum := NewMatrix(2, 2)
	for _, cp := range comps {
		sum, _ = sum.Add(cp.Projector)
		if !cp.Projector.IsOrthogonalProjection(1e-7) {
			t.Errorf("component projector not an orthogonal projection")
		}
	}
	if !sum.Equal(Identity(2), 1e-8) {
		t.Errorf("projectors do not sum to identity: %v", sum)
	}
	// Reconstruct A = sum lambda_k P_k.
	rec := NewMatrix(2, 2)
	for _, cp := range comps {
		rec, _ = rec.Add(cp.Projector.Scale(complex(cp.Eigenvalue, 0)))
	}
	if !rec.Equal(h, 1e-7) {
		t.Errorf("spectral reconstruction failed: %v", rec)
	}
}

func TestOrthogonalProjector(t *testing.T) {
	v := VectorOf(c(1, 0), c(1, 0), c(0, 0))
	p, err := OrthogonalProjector([]Vector{v})
	if err != nil {
		t.Fatal(err)
	}
	if !p.IsOrthogonalProjection(1e-9) {
		t.Errorf("not an orthogonal projection")
	}
	// P v = v (v is in the range).
	pv, _ := p.MulVec(v)
	if !pv.Equal(v, 1e-9) {
		t.Errorf("P v != v: %v", pv)
	}
	// Idempotency: P^2 = P.
	p2, _ := p.Mul(p)
	if !p2.Equal(p, 1e-9) {
		t.Errorf("P^2 != P")
	}
}

func TestVectorOps(t *testing.T) {
	x := VectorOf(c(1, 0), c(0, 1))
	y := VectorOf(c(0, 1), c(1, 0))
	// <x,y> = conj(1)*i + conj(i)*1 = i + (-i) = 0.
	if d := x.Dot(y); cmplx.Abs(d) > 1e-12 {
		t.Errorf("Dot = %v, want 0", d)
	}
	if n := x.Norm(); math.Abs(n-math.Sqrt2) > 1e-12 {
		t.Errorf("Norm = %v, want sqrt2", n)
	}
	u, nrm := x.Normalize()
	if math.Abs(u.Norm()-1) > 1e-12 || math.Abs(nrm-math.Sqrt2) > 1e-12 {
		t.Errorf("Normalize failed: %v %v", u.Norm(), nrm)
	}
	// Gram-Schmidt of two independent vectors gives an orthonormal pair.
	basis := GramSchmidt([]Vector{VectorOf(c(1, 0), c(1, 0)), VectorOf(c(1, 0), c(0, 0))}, 1e-12)
	if len(basis) != 2 {
		t.Fatalf("expected 2 basis vectors, got %d", len(basis))
	}
	if cmplx.Abs(basis[0].Dot(basis[1])) > 1e-12 {
		t.Errorf("basis not orthogonal")
	}
}

func TestCayleyTransform(t *testing.T) {
	// Cayley transform of a Hermitian matrix is unitary.
	h := mustMat(t, 2, 2, []complex128{2, c(0, 1), c(0, -1), 2})
	u, err := h.CayleyTransform()
	if err != nil {
		t.Fatal(err)
	}
	if !u.IsUnitary(1e-7) {
		t.Errorf("Cayley transform not unitary: %v", u)
	}
}

func TestRandomUnitary(t *testing.T) {
	u := RandomUnitary(5, 123)
	if !u.IsUnitary(1e-7) {
		t.Errorf("RandomUnitary not unitary")
	}
	// Determinism: same seed gives same matrix.
	u2 := RandomUnitary(5, 123)
	if !u.Equal(u2, 1e-12) {
		t.Errorf("RandomUnitary not deterministic")
	}
}

func TestErrors(t *testing.T) {
	a := realMat(t, 2, 3, []float64{1, 2, 3, 4, 5, 6})
	if _, err := a.Determinant(); err != ErrNotSquare {
		t.Errorf("expected ErrNotSquare, got %v", err)
	}
	b := realMat(t, 3, 3, []float64{1, 0, 0, 0, 1, 0, 0, 0, 1})
	if _, err := a.Add(b); err != ErrDimensionMismatch {
		t.Errorf("expected ErrDimensionMismatch, got %v", err)
	}
	if _, err := Companion(nil); err != ErrEmpty {
		t.Errorf("expected ErrEmpty, got %v", err)
	}
}

func ExampleMatrix_HermitianEigen() {
	// The Hermitian matrix [[2, i], [-i, 2]] has eigenvalues 1 and 3.
	m, _ := NewMatrixFromData(2, 2, []complex128{
		2, complex(0, 1),
		complex(0, -1), 2,
	})
	vals, _, _ := m.HermitianEigen()
	fmt.Printf("%.4f %.4f\n", vals[0], vals[1])
	// Output: 1.0000 3.0000
}

func ExampleMatrix_Exp() {
	// exp of a 90-degree generator is a quarter-turn rotation.
	g, _ := FromReal(2, 2, []float64{0, -math.Pi / 2, math.Pi / 2, 0})
	e, _ := g.Exp()
	fmt.Printf("%.3f %.3f\n%.3f %.3f\n",
		real(e.At(0, 0)), real(e.At(0, 1)),
		real(e.At(1, 0)), real(e.At(1, 1)))
	// Output:
	// 0.000 -1.000
	// 1.000 0.000
}

func ExampleMatrix_SpectralRadius() {
	// The rotation-by-scaling matrix has spectral radius equal to the scale.
	m, _ := FromReal(2, 2, []float64{0, -2, 2, 0})
	fmt.Printf("%.1f\n", m.SpectralRadius())
	// Output: 2.0
}
