package randommatrix

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"
)

const tol = 1e-9

func approx(a, b, t float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		return false
	}
	return math.Abs(a-b) <= t
}

func TestCatalanAndBinomial(t *testing.T) {
	cats := []struct {
		n    int
		want float64
	}{
		{0, 1}, {1, 1}, {2, 2}, {3, 5}, {4, 14}, {5, 42}, {6, 132},
	}
	for _, c := range cats {
		if got := CatalanNumber(c.n); got != c.want {
			t.Errorf("CatalanNumber(%d)=%v want %v", c.n, got, c.want)
		}
	}
	binoms := []struct {
		n, k int
		want float64
	}{
		{5, 2, 10}, {6, 3, 20}, {10, 0, 1}, {10, 10, 1}, {4, 5, 0}, {7, 2, 21},
	}
	for _, b := range binoms {
		if got := BinomialCoefficient(b.n, b.k); got != b.want {
			t.Errorf("BinomialCoefficient(%d,%d)=%v want %v", b.n, b.k, got, b.want)
		}
	}
	if NonCrossingPartitionsCount(4) != 14 {
		t.Errorf("NonCrossingPartitionsCount(4) should be Catalan(4)=14")
	}
}

func TestSemicircleLaw(t *testing.T) {
	if !approx(SemicircleRadius(1), 2, tol) {
		t.Errorf("radius for variance 1 should be 2")
	}
	if !approx(SemicircleVariance(2), 1, tol) {
		t.Errorf("variance for radius 2 should be 1")
	}
	tests := []struct {
		name string
		got  float64
		want float64
		tol  float64
	}{
		{"density@0,R=2", SemicircleDensity(0, 2), 1 / math.Pi, tol},
		{"peak,R=2", SemicirclePeak(2), 1 / math.Pi, tol},
		{"density outside", SemicircleDensity(3, 2), 0, tol},
		{"cdf@0", SemicircleCDF(0, 2), 0.5, tol},
		{"cdf@R", SemicircleCDF(2, 2), 1, tol},
		{"cdf@-R", SemicircleCDF(-2, 2), 0, tol},
		{"moment0", SemicircleMoment(0, 2), 1, tol},
		{"moment2(var)", SemicircleMoment(2, 2), 1, tol},
		{"moment4", SemicircleMoment(4, 2), 2, tol},
		{"moment6", SemicircleMoment(6, 2), 5, tol},
		{"odd moment", SemicircleMoment(3, 2), 0, tol},
	}
	for _, tc := range tests {
		if !approx(tc.got, tc.want, tc.tol) {
			t.Errorf("%s: got %v want %v", tc.name, tc.got, tc.want)
		}
	}
	// Density integrates to one over the support.
	integral := simpson(func(x float64) float64 { return SemicircleDensity(x, 2) }, -2, 2, 4000)
	if !approx(integral, 1, 1e-3) {
		t.Errorf("semicircle density integral = %v want ~1", integral)
	}
}

func TestSemicircleStieltjes(t *testing.T) {
	// For radius 2 (variance 1) the Stieltjes transform satisfies
	// sigma2 m^2 + z m + 1 = 0, i.e. m^2 + z m + 1 = 0.
	z := complex(0.3, 0.7)
	m := SemicircleStieltjes(z, 2)
	res := m*m + z*m + 1
	if cmplx.Abs(res) > 1e-9 {
		t.Errorf("Stieltjes does not satisfy self-consistency: residual %v", cmplx.Abs(res))
	}
	if imag(m) <= 0 {
		t.Errorf("for Im z>0 the Stieltjes transform must have Im m>0, got %v", m)
	}
}

func TestMarchenkoPastur(t *testing.T) {
	c, s2 := 0.5, 1.0
	lo, hi := MarchenkoPasturSupport(c, s2)
	if !approx(lo, (1-math.Sqrt(0.5))*(1-math.Sqrt(0.5)), tol) {
		t.Errorf("lower edge wrong: %v", lo)
	}
	if !approx(hi, (1+math.Sqrt(0.5))*(1+math.Sqrt(0.5)), tol) {
		t.Errorf("upper edge wrong: %v", hi)
	}
	if !approx(MarchenkoPasturMean(c, s2), 1, tol) {
		t.Errorf("mean should be sigma2=1")
	}
	if !approx(MarchenkoPasturVariance(c, s2), 0.5, tol) {
		t.Errorf("variance should be c*sigma2^2=0.5")
	}
	if !approx(MarchenkoPasturMoment(1, c, s2), 1, tol) {
		t.Errorf("first moment should be 1")
	}
	if !approx(MarchenkoPasturMoment(2, c, s2), 1.5, tol) {
		t.Errorf("second moment should be 1+c=1.5")
	}
	if MarchenkoPasturHasAtom(0.5) {
		t.Errorf("c<1 must not have atom")
	}
	if !MarchenkoPasturHasAtom(2.0) {
		t.Errorf("c>1 must have atom")
	}
	if !approx(MarchenkoPasturAtomMass(2.0), 0.5, tol) {
		t.Errorf("atom mass for c=2 should be 0.5")
	}
	// CDF over full bulk integrates to one (no atom for c<1).
	if got := MarchenkoPasturCDF(hi, c, s2); !approx(got, 1, 5e-3) {
		t.Errorf("MP CDF at upper edge = %v want ~1", got)
	}
	// Density integrates to first moment 1.
	moment0 := simpson(func(x float64) float64 { return MarchenkoPasturDensity(x, c, s2) }, lo, hi, 6000)
	if !approx(moment0, 1, 5e-3) {
		t.Errorf("MP density mass = %v want ~1", moment0)
	}
}

func TestMPStieltjes(t *testing.T) {
	// The MP Stieltjes transform solves c z m^2 + (z + c - 1) m + 1 = 0.
	c := 0.4
	z := complex(0.5, 0.6)
	m := MarchenkoPasturStieltjes(z, c)
	res := complex(c, 0)*z*m*m + (z+complex(c-1, 0))*m + 1
	if cmplx.Abs(res) > 1e-9 {
		t.Errorf("MP Stieltjes residual too large: %v", cmplx.Abs(res))
	}
}

func TestEigSymmetric(t *testing.T) {
	tests := []struct {
		name string
		mat  *Matrix
		want []float64
	}{
		{"diagonal", Diag([]float64{2, 3, -1}), []float64{-1, 2, 3}},
		{"2x2 sym", NewMatrixFromRows([][]float64{{2, 1}, {1, 2}}), []float64{1, 3}},
		{"swap", NewMatrixFromRows([][]float64{{0, 1}, {1, 0}}), []float64{-1, 1}},
		{"3x3", NewMatrixFromRows([][]float64{{2, 0, 0}, {0, 3, 4}, {0, 4, 9}}), []float64{1, 2, 11}},
	}
	for _, tc := range tests {
		got, err := EigenvaluesSymmetric(tc.mat)
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		if len(got) != len(tc.want) {
			t.Fatalf("%s: got %d eigenvalues want %d", tc.name, len(got), len(tc.want))
		}
		for i := range got {
			if !approx(got[i], tc.want[i], 1e-8) {
				t.Errorf("%s: eig[%d]=%v want %v", tc.name, i, got[i], tc.want[i])
			}
		}
	}
	// Eigenvector correctness: A v = lambda v.
	A := NewMatrixFromRows([][]float64{{2, 1}, {1, 2}})
	e, _ := EigSymmetric(A)
	for j := 0; j < 2; j++ {
		v := e.Vectors.Col(j)
		Av, _ := A.MulVec(v)
		for i := range v {
			if !approx(Av[i], e.Values[j]*v[i], 1e-8) {
				t.Errorf("eigenvector check failed at column %d", j)
			}
		}
	}
}

func TestEigHermitian(t *testing.T) {
	// Pauli Z.
	z := NewCMatrix(2, 2)
	z.Set(0, 0, 1)
	z.Set(1, 1, -1)
	vals, err := EigenvaluesHermitian(z)
	if err != nil {
		t.Fatal(err)
	}
	if len(vals) != 2 || !approx(vals[0], -1, 1e-8) || !approx(vals[1], 1, 1e-8) {
		t.Errorf("Pauli Z eigenvalues = %v want [-1 1]", vals)
	}
	// Pauli Y.
	y := NewCMatrix(2, 2)
	y.Set(0, 1, complex(0, -1))
	y.Set(1, 0, complex(0, 1))
	if !y.IsHermitian(1e-12) {
		t.Errorf("Pauli Y should be Hermitian")
	}
	vy, _ := EigenvaluesHermitian(y)
	if !approx(vy[0], -1, 1e-8) || !approx(vy[1], 1, 1e-8) {
		t.Errorf("Pauli Y eigenvalues = %v want [-1 1]", vy)
	}
}

func TestMatrixOps(t *testing.T) {
	a := NewMatrixFromRows([][]float64{{1, 2}, {3, 4}})
	if a.Trace() != 5 {
		t.Errorf("trace = %v want 5", a.Trace())
	}
	at := a.Transpose()
	if at.At(0, 1) != 3 || at.At(1, 0) != 2 {
		t.Errorf("transpose wrong")
	}
	b := NewMatrixFromRows([][]float64{{0, 1}, {1, 0}})
	p, _ := a.Mul(b)
	if p.At(0, 0) != 2 || p.At(0, 1) != 1 {
		t.Errorf("product wrong: %v", p)
	}
	if !approx(a.FrobeniusNorm(), math.Sqrt(30), tol) {
		t.Errorf("frobenius wrong: %v", a.FrobeniusNorm())
	}
	sym := a.Symmetrize()
	if !sym.IsSymmetric(tol) {
		t.Errorf("symmetrize did not produce symmetric matrix")
	}
	if sym.At(0, 1) != 2.5 {
		t.Errorf("symmetric part wrong: %v", sym.At(0, 1))
	}
}

func TestEnsemblesStructure(t *testing.T) {
	goe := GOE(20, 12345)
	if !goe.IsSymmetric(1e-12) {
		t.Errorf("GOE must be symmetric")
	}
	gue := GUE(15, 999)
	if !gue.IsHermitian(1e-12) {
		t.Errorf("GUE must be Hermitian")
	}
	gse := GSE(4, 7)
	if !gse.IsHermitian(1e-12) {
		t.Errorf("GSE must be Hermitian")
	}
	// GSE eigenvalues are doubly degenerate.
	ev, _ := EigenvaluesHermitian(gse)
	if len(ev) != 8 {
		t.Fatalf("GSE(4) embeds to 8 eigenvalues, got %d", len(ev))
	}
	// Determinism: same seed => same spectrum.
	s1, _ := SampleGOE(30, 2024, math.Sqrt(30))
	s2, _ := SampleGOE(30, 2024, math.Sqrt(30))
	if len(s1) != len(s2) {
		t.Fatal("length mismatch")
	}
	for i := range s1 {
		if s1[i] != s2[i] {
			t.Fatalf("seeded GOE not reproducible at %d", i)
		}
	}
	// Different seed => different spectrum (almost surely).
	s3, _ := SampleGOE(30, 777, math.Sqrt(30))
	same := true
	for i := range s1 {
		if s1[i] != s3[i] {
			same = false
			break
		}
	}
	if same {
		t.Errorf("different seeds produced identical spectra")
	}
}

func TestSemicircleConvergence(t *testing.T) {
	// The GOE spectrum, rescaled by sqrt(n), should roughly fill [-2, 2] and
	// have second moment close to 1.
	n := 400
	vals, _ := SampleGOE(n, 555, math.Sqrt(float64(n)))
	if lo := vals[0]; lo < -2.6 || lo > -1.4 {
		t.Errorf("smallest rescaled eigenvalue out of range: %v", lo)
	}
	if hi := vals[len(vals)-1]; hi > 2.6 || hi < 1.4 {
		t.Errorf("largest rescaled eigenvalue out of range: %v", hi)
	}
	m2 := SpectralMoment(vals, 2)
	if !approx(m2, 1, 0.15) {
		t.Errorf("second spectral moment %v not close to semicircle value 1", m2)
	}
}

func TestWishartConvergence(t *testing.T) {
	// Sample covariance eigenvalues should sit inside the MP support.
	n, p := 200, 400
	c := float64(n) / float64(p)
	vals, _ := SampleWishartReal(n, p, 42)
	lo, hi := MarchenkoPasturSupport(c, 1)
	if vals[0] < lo-0.2 {
		t.Errorf("smallest eigenvalue %v below MP support %v", vals[0], lo)
	}
	if vals[len(vals)-1] > hi+0.2 {
		t.Errorf("largest eigenvalue %v above MP support %v", vals[len(vals)-1], hi)
	}
	if !approx(SpectralMean(vals), 1, 0.05) {
		t.Errorf("mean eigenvalue %v not close to 1", SpectralMean(vals))
	}
}

func TestSpacingLaws(t *testing.T) {
	if !approx(PoissonSpacingCDF(math.Ln2), 0.5, tol) {
		t.Errorf("Poisson median should be ln2")
	}
	if !approx(PoissonSpacingDensity(0), 1, tol) {
		t.Errorf("Poisson density at 0 should be 1")
	}
	// GOE surmise matches the closed form.
	for _, s := range []float64{0.2, 0.5, 1.0, 1.7} {
		if !approx(WignerSurmise(s, 1), GOESpacingDensity(s), 1e-10) {
			t.Errorf("WignerSurmise(beta=1) != GOESpacingDensity at s=%v", s)
		}
		if !approx(WignerSurmise(s, 2), GUESpacingDensity(s), 1e-10) {
			t.Errorf("WignerSurmise(beta=2) != GUESpacingDensity at s=%v", s)
		}
		if !approx(WignerSurmise(s, 4), GSESpacingDensity(s), 1e-9) {
			t.Errorf("WignerSurmise(beta=4) != GSESpacingDensity at s=%v", s)
		}
	}
	// GOE closed-form CDF.
	if !approx(GOESpacingCDF(1), 1-math.Exp(-math.Pi/4), tol) {
		t.Errorf("GOE CDF closed form wrong")
	}
	// Each surmise integrates to one and has unit mean.
	for _, beta := range []float64{1, 2, 4} {
		mass := simpson(func(s float64) float64 { return WignerSurmise(s, beta) }, 0, 12, 8000)
		if !approx(mass, 1, 2e-3) {
			t.Errorf("beta=%v surmise mass %v not ~1", beta, mass)
		}
		mean := simpson(func(s float64) float64 { return s * WignerSurmise(s, beta) }, 0, 12, 8000)
		if !approx(mean, 1, 2e-3) {
			t.Errorf("beta=%v surmise mean %v not ~1", beta, mean)
		}
	}
}

func TestSpacingRatios(t *testing.T) {
	eigs := []float64{0, 1, 3, 6, 10}
	gaps := NearestNeighbourSpacings(eigs)
	wantGaps := []float64{1, 2, 3, 4}
	for i := range gaps {
		if gaps[i] != wantGaps[i] {
			t.Errorf("gap %d = %v want %v", i, gaps[i], wantGaps[i])
		}
	}
	r := ConsecutiveSpacingRatios(eigs)
	// min/max of (1,2),(2,3),(3,4)
	want := []float64{0.5, 2.0 / 3.0, 0.75}
	for i := range r {
		if !approx(r[i], want[i], tol) {
			t.Errorf("ratio %d = %v want %v", i, r[i], want[i])
		}
		if r[i] < 0 || r[i] > 1 {
			t.Errorf("restricted ratio out of [0,1]: %v", r[i])
		}
	}
	if !approx(PoissonMeanRatio(), 2*math.Ln2-1, tol) {
		t.Errorf("Poisson mean ratio wrong")
	}
	// RatioSurmiseDensity integrates to one on [0,1].
	for _, beta := range []float64{1, 2, 4} {
		mass := simpson(func(x float64) float64 { return RatioSurmiseDensity(x, beta) }, 0, 1, 4000)
		if !approx(mass, 1, 5e-3) {
			t.Errorf("ratio surmise beta=%v mass %v not ~1", beta, mass)
		}
	}
}

func TestFreeProbability(t *testing.T) {
	// Semicircle of variance 1 has moments 0,1,0,2,0,5,...
	moments := []float64{0, 1, 0, 2, 0, 5}
	kappa := FreeCumulantsFromMoments(moments)
	want := []float64{0, 1, 0, 0, 0, 0}
	for i := range want {
		if !approx(kappa[i], want[i], 1e-9) {
			t.Errorf("free cumulant %d = %v want %v", i, kappa[i], want[i])
		}
	}
	// Round-trip.
	back := MomentsFromFreeCumulants(kappa)
	for i := range moments {
		if !approx(back[i], moments[i], 1e-9) {
			t.Errorf("roundtrip moment %d = %v want %v", i, back[i], moments[i])
		}
	}
	// Free convolution of two unit semicircles is a semicircle of variance 2:
	// second moment 2, fourth moment 2*(2)^2 = 8.
	conv := FreeConvolutionMoments(moments, moments)
	if !approx(conv[1], 2, 1e-9) {
		t.Errorf("free-convolved second moment = %v want 2", conv[1])
	}
	if !approx(conv[3], 8, 1e-8) {
		t.Errorf("free-convolved fourth moment = %v want 8", conv[3])
	}
	if !approx(FreeConvolutionSemicircleVariance(1, 1), 2, tol) {
		t.Errorf("semicircle free conv variance wrong")
	}
	if !approx(FreeConvolutionSemicircleRadius(2, 2), 2*math.Sqrt2, tol) {
		t.Errorf("semicircle free conv radius wrong")
	}
}

func TestSpectralMomentsAndTransforms(t *testing.T) {
	eigs := []float64{-1, 0, 1}
	if !approx(SpectralMoment(eigs, 1), 0, tol) {
		t.Errorf("mean should be 0")
	}
	if !approx(SpectralMoment(eigs, 2), 2.0/3.0, tol) {
		t.Errorf("second moment should be 2/3")
	}
	if !approx(SpectralVariance(eigs), 2.0/3.0, tol) {
		t.Errorf("variance should be 2/3")
	}
	z := complex(0, 1)
	st := EmpiricalStieltjes(eigs, z)
	ca := EmpiricalCauchy(eigs, z)
	if cmplx.Abs(st+ca) > 1e-12 {
		t.Errorf("Cauchy should be negative of Stieltjes")
	}
}

func TestHistogram(t *testing.T) {
	data := []float64{0.1, 0.2, 0.9, 1.5, 1.6, 1.7, 1.8}
	h := NewHistogram(data, 2, 0, 2)
	if h.NumBins() != 2 {
		t.Fatalf("bins = %d", h.NumBins())
	}
	if h.Counts[0] != 3 || h.Counts[1] != 4 {
		t.Errorf("counts = %v want [3 4]", h.Counts)
	}
	if h.Total != 7 {
		t.Errorf("total = %d want 7", h.Total)
	}
	// Density integrates to one.
	dens := h.Density()
	var mass float64
	for _, d := range dens {
		mass += d * h.BinWidth()
	}
	if !approx(mass, 1, tol) {
		t.Errorf("histogram density mass = %v want 1", mass)
	}
	if h.ModeBin() != 1 {
		t.Errorf("mode bin should be 1")
	}
}

func TestTracyWidom(t *testing.T) {
	if !TracyWidomSupportedBeta(1) || !TracyWidomSupportedBeta(2) || !TracyWidomSupportedBeta(4) {
		t.Errorf("beta 1,2,4 must be supported")
	}
	if TracyWidomSupportedBeta(3) {
		t.Errorf("beta 3 must not be supported")
	}
	// Tabulated moments.
	if !approx(TracyWidomMean(1), -1.2065335745820, 1e-9) {
		t.Errorf("TW1 mean wrong: %v", TracyWidomMean(1))
	}
	if !approx(TracyWidomVariance(2), 0.8131947928, 1e-3) {
		t.Errorf("TW2 variance wrong: %v", TracyWidomVariance(2))
	}
	// Approximation mean matches tabulated mean by construction.
	k, th, sh, _ := TracyWidomApproxParams(2)
	if !approx(k*th-sh, TracyWidomMean(2), 1e-9) {
		t.Errorf("approx mean k*theta-shift = %v want %v", k*th-sh, TracyWidomMean(2))
	}
	// CDF is a valid distribution function.
	if !approx(TracyWidomCDF(-8, 2), 0, 1e-3) {
		t.Errorf("TW2 CDF far left should be ~0")
	}
	if !approx(TracyWidomCDF(6, 2), 1, 1e-3) {
		t.Errorf("TW2 CDF far right should be ~1")
	}
	// Reference values (Prahofer-Spohn / Bornemann tables), matched to the
	// accuracy of the shifted-gamma approximation.
	refs := []struct {
		s, want, tol float64
	}{
		{-3.0, 0.0803, 0.01},
		{-2.0, 0.4132, 0.01},
		{-1.0, 0.8069, 0.01},
		{0.0, 0.9694, 0.01},
		{1.0, 0.9975, 0.01},
	}
	for _, r := range refs {
		if got := TracyWidomCDF(r.s, 2); !approx(got, r.want, r.tol) {
			t.Errorf("TW2 CDF(%v)=%v want %v +/- %v", r.s, got, r.want, r.tol)
		}
	}
	// Monotonicity.
	prev := -1.0
	for s := -6.0; s <= 4.0; s += 0.25 {
		v := TracyWidomCDF(s, 1)
		if v < prev-1e-12 {
			t.Errorf("TW1 CDF not monotone at s=%v", s)
		}
		prev = v
	}
	// Density integrates to the CDF (fundamental theorem check).
	area := simpson(func(s float64) float64 { return TracyWidomDensity(s, 2) }, -8, 6, 8000)
	if !approx(area, 1, 5e-3) {
		t.Errorf("TW2 density mass %v not ~1", area)
	}
}

func TestSoftEdgeScaling(t *testing.T) {
	if !approx(TracyWidomEdgeCenter(64), 16, tol) {
		t.Errorf("edge center for n=64 should be 2*8=16")
	}
	if !approx(TracyWidomEdgeScale(64), math.Pow(64, -1.0/6.0), tol) {
		t.Errorf("edge scale wrong")
	}
	// Rescaling the center gives 0.
	if !approx(RescaleLargestEigenvalue(16, 64), 0, 1e-9) {
		t.Errorf("rescaling the center should give 0")
	}
	// Wishart Johnstone scaling for small n, p.
	n, p := 10, 20
	a := math.Sqrt(float64(n-1)) + math.Sqrt(float64(p))
	if !approx(WishartTracyWidomCenter(n, p), a*a, tol) {
		t.Errorf("Wishart TW center wrong")
	}
}

func TestErrors(t *testing.T) {
	rect := NewMatrix(2, 3)
	if _, err := EigenvaluesSymmetric(rect); err == nil {
		t.Errorf("expected error for non-square matrix")
	}
	a := NewMatrix(2, 2)
	b := NewMatrix(3, 3)
	if _, err := a.Plus(b); err == nil {
		t.Errorf("expected dimension mismatch error")
	}
}

func ExampleSemicircleMoment() {
	// Even moments of the semicircle law of radius 2 are the Catalan numbers.
	fmt.Printf("%.0f %.0f %.0f %.0f\n",
		SemicircleMoment(0, 2), SemicircleMoment(2, 2),
		SemicircleMoment(4, 2), SemicircleMoment(6, 2))
	// Output: 1 1 2 5
}

func ExampleGOE() {
	// GOE matrices are real symmetric, and generation is fully seeded.
	m := GOE(3, 42)
	fmt.Println(m.IsSymmetric(1e-12))
	// Output: true
}

func ExampleWignerSurmise() {
	// The Wigner surmise for the GUE (beta=2) at unit spacing.
	fmt.Printf("%.4f\n", WignerSurmise(1, 2))
	// Output: 0.9076
}
