package powerseries

import (
	"math"
	"testing"
)

func TestNumberSequences(t *testing.T) {
	approxSlice(t, "Catalan", CatalanNumbers(10), []float64{1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862})
	approxSlice(t, "Fibonacci", FibonacciNumbers(11), []float64{0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55})
	approxSlice(t, "Partition", PartitionNumbers(11), []float64{1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42})
	approxSlice(t, "Bell", BellNumbers(9), []float64{1, 1, 2, 5, 15, 52, 203, 877, 4140})
	approxSlice(t, "Derangement", DerangementNumbers(8), []float64{1, 0, 1, 2, 9, 44, 265, 1854})
	approxSlice(t, "Motzkin", MotzkinNumbers(9), []float64{1, 1, 2, 4, 9, 21, 51, 127, 323})
	approxSlice(t, "Bernoulli", BernoulliNumbers(9),
		[]float64{1, -0.5, 1.0 / 6, 0, -1.0 / 30, 0, 1.0 / 42, 0, -1.0 / 30})
	if math.Abs(Catalan(6)-132) > tol || math.Abs(Bell(6)-203) > tol {
		t.Errorf("scalar accessors wrong")
	}
	if math.Abs(PartitionCount(10)-42) > tol || math.Abs(Fibonacci(10)-55) > tol {
		t.Errorf("scalar accessors wrong")
	}
	if math.Abs(Bernoulli(4)-(-1.0/30)) > tol || math.Abs(Derangement(4)-9) > tol {
		t.Errorf("scalar accessors wrong")
	}
}

func TestGeneratingFunctions(t *testing.T) {
	n := 12
	approxSlice(t, "GeometricGF", GeometricGF(5).Coeffs(), []float64{1, 1, 1, 1, 1})
	approxSlice(t, "GeometricRatioGF", GeometricRatioGF(3, 4).Coeffs(), []float64{1, 3, 9, 27})
	// The Catalan OGF built with the series square root matches the numbers.
	approxSlice(t, "CatalanGF", CatalanGF(n).Coeffs(), CatalanNumbers(n))
	// The Fibonacci OGF matches the numbers.
	approxSlice(t, "FibonacciGF", FibonacciGF(n).Coeffs(), FibonacciNumbers(n))
	// The partition OGF product matches the pentagonal recurrence.
	approxSlice(t, "PartitionGF", PartitionGF(n).Coeffs(), PartitionNumbers(n))
	// (1+x)^alpha with alpha=4 is a polynomial with binomial coefficients.
	approxSlice(t, "BinomialGF", BinomialGF(4, 6).Coeffs(), []float64{1, 4, 6, 4, 1, 0})
	// e^x exponential generating function.
	approxSlice(t, "ExpGF", ExpGF(4).Coeffs(), []float64{1, 1, 0.5, 1.0 / 6})
	// Harmonic numbers as an OGF.
	approxSlice(t, "HarmonicGF", HarmonicGF(6).Coeffs(),
		[]float64{0, 1, 1.5, 1 + 0.5 + 1.0/3, 1 + 0.5 + 1.0/3 + 0.25, 1 + 0.5 + 1.0/3 + 0.25 + 0.2})
}

func TestGeneratingFunctionsViaEGF(t *testing.T) {
	n := 9
	// Recovering the integer sequence from an exponential generating function
	// multiplies the degree-k coefficient by k!.
	approxSlice(t, "BernoulliGF", BernoulliGF(n).SequenceFromEGF(), BernoulliNumbers(n))
	approxSlice(t, "BellGF", BellGF(n).SequenceFromEGF(), BellNumbers(n))
	approxSlice(t, "DerangementGF", DerangementGF(n).SequenceFromEGF(), DerangementNumbers(n))
}

func TestOGFtoEGFRoundTrip(t *testing.T) {
	seq := []float64{1, 3, 5, 7, 9, 11}
	ogf := OGFFromSequence(seq)
	egf := OGFtoEGF(ogf)
	// The EGF coefficient of degree k is seq[k]/k!.
	for k := range seq {
		if math.Abs(egf.Coeff(k)-seq[k]/powerseriesFactorial(k)) > tol {
			t.Errorf("OGFtoEGF[%d] wrong", k)
		}
	}
	// EGFtoOGF is the inverse.
	approxSlice(t, "round trip", EGFtoOGF(egf).Coeffs(), seq)
	// EGFFromSequence then SequenceFromEGF recovers the sequence.
	approxSlice(t, "EGF seq", EGFFromSequence(seq).SequenceFromEGF(), seq)
}

func TestLagrangeInversion(t *testing.T) {
	n := 9
	// w = t/(1-w) => w = t + w^2, giving shifted Catalan numbers.
	w := LagrangeInversion(GeometricGF(n), n)
	cat := CatalanNumbers(n)
	want := make([]float64, n)
	for m := 1; m < n; m++ {
		want[m] = cat[m-1]
	}
	approxSlice(t, "LagrangeInversion", w.Coeffs(), want)

	// The tree function T(t) = t·e^{T} solves w = t·phi(w) with phi = e^w and has
	// [t^m] T = m^{m-1}/m! (Cayley). Check via the exponential coefficients.
	tree := LagrangeInversion(ExpGF(n), n)
	treeSeq := tree.SequenceFromEGF() // multiplies by m!, giving m^{m-1}.
	wantTree := []float64{0, 1, 2, 9, 64, 625, 7776, 117649, 2097152}
	approxSlice(t, "tree function", treeSeq, wantTree)

	// LagrangeInversionApply with H(x)=x reproduces the plain inversion series.
	applied := LagrangeInversionApply(GeometricGF(n), Ident(n), n)
	approxSlice(t, "LagrangeInversionApply", applied.Coeffs(), want)
}

// BenchmarkReversion exercises the heaviest routine, the O(n^3) Lagrange-based
// compositional inverse of a dense truncated series.
func BenchmarkReversion(b *testing.B) {
	n := 64
	f := Ident(n).Sin() // odd dense series with unit linear term
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = f.Reversion()
	}
}
