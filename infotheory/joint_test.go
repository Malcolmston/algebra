package infotheory

import (
	"math"
	"testing"
)

func TestMarginals(t *testing.T) {
	j := [][]float64{{0.1, 0.2}, {0.3, 0.4}}
	px := MarginalX(j)
	py := MarginalY(j)
	if !approx(px[0], 0.3, tol) || !approx(px[1], 0.7, tol) {
		t.Errorf("MarginalX = %v, want [0.3 0.7]", px)
	}
	if !approx(py[0], 0.4, tol) || !approx(py[1], 0.6, tol) {
		t.Errorf("MarginalY = %v, want [0.4 0.6]", py)
	}
}

func TestJointIndependent(t *testing.T) {
	// Independent joint: p_ij = px_i py_j => I = 0, H(X,Y)=H(X)+H(Y).
	j := [][]float64{{0.25, 0.25}, {0.25, 0.25}}
	if got := MutualInformation(j); !approx(got, 0, tol) {
		t.Errorf("MutualInformation(independent) = %v, want 0", got)
	}
	if got := JointEntropy(j); !approx(got, 2, tol) {
		t.Errorf("JointEntropy(uniform 2x2) = %v, want 2", got)
	}
	if got := ConditionalEntropyYgivenX(j); !approx(got, 1, tol) {
		t.Errorf("H(Y|X) independent = %v, want 1", got)
	}
}

func TestJointPerfectlyCorrelated(t *testing.T) {
	// Diagonal joint: X == Y, both fair bits. H(X)=H(Y)=1, H(X,Y)=1, I=1.
	j := [][]float64{{0.5, 0}, {0, 0.5}}
	if got := JointEntropy(j); !approx(got, 1, tol) {
		t.Errorf("JointEntropy(diagonal) = %v, want 1", got)
	}
	if got := MutualInformation(j); !approx(got, 1, tol) {
		t.Errorf("MutualInformation(diagonal) = %v, want 1", got)
	}
	if got := ConditionalEntropyYgivenX(j); !approx(got, 0, tol) {
		t.Errorf("H(Y|X) diagonal = %v, want 0", got)
	}
	if got := VariationOfInformation(j); !approx(got, 0, tol) {
		t.Errorf("VI(diagonal) = %v, want 0", got)
	}
	if got := NormalizedMutualInformation(j); !approx(got, 1, tol) {
		t.Errorf("NMI(diagonal) = %v, want 1", got)
	}
}

func TestMutualInformationIdentity(t *testing.T) {
	// I(X;Y) = H(Y) - H(Y|X) must hold for any joint.
	j := [][]float64{{0.1, 0.2, 0.05}, {0.15, 0.1, 0.4}}
	i := MutualInformation(j)
	alt := Entropy(MarginalY(j)) - ConditionalEntropyYgivenX(j)
	if !approx(i, alt, 1e-9) {
		t.Errorf("I=%v but H(Y)-H(Y|X)=%v", i, alt)
	}
	// VI = H(X|Y) + H(Y|X).
	vi := VariationOfInformation(j)
	altVI := ConditionalEntropyXgivenY(j) + ConditionalEntropyYgivenX(j)
	if !approx(vi, altVI, 1e-9) {
		t.Errorf("VI=%v but H(X|Y)+H(Y|X)=%v", vi, altVI)
	}
}

func TestKLDivergence(t *testing.T) {
	p := []float64{0.5, 0.5}
	q := []float64{0.25, 0.75}
	if got := KLDivergence(p, p); !approx(got, 0, tol) {
		t.Errorf("KL(p||p) = %v, want 0", got)
	}
	if got := KLDivergence(p, q); !approx(got, 0.2075187496394219, 1e-9) {
		t.Errorf("KL(p||q) = %v, want 0.20751...", got)
	}
	// KL is infinite when support of p exceeds support of q.
	if got := KLDivergence([]float64{0.5, 0.5}, []float64{1, 0}); !math.IsInf(got, 1) {
		t.Errorf("KL with zero q = %v, want +Inf", got)
	}
	if got := KLDivergence([]float64{0.5, 0.5}, []float64{0.5}); !math.IsNaN(got) {
		t.Errorf("KL mismatched lengths = %v, want NaN", got)
	}
}

func TestCrossEntropy(t *testing.T) {
	p := []float64{0.5, 0.5}
	q := []float64{0.25, 0.75}
	// H(p,q) = H(p) + KL(p||q).
	if got := CrossEntropy(p, q); !approx(got, 1.207518749639422, 1e-9) {
		t.Errorf("CrossEntropy = %v, want 1.20751...", got)
	}
	if got := CrossEntropy(p, p); !approx(got, Entropy(p), tol) {
		t.Errorf("CrossEntropy(p,p) = %v, want H(p) = %v", got, Entropy(p))
	}
}

func TestJensenShannon(t *testing.T) {
	if got := JensenShannonDivergence([]float64{1, 0}, []float64{0, 1}); !approx(got, 1, tol) {
		t.Errorf("JSD(disjoint) = %v, want 1", got)
	}
	p := []float64{0.4, 0.6}
	if got := JensenShannonDivergence(p, p); !approx(got, 0, tol) {
		t.Errorf("JSD(p,p) = %v, want 0", got)
	}
	// JSD is symmetric.
	q := []float64{0.7, 0.3}
	if !approx(JensenShannonDivergence(p, q), JensenShannonDivergence(q, p), tol) {
		t.Error("JSD not symmetric")
	}
	if got := JensenShannonDistance([]float64{1, 0}, []float64{0, 1}); !approx(got, 1, tol) {
		t.Errorf("JS distance(disjoint) = %v, want 1", got)
	}
}
