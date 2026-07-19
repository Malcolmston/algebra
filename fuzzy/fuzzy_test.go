package fuzzy

import (
	"fmt"
	"math"
	"testing"
)

const tol = 1e-9

func approx(a, b, t float64) bool { return math.Abs(a-b) <= t }

func TestTriangular(t *testing.T) {
	tests := []struct {
		name       string
		x, a, b, c float64
		want       float64
	}{
		{"left-foot", 0, 0, 5, 10, 0},
		{"peak", 5, 0, 5, 10, 1},
		{"rise-mid", 2.5, 0, 5, 10, 0.5},
		{"fall-mid", 7.5, 0, 5, 10, 0.5},
		{"right-foot", 10, 0, 5, 10, 0},
		{"outside-left", -1, 0, 5, 10, 0},
		{"outside-right", 11, 0, 5, 10, 0},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := TriangularAt(tc.x, tc.a, tc.b, tc.c)
			if !approx(got, tc.want, tol) {
				t.Fatalf("TriangularAt(%v)=%v want %v", tc.x, got, tc.want)
			}
		})
	}
}

func TestTrapezoidal(t *testing.T) {
	f := Trapezoidal(1, 3, 6, 8)
	tests := []struct {
		x, want float64
	}{
		{0, 0}, {1, 0}, {2, 0.5}, {3, 1}, {4.5, 1}, {6, 1}, {7, 0.5}, {8, 0}, {9, 0},
	}
	for _, tc := range tests {
		if got := f(tc.x); !approx(got, tc.want, tol) {
			t.Errorf("Trapezoidal(%v)=%v want %v", tc.x, got, tc.want)
		}
	}
}

func TestGaussianAndBell(t *testing.T) {
	if got := GaussianAt(0, 0, 1); !approx(got, 1, tol) {
		t.Errorf("GaussianAt center = %v want 1", got)
	}
	if got := GaussianAt(1, 0, 1); !approx(got, math.Exp(-0.5), tol) {
		t.Errorf("GaussianAt = %v want %v", got, math.Exp(-0.5))
	}
	// bell at center is 1, at |x-c|=a is 0.5
	if got := BellAt(2, 2, 4, 2); !approx(got, 1, tol) {
		t.Errorf("BellAt center = %v want 1", got)
	}
	if got := BellAt(4, 2, 4, 2); !approx(got, 0.5, tol) {
		t.Errorf("BellAt half point = %v want 0.5", got)
	}
}

func TestSigmoid(t *testing.T) {
	if got := SigmoidAt(0, 1, 0); !approx(got, 0.5, tol) {
		t.Errorf("SigmoidAt inflection = %v want 0.5", got)
	}
	if got := SigmoidAt(100, 1, 0); !approx(got, 1, 1e-6) {
		t.Errorf("SigmoidAt far right = %v want ~1", got)
	}
}

func TestSZShape(t *testing.T) {
	if got := SShapeAt(0, 0, 4); got != 0 {
		t.Errorf("SShape foot = %v want 0", got)
	}
	if got := SShapeAt(2, 0, 4); !approx(got, 0.5, tol) {
		t.Errorf("SShape mid = %v want 0.5", got)
	}
	if got := SShapeAt(4, 0, 4); !approx(got, 1, tol) {
		t.Errorf("SShape shoulder = %v want 1", got)
	}
	// Z is complement of S
	if got := ZShapeAt(2, 0, 4); !approx(got, 0.5, tol) {
		t.Errorf("ZShape mid = %v want 0.5", got)
	}
}

func TestTNorms(t *testing.T) {
	a, b := 0.4, 0.7
	tests := []struct {
		name string
		got  float64
		want float64
	}{
		{"min", TNormMin(a, b), 0.4},
		{"product", TNormProduct(a, b), 0.28},
		{"lukasiewicz", TNormLukasiewicz(a, b), 0.1},
		{"lukasiewicz-zero", TNormLukasiewicz(0.2, 0.3), 0},
		{"drastic", TNormDrastic(a, b), 0},
		{"drastic-one", TNormDrastic(1, b), 0.7},
		{"hamacher1==product", TNormHamacher(1)(a, b), 0.28},
		{"einstein==hamacher2", TNormEinstein(a, b), TNormHamacher(2)(a, b)},
		{"nilmin", TNormNilpotentMin(a, b), 0.4},
		{"nilmin-zero", TNormNilpotentMin(0.2, 0.3), 0},
	}
	for _, tc := range tests {
		if !approx(tc.got, tc.want, 1e-9) {
			t.Errorf("%s = %v want %v", tc.name, tc.got, tc.want)
		}
	}
}

func TestTConorms(t *testing.T) {
	a, b := 0.4, 0.7
	tests := []struct {
		name string
		got  float64
		want float64
	}{
		{"max", TConormMax(a, b), 0.7},
		{"prob", TConormProbabilistic(a, b), 0.82},
		{"lukasiewicz", TConormLukasiewicz(a, b), 1},
		{"lukasiewicz-sum", TConormLukasiewicz(0.2, 0.3), 0.5},
		{"drastic", TConormDrastic(a, b), 1},
		{"drastic-zero", TConormDrastic(0, b), 0.7},
	}
	for _, tc := range tests {
		if !approx(tc.got, tc.want, 1e-9) {
			t.Errorf("%s = %v want %v", tc.name, tc.got, tc.want)
		}
	}
	// De Morgan duality: dual conorm of min t-norm is max
	dual := DualTConorm(TNormMin)
	if !approx(dual(a, b), TConormMax(a, b), tol) {
		t.Errorf("DualTConorm(min) = %v want max %v", dual(a, b), TConormMax(a, b))
	}
	// probabilistic sum is dual of product
	dp := DualTConorm(TNormProduct)
	if !approx(dp(a, b), TConormProbabilistic(a, b), tol) {
		t.Errorf("dual product = %v want prob sum %v", dp(a, b), TConormProbabilistic(a, b))
	}
}

func TestComplements(t *testing.T) {
	if got := ComplementStandard(0.3); !approx(got, 0.7, tol) {
		t.Errorf("standard complement = %v want 0.7", got)
	}
	// Sugeno with lambda 0 == standard
	if got := ComplementSugeno(0)(0.3); !approx(got, 0.7, tol) {
		t.Errorf("sugeno lambda0 = %v want 0.7", got)
	}
	// Yager with w 1 == standard
	if got := ComplementYager(1)(0.3); !approx(got, 0.7, tol) {
		t.Errorf("yager w1 = %v want 0.7", got)
	}
	// equilibrium of standard complement is 0.5
	if got := ComplementYager(2)(math.Sqrt(0.5)); !approx(got, math.Sqrt(0.5), 1e-9) {
		t.Errorf("yager equilibrium = %v", got)
	}
}

func TestSetOpsAndCardinality(t *testing.T) {
	xs := []float64{1, 2, 3, 4}
	a, _ := NewSet(xs, []float64{0.2, 0.8, 1.0, 0.4})
	b, _ := NewSet(xs, []float64{0.5, 0.5, 0.5, 0.5})

	if got := a.Cardinality(); !approx(got, 2.4, tol) {
		t.Errorf("cardinality = %v want 2.4", got)
	}
	if got := a.Height(); !approx(got, 1.0, tol) {
		t.Errorf("height = %v want 1", got)
	}
	if !a.IsNormal() {
		t.Error("a should be normal")
	}
	un, err := a.Union(b)
	if err != nil {
		t.Fatal(err)
	}
	wantUn := []float64{0.5, 0.8, 1.0, 0.5}
	for i, w := range wantUn {
		if !approx(un.Mu[i], w, tol) {
			t.Errorf("union[%d] = %v want %v", i, un.Mu[i], w)
		}
	}
	in, _ := a.Intersection(b)
	wantIn := []float64{0.2, 0.5, 0.5, 0.4}
	for i, w := range wantIn {
		if !approx(in.Mu[i], w, tol) {
			t.Errorf("intersection[%d] = %v want %v", i, in.Mu[i], w)
		}
	}
	comp := a.Complement()
	wantC := []float64{0.8, 0.2, 0.0, 0.6}
	for i, w := range wantC {
		if !approx(comp.Mu[i], w, tol) {
			t.Errorf("complement[%d] = %v want %v", i, comp.Mu[i], w)
		}
	}
}

func TestAlphaCutSupportCore(t *testing.T) {
	xs := []float64{1, 2, 3, 4, 5}
	s, _ := NewSet(xs, []float64{0.1, 0.5, 1.0, 0.5, 0.1})
	cut := s.AlphaCut(0.5)
	want := []float64{2, 3, 4}
	if len(cut) != len(want) {
		t.Fatalf("alpha cut = %v want %v", cut, want)
	}
	for i := range want {
		if cut[i] != want[i] {
			t.Errorf("alpha cut[%d] = %v want %v", i, cut[i], want[i])
		}
	}
	strong := s.StrongAlphaCut(0.5)
	if len(strong) != 1 || strong[0] != 3 {
		t.Errorf("strong alpha cut = %v want [3]", strong)
	}
	if core := s.Core(); len(core) != 1 || core[0] != 3 {
		t.Errorf("core = %v want [3]", core)
	}
	if sup := s.Support(); len(sup) != 5 {
		t.Errorf("support = %v want all 5", sup)
	}
}

func TestHedges(t *testing.T) {
	xs := []float64{1, 2, 3}
	s, _ := NewSet(xs, []float64{0.25, 0.5, 1.0})
	very := s.Very()
	if !approx(very.Mu[0], 0.0625, tol) || !approx(very.Mu[1], 0.25, tol) {
		t.Errorf("very = %v", very.Mu)
	}
	som := s.Somewhat()
	if !approx(som.Mu[0], 0.5, tol) || !approx(som.Mu[1], math.Sqrt(0.5), tol) {
		t.Errorf("somewhat = %v", som.Mu)
	}
	// intensify pushes 0.25 down and keeps 1 at 1
	inten := s.Intensify()
	if !approx(inten.Mu[0], 0.125, tol) {
		t.Errorf("intensify low = %v want 0.125", inten.Mu[0])
	}
	if !approx(inten.Mu[2], 1, tol) {
		t.Errorf("intensify high = %v want 1", inten.Mu[2])
	}
}

func TestConvexity(t *testing.T) {
	xs := []float64{1, 2, 3, 4, 5}
	convex, _ := NewSet(xs, []float64{0.1, 0.6, 1.0, 0.6, 0.1})
	if !convex.IsConvex() {
		t.Error("expected convex")
	}
	nonconvex, _ := NewSet(xs, []float64{0.1, 0.8, 0.2, 0.9, 0.1})
	if nonconvex.IsConvex() {
		t.Error("expected non-convex")
	}
}

func TestRelationComposition(t *testing.T) {
	// classic max-min composition example
	x := []float64{1, 2}
	y := []float64{1, 2, 3}
	z := []float64{1, 2}
	r, _ := NewRelation(x, y, [][]float64{
		{0.3, 0.8, 0.5},
		{0.9, 0.2, 0.6},
	})
	s, _ := NewRelation(y, z, [][]float64{
		{0.4, 0.7},
		{0.6, 0.3},
		{0.5, 0.9},
	})
	comp, err := r.MaxMinComposition(s)
	if err != nil {
		t.Fatal(err)
	}
	// (1,1): max(min(.3,.4),min(.8,.6),min(.5,.5)) = max(.3,.6,.5)=.6
	if !approx(comp.M[0][0], 0.6, tol) {
		t.Errorf("comp[0][0] = %v want 0.6", comp.M[0][0])
	}
	// (1,2): max(min(.3,.7),min(.8,.3),min(.5,.9)) = max(.3,.3,.5)=.5
	if !approx(comp.M[0][1], 0.5, tol) {
		t.Errorf("comp[0][1] = %v want 0.5", comp.M[0][1])
	}
	// (2,1): max(min(.9,.4),min(.2,.6),min(.6,.5)) = max(.4,.2,.5)=.5
	if !approx(comp.M[1][0], 0.5, tol) {
		t.Errorf("comp[1][0] = %v want 0.5", comp.M[1][0])
	}
	// (2,2): max(min(.9,.7),min(.2,.3),min(.6,.9)) = max(.7,.2,.6)=.7
	if !approx(comp.M[1][1], 0.7, tol) {
		t.Errorf("comp[1][1] = %v want 0.7", comp.M[1][1])
	}
}

func TestMaxProductComposition(t *testing.T) {
	x := []float64{1}
	y := []float64{1, 2}
	z := []float64{1}
	r, _ := NewRelation(x, y, [][]float64{{0.5, 0.8}})
	s, _ := NewRelation(y, z, [][]float64{{0.4}, {0.5}})
	comp, _ := r.MaxProductComposition(s)
	// max(0.5*0.4, 0.8*0.5) = max(0.2,0.4)=0.4
	if !approx(comp.M[0][0], 0.4, tol) {
		t.Errorf("max-product = %v want 0.4", comp.M[0][0])
	}
}

func TestRelationProperties(t *testing.T) {
	u := []float64{1, 2, 3}
	r, _ := NewRelation(u, u, [][]float64{
		{1.0, 0.8, 0.3},
		{0.8, 1.0, 0.4},
		{0.3, 0.4, 1.0},
	})
	if !r.IsReflexive(tol) {
		t.Error("expected reflexive")
	}
	if !r.IsSymmetric(tol) {
		t.Error("expected symmetric")
	}
}

func TestDefuzzification(t *testing.T) {
	// symmetric triangle centered at 5 -> centroid, bisector, MOM all 5
	xs := Linspace(0, 10, 11)
	s := FromMF(Triangular(0, 5, 10), xs)
	c, err := s.Centroid()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(c, 5, 1e-9) {
		t.Errorf("centroid = %v want 5", c)
	}
	mom, _ := s.MeanOfMaxima()
	if !approx(mom, 5, tol) {
		t.Errorf("MOM = %v want 5", mom)
	}
	bis, _ := s.Bisector()
	if !approx(bis, 5, tol) {
		t.Errorf("bisector = %v want 5", bis)
	}

	// plateau set for SOM/LOM
	p, _ := NewSet([]float64{0, 1, 2, 3, 4}, []float64{0.2, 1.0, 1.0, 1.0, 0.2})
	som, _ := p.SmallestOfMaxima()
	lom, _ := p.LargestOfMaxima()
	if !approx(som, 1, tol) || !approx(lom, 3, tol) {
		t.Errorf("SOM=%v LOM=%v want 1,3", som, lom)
	}

	// zero area error
	zero, _ := NewSet([]float64{0, 1}, []float64{0, 0})
	if _, err := zero.Centroid(); err != ErrNoArea {
		t.Errorf("expected ErrNoArea got %v", err)
	}
}

func TestCentroidTrapz(t *testing.T) {
	// trapezoid symmetric about 5
	s := FromMF(Trapezoidal(1, 4, 6, 9), Linspace(0, 10, 101))
	c, err := s.CentroidTrapz()
	if err != nil {
		t.Fatal(err)
	}
	if !approx(c, 5, 1e-6) {
		t.Errorf("trapz centroid = %v want 5", c)
	}
}

func TestMamdaniTipping(t *testing.T) {
	// classic tipping problem, single input for a clean assertion
	service := NewVariable("service", 0, 10)
	service.AddTerm("poor", Triangular(0, 0, 5))
	service.AddTerm("good", Triangular(0, 5, 10))
	service.AddTerm("excellent", Triangular(5, 10, 10))

	tip := NewVariable("tip", 0, 30)
	tip.AddTerm("low", Triangular(0, 5, 10))
	tip.AddTerm("medium", Triangular(10, 15, 20))
	tip.AddTerm("high", Triangular(20, 25, 30))

	m := NewMamdani(tip)
	m.AddInput(service)
	m.AddRule(OpAnd, "low", 1, If("service", "poor"))
	m.AddRule(OpAnd, "medium", 1, If("service", "good"))
	m.AddRule(OpAnd, "high", 1, If("service", "excellent"))

	// excellent service -> high tip, near 25
	out, err := m.Infer(map[string]float64{"service": 10})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(out, 25, 0.5) {
		t.Errorf("tip for excellent = %v want ~25", out)
	}
	// good service (5) -> medium tip, near 15
	out2, _ := m.Infer(map[string]float64{"service": 5})
	if !approx(out2, 15, 0.5) {
		t.Errorf("tip for good = %v want ~15", out2)
	}
}

func TestSugeno(t *testing.T) {
	x := NewVariable("x", 0, 10)
	x.AddTerm("small", Triangular(0, 0, 10))
	x.AddTerm("large", Triangular(0, 10, 10))

	s := NewSugeno()
	s.AddInput(x)
	// z = x for small, z = 2x for large
	s.AddRule(OpAnd, SugenoConsequent{Coeffs: map[string]float64{"x": 1}}, 1, If("x", "small"))
	s.AddRule(OpAnd, SugenoConsequent{Coeffs: map[string]float64{"x": 2}}, 1, If("x", "large"))

	// at x=0: only small fires (w=1), z=0
	out, err := s.Infer(map[string]float64{"x": 0})
	if err != nil {
		t.Fatal(err)
	}
	if !approx(out, 0, tol) {
		t.Errorf("sugeno x=0 = %v want 0", out)
	}
	// at x=5: small w=0.5 z=5, large w=0.5 z=10 -> (0.5*5+0.5*10)/1 = 7.5
	out2, _ := s.Infer(map[string]float64{"x": 5})
	if !approx(out2, 7.5, tol) {
		t.Errorf("sugeno x=5 = %v want 7.5", out2)
	}
}

func TestCartesianProductAndComposeSet(t *testing.T) {
	a, _ := NewSet([]float64{1, 2}, []float64{1.0, 0.5})
	b, _ := NewSet([]float64{1, 2, 3}, []float64{0.2, 0.8, 1.0})
	rel := CartesianProduct(a, b, TNormMin)
	if !approx(rel.M[0][2], 1.0, tol) {
		t.Errorf("cartesian[0][2] = %v want 1.0", rel.M[0][2])
	}
	if !approx(rel.M[1][2], 0.5, tol) {
		t.Errorf("cartesian[1][2] = %v want 0.5", rel.M[1][2])
	}
	// compose a with the relation
	out, err := rel.ComposeSet(a)
	if err != nil {
		t.Fatal(err)
	}
	if len(out.X) != 3 {
		t.Errorf("composed set universe len = %d want 3", len(out.X))
	}
}

func TestTransitiveClosure(t *testing.T) {
	u := []float64{1, 2, 3}
	r, _ := NewRelation(u, u, [][]float64{
		{1.0, 0.7, 0.0},
		{0.7, 1.0, 0.6},
		{0.0, 0.6, 1.0},
	})
	tc, err := r.TransitiveClosure()
	if err != nil {
		t.Fatal(err)
	}
	if !tc.IsMaxMinTransitive(1e-9) {
		t.Error("transitive closure should be max-min transitive")
	}
	// closure (1,3) should be at least min(0.7,0.6)=0.6
	if tc.M[0][2] < 0.6-tol {
		t.Errorf("closure[0][2] = %v want >= 0.6", tc.M[0][2])
	}
}

func TestSubsethood(t *testing.T) {
	xs := []float64{1, 2, 3}
	a, _ := NewSet(xs, []float64{0.2, 0.4, 0.6})
	b, _ := NewSet(xs, []float64{0.5, 0.8, 1.0})
	if !a.IsSubset(b, tol) {
		t.Error("a should be subset of b")
	}
	deg, _ := a.DegreeOfSubsethood(b)
	if !approx(deg, 1, tol) {
		t.Errorf("degree of subsethood = %v want 1", deg)
	}
}

func ExampleMamdaniSystem_Infer() {
	// A one-input Mamdani controller mapping temperature to fan speed.
	temp := NewVariable("temp", 0, 40)
	temp.AddTerm("cold", Trapezoidal(0, 0, 10, 20))
	temp.AddTerm("warm", Triangular(10, 20, 30))
	temp.AddTerm("hot", Trapezoidal(20, 30, 40, 40))

	fan := NewVariable("fan", 0, 100)
	fan.AddTerm("slow", Triangular(0, 0, 50))
	fan.AddTerm("medium", Triangular(0, 50, 100))
	fan.AddTerm("fast", Triangular(50, 100, 100))

	m := NewMamdani(fan)
	m.AddInput(temp)
	m.AddRule(OpAnd, "slow", 1, If("temp", "cold"))
	m.AddRule(OpAnd, "medium", 1, If("temp", "warm"))
	m.AddRule(OpAnd, "fast", 1, If("temp", "hot"))

	speed, _ := m.Infer(map[string]float64{"temp": 20})
	fmt.Printf("%.1f\n", speed)
	// Output: 50.0
}

func ExampleSet_Centroid() {
	s := FromMF(Triangular(0, 5, 10), Linspace(0, 10, 11))
	c, _ := s.Centroid()
	fmt.Printf("%.1f\n", c)
	// Output: 5.0
}
