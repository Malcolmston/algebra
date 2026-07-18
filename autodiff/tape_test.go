package autodiff

import (
	"math"
	"testing"
)

func TestReverseGradient(t *testing.T) {
	// f(x,y) = x*y + sin(x). df/dx = y + cos(x); df/dy = x. At (1,2).
	f := func(t *Tape, v []Var) Var {
		return v[0].Mul(v[1]).Add(v[0].Sin())
	}
	g := GradientReverse(f, []float64{1, 2})
	want := []float64{2 + math.Cos(1), 1}
	for i := range want {
		if !closeTo(g[i], want[i], tol) {
			t.Fatalf("grad[%d]=%v want %v", i, g[i], want[i])
		}
	}
}

func TestReverseVsForward(t *testing.T) {
	// Cross-check reverse mode against forward-mode Gradient on a nontrivial
	// function f(x) = exp(x0*x1) + log(x2) + x0²·x3 - tanh(x1).
	fwd := func(v []Dual) Dual {
		return Exp(v[0].Mul(v[1])).
			Add(Log(v[2])).
			Add(v[0].Mul(v[0]).Mul(v[3])).
			Sub(Tanh(v[1]))
	}
	rev := func(t *Tape, v []Var) Var {
		return v[0].Mul(v[1]).Exp().
			Add(v[2].Log()).
			Add(v[0].Mul(v[0]).Mul(v[3])).
			Sub(v[1].Tanh())
	}
	x := []float64{0.4, 1.1, 2.5, -0.7}
	gf := Gradient(fwd, x)
	gr := GradientReverse(rev, x)
	for i := range gf {
		if !closeTo(gf[i], gr[i], 1e-9) {
			t.Fatalf("mismatch at %d: forward %v reverse %v", i, gf[i], gr[i])
		}
	}
}

func TestReverseUnaryOps(t *testing.T) {
	// f(x) = sqrt(x) with x=9: value 3, derivative 1/(2*3)=1/6.
	tp := NewTape()
	x := tp.Variable(9)
	y := x.Sqrt()
	if !closeTo(y.Value(), 3, tol) {
		t.Fatalf("sqrt value %v", y.Value())
	}
	g := tp.Backward(y)
	if !closeTo(g[0], 1.0/6.0, tol) {
		t.Fatalf("sqrt grad %v want %v", g[0], 1.0/6.0)
	}

	// Chained: g(x) = sigmoid(3x) at x=0.5. d/dx = 3*s*(1-s), s=sigmoid(1.5).
	tp2 := NewTape()
	a := tp2.Variable(0.5)
	out := a.Scale(3).Sigmoid()
	s := 1 / (1 + math.Exp(-1.5))
	if !closeTo(out.Value(), s, tol) {
		t.Fatalf("sigmoid value %v want %v", out.Value(), s)
	}
	grad := tp2.Backward(out)
	if !closeTo(grad[0], 3*s*(1-s), tol) {
		t.Fatalf("sigmoid grad %v want %v", grad[0], 3*s*(1-s))
	}
}

func TestTapeNumVars(t *testing.T) {
	tp := NewTape()
	tp.Variable(1)
	tp.Constant(2)
	tp.Variable(3)
	if tp.NumVars() != 2 {
		t.Fatalf("NumVars=%d want 2", tp.NumVars())
	}
}

// BenchmarkHessian exercises the heaviest routine: a full n×n Hessian assembly
// via hyper-dual numbers on a moderately sized quadratic-plus-transcendental
// objective.
func BenchmarkHessian(b *testing.B) {
	const n = 12
	f := func(v []HyperDual) HyperDual {
		acc := HyperConstant(0)
		for i := 0; i < n; i++ {
			acc = acc.Add(v[i].Mul(v[i]).Scale(0.5))
			if i+1 < n {
				acc = acc.Add(v[i].Mul(v[i+1]).Sin())
			}
		}
		return acc
	}
	x := make([]float64, n)
	for i := range x {
		x[i] = 0.1 * float64(i+1)
	}
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = Hessian(f, x)
	}
}
