package dynamical

import (
	"math"
	"math/cmplx"
	"testing"
)

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func TestLogisticBasics(t *testing.T) {
	if got := Logistic(4, 0.5); !approx(got, 1, 1e-12) {
		t.Errorf("Logistic(4,0.5)=%v want 1", got)
	}
	if got := LogisticDeriv(4, 0.5); !approx(got, 0, 1e-12) {
		t.Errorf("LogisticDeriv(4,0.5)=%v want 0", got)
	}
	f := LogisticMap(3.2)
	if got := f(0.5); !approx(got, 0.8, 1e-12) {
		t.Errorf("LogisticMap(3.2)(0.5)=%v want 0.8", got)
	}
}

func TestLogisticFixedPoints(t *testing.T) {
	r := 2.5
	pts := LogisticFixedPoints(r)
	want := []float64{0, 1 - 1/r}
	if len(pts) != 2 {
		t.Fatalf("got %d fixed points want 2", len(pts))
	}
	for i := range want {
		if !approx(pts[i], want[i], 1e-12) {
			t.Errorf("fixed point %d = %v want %v", i, pts[i], want[i])
		}
		if !IsFixedPoint(LogisticMap(r), pts[i], 1e-12) {
			t.Errorf("point %v not fixed", pts[i])
		}
	}
	// Multiplier at 0 is r; at 1-1/r is 2-r. For r=2.5: 2.5 (unstable), -0.5 (stable).
	st := LogisticStability(r, 1e-9)
	if st[0] != Unstable {
		t.Errorf("origin stability = %v want unstable", st[0])
	}
	if st[1] != Stable {
		t.Errorf("nonzero fixed point stability = %v want stable", st[1])
	}
}

func TestClassifyStability(t *testing.T) {
	if ClassifyStability(0, 1e-9) != SuperStable {
		t.Error("zero multiplier should be superstable")
	}
	if ClassifyStability(1, 1e-9) != Neutral {
		t.Error("unit multiplier should be neutral")
	}
	if ClassifyStability(0.5, 1e-9) != Stable {
		t.Error("|m|<1 should be stable")
	}
	if ClassifyStability(2, 1e-9) != Unstable {
		t.Error("|m|>1 should be unstable")
	}
	if Stable.String() != "stable" {
		t.Errorf("Stable.String()=%q", Stable.String())
	}
}

func TestTentFixedPoints(t *testing.T) {
	mu := 1.9
	pts := TentFixedPoints(mu)
	want := mu / (1 + mu)
	if len(pts) != 2 || !approx(pts[1], want, 1e-12) {
		t.Fatalf("tent fixed points %v want [0 %v]", pts, want)
	}
	if !IsFixedPoint(TentMap(mu), pts[1], 1e-12) {
		t.Errorf("tent point %v not fixed", pts[1])
	}
}

func TestFixedPoint1DNewton(t *testing.T) {
	r := 3.2
	f := LogisticMap(r)
	df := func(x float64) float64 { return LogisticDeriv(r, x) }
	root, _, ok := FixedPoint1D(f, df, 0.9, 100, 1e-12)
	if !ok || !approx(root, 1-1/r, 1e-9) {
		t.Errorf("FixedPoint1D=%v ok=%v want %v", root, ok, 1-1/r)
	}
}

func TestLyapunovTentExact(t *testing.T) {
	mu := 1.9
	f := TentMap(mu)
	df := func(x float64) float64 {
		if x < 0.5 {
			return mu
		}
		return -mu
	}
	// The tent map has constant |slope|=mu, so the Lyapunov exponent is
	// exactly log(mu) for any orbit.
	le := Lyapunov(f, df, 0.1234, 100, 5000)
	if !approx(le, math.Log(mu), 1e-9) {
		t.Errorf("tent Lyapunov=%v want %v", le, math.Log(mu))
	}
}

func TestLyapunovLogisticR4(t *testing.T) {
	// Known exact ergodic value for the logistic map at r=4 is log(2).
	f := LogisticMap(4)
	df := func(x float64) float64 { return LogisticDeriv(4, x) }
	le := Lyapunov(f, df, 0.1, 1000, 400000)
	if !approx(le, math.Log(2), 0.02) {
		t.Errorf("logistic r=4 Lyapunov=%v want ~%v", le, math.Log(2))
	}
}

func TestLyapunovSeparationLogistic(t *testing.T) {
	f := LogisticMap(4)
	le := LyapunovSeparation(f, 0.1, 1e-9, 1000, 100000)
	if !approx(le, math.Log(2), 0.05) {
		t.Errorf("separation Lyapunov=%v want ~%v", le, math.Log(2))
	}
}

func TestLyapunovHenon(t *testing.T) {
	le := LyapunovHenon(1.4, 0.3, 1000, 100000)
	if !approx(le, 0.419, 0.03) {
		t.Errorf("Henon largest Lyapunov=%v want ~0.419", le)
	}
}

func TestHenonFixedPoints(t *testing.T) {
	a, b := 1.4, 0.3
	pts := HenonFixedPoints(a, b)
	if len(pts) != 2 {
		t.Fatalf("got %d Henon fixed points want 2", len(pts))
	}
	f := HenonMap(a, b)
	for _, p := range pts {
		q := f(p)
		if q.Sub(p).Norm() > 1e-12 {
			t.Errorf("Henon point %v maps to %v", p, q)
		}
	}
}

func TestOrbitAndIterate(t *testing.T) {
	f := LogisticMap(2)
	orb := Orbit(f, 0.1, 5)
	if len(orb) != 6 {
		t.Fatalf("orbit length %d want 6", len(orb))
	}
	if !approx(orb[0], 0.1, 1e-15) {
		t.Errorf("orbit[0]=%v", orb[0])
	}
	if !approx(NthIterate(f, 0.1, 5), orb[5], 1e-15) {
		t.Errorf("NthIterate disagrees with Orbit")
	}
	// Logistic r=2 converges to fixed point 0.5.
	tail := OrbitTransient(f, 0.1, 100, 0)
	if !approx(tail[0], 0.5, 1e-9) {
		t.Errorf("r=2 attractor=%v want 0.5", tail[0])
	}
}

func TestPeriodDetection(t *testing.T) {
	// Rigid rotation by 1/4 has period 4.
	rot := CircleRotation(0.25)
	if p := DetectPeriod(rot, 0, 0, 12, 1e-9); p != 4 {
		t.Errorf("rotation period=%d want 4", p)
	}
	if !IsPeriodicPoint(rot, 0, 4, 1e-9) {
		t.Error("0 should be period-4 point of rotation")
	}
	if IsPeriodicPoint(rot, 0, 2, 1e-9) {
		t.Error("0 should NOT be period-2 point of rotation by 1/4")
	}
	// Doubling map: 1/3 -> 2/3 -> 1/3 has period 2.
	if p := DetectPeriod(DoublingMap(), 1.0/3.0, 0, 8, 1e-9); p != 2 {
		t.Errorf("doubling period=%d want 2", p)
	}
}

func TestCobweb(t *testing.T) {
	f := LogisticMap(2)
	segs := Cobweb(f, 0.1, 3)
	if len(segs) != 6 {
		t.Fatalf("cobweb segments=%d want 6", len(segs))
	}
	// First vertical segment goes from (x0,x0) to (x0,f(x0)).
	s := segs[0]
	if !approx(s.X0, 0.1, 1e-15) || !approx(s.Y0, 0.1, 1e-15) ||
		!approx(s.X1, 0.1, 1e-15) || !approx(s.Y1, f(0.1), 1e-15) {
		t.Errorf("first cobweb segment = %+v", s)
	}
}

func TestBifurcationConverges(t *testing.T) {
	// Below the first bifurcation the logistic map settles on 1-1/r.
	pts := LogisticBifurcation(2.5, 2.5, 1, 0.2, 2000, 10)
	if len(pts) != 1 {
		t.Fatalf("got %d bifurcation points", len(pts))
	}
	want := 1 - 1/2.5
	for _, v := range pts[0].Values {
		if !approx(v, want, 1e-6) {
			t.Errorf("bifurcation sample=%v want %v", v, want)
		}
	}
}

func TestNewtonComplexCubeRoots(t *testing.T) {
	f := func(z complex128) complex128 { return z*z*z - 1 }
	df := func(z complex128) complex128 { return 3 * z * z }
	res := NewtonComplex(f, df, complex(1, 0.1), 100, 1e-12)
	if !res.Converged || cmplx.Abs(res.Root-1) > 1e-9 {
		t.Errorf("Newton cube root=%v converged=%v want 1", res.Root, res.Converged)
	}
	// Start near the primitive root e^{2pi i/3}.
	w := cmplx.Exp(complex(0, 2*math.Pi/3))
	res2 := NewtonComplex(f, df, w*complex(1.01, 0), 100, 1e-12)
	if cmplx.Abs(res2.Root-w) > 1e-9 {
		t.Errorf("Newton root=%v want %v", res2.Root, w)
	}
}

func TestNewtonReal(t *testing.T) {
	f := func(x float64) float64 { return x*x - 2 }
	df := func(x float64) float64 { return 2 * x }
	root, _, ok := NewtonReal(f, df, 1, 100, 1e-14)
	if !ok || !approx(root, math.Sqrt2, 1e-12) {
		t.Errorf("Newton sqrt2=%v ok=%v", root, ok)
	}
}

func TestNewtonBasin(t *testing.T) {
	f := func(z complex128) complex128 { return z*z*z - 1 }
	df := func(z complex128) complex128 { return 3 * z * z }
	roots := []complex128{
		complex(1, 0),
		cmplx.Exp(complex(0, 2*math.Pi/3)),
		cmplx.Exp(complex(0, -2*math.Pi/3)),
	}
	grid := NewtonBasin(f, df, roots, -1, 1, -1, 1, 11, 11, 100, 1e-12)
	if len(grid) != 11 || len(grid[0]) != 11 {
		t.Fatalf("grid shape %dx%d", len(grid), len(grid[0]))
	}
	// A point near the real root 1 must land in basin 0.
	if grid[5][10] != 0 {
		t.Errorf("point near root 1 classified as %d want 0", grid[5][10])
	}
}

func TestRK4Decay(t *testing.T) {
	// dx/dt = -x has solution x(t)=x0*exp(-t); RK4 should be very accurate.
	field := func(s Vec3) Vec3 { return s.Scale(-1) }
	s0 := Vec3{1, 2, 3}
	traj := Integrate3D(field, s0, 1e-3, 1000) // integrate to t=1
	last := traj[len(traj)-1]
	e := math.Exp(-1)
	if !approx(last.X, s0.X*e, 1e-8) || !approx(last.Y, s0.Y*e, 1e-8) || !approx(last.Z, s0.Z*e, 1e-8) {
		t.Errorf("RK4 decay end=%+v want %v*exp(-1)", last, s0)
	}
}

func TestLorenzFixedPoints(t *testing.T) {
	p := DefaultLorenz()
	fps := p.LorenzFixedPoints()
	if len(fps) != 3 {
		t.Fatalf("got %d Lorenz fixed points want 3", len(fps))
	}
	field := p.Field()
	for _, fp := range fps {
		if field(fp).Norm() > 1e-9 {
			t.Errorf("Lorenz fixed point %v has nonzero field %v", fp, field(fp))
		}
	}
	c := math.Sqrt(p.Beta * (p.Rho - 1))
	if !approx(fps[1].Z, p.Rho-1, 1e-12) || !approx(fps[1].X, c, 1e-12) {
		t.Errorf("C+ = %v", fps[1])
	}
}

func TestLorenzTrajectoryBounded(t *testing.T) {
	p := DefaultLorenz()
	traj := LorenzTrajectory(p, Vec3{1, 1, 1}, 0.01, 5000)
	last := traj[len(traj)-1]
	if math.IsNaN(last.Norm()) || last.Norm() > 200 {
		t.Errorf("Lorenz trajectory escaped: %v", last)
	}
}

func TestRosslerTrajectoryBounded(t *testing.T) {
	p := DefaultRossler()
	traj := RosslerTrajectory(p, Vec3{1, 1, 1}, 0.02, 5000)
	last := traj[len(traj)-1]
	if math.IsNaN(last.Norm()) || last.Norm() > 500 {
		t.Errorf("Rossler trajectory escaped: %v", last)
	}
}

func TestVec3Ops(t *testing.T) {
	a := Vec3{1, 2, 2}
	if !approx(a.Norm(), 3, 1e-12) {
		t.Errorf("norm=%v want 3", a.Norm())
	}
	if a.Add(Vec3{1, 1, 1}) != (Vec3{2, 3, 3}) {
		t.Error("Add wrong")
	}
	if a.Sub(Vec3{1, 2, 2}) != (Vec3{0, 0, 0}) {
		t.Error("Sub wrong")
	}
	if a.Dot(Vec3{1, 0, 0}) != 1 {
		t.Error("Dot wrong")
	}
}

func BenchmarkNewtonBasin(b *testing.B) {
	f := func(z complex128) complex128 { return z*z*z - 1 }
	df := func(z complex128) complex128 { return 3 * z * z }
	roots := []complex128{
		complex(1, 0),
		cmplx.Exp(complex(0, 2*math.Pi/3)),
		cmplx.Exp(complex(0, -2*math.Pi/3)),
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		NewtonBasin(f, df, roots, -2, 2, -2, 2, 64, 64, 50, 1e-10)
	}
}
