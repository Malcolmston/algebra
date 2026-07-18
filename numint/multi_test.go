package numint

import (
	"math"
	"math/rand"
	"testing"
)

func numintmvClose(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.IsNaN(got) || math.Abs(got-want) > tol {
		t.Errorf("%s = %v, want %v (tol %v)", name, got, want, tol)
	}
}

func TestRectangleDoubleIntegrals(t *testing.T) {
	// integral of x*y over [0,1]x[0,1] = 1/4
	xy := func(x, y float64) float64 { return x * y }
	numintmvClose(t, "DoubleSimpson", DoubleSimpson(xy, 0, 1, 0, 1, 8, 8), 0.25, 1e-12)
	numintmvClose(t, "DoubleBoole", DoubleBoole(xy, 0, 1, 0, 1, 8, 8), 0.25, 1e-12)
	numintmvClose(t, "DoubleGaussLegendre", DoubleGaussLegendre(xy, 0, 1, 0, 1, 4), 0.25, 1e-12)
	numintmvClose(t, "DoubleTrapezoid", DoubleTrapezoid(xy, 0, 1, 0, 1, 200, 200), 0.25, 1e-4)
	numintmvClose(t, "DoubleMidpoint", DoubleMidpoint(xy, 0, 1, 0, 1, 200, 200), 0.25, 1e-4)
	numintmvClose(t, "DoubleIntegral", DoubleIntegral(xy, 0, 1, 0, 1), 0.25, 1e-9)

	// integral of x^2*y^2 over [0,1]^2 = 1/9
	x2y2 := func(x, y float64) float64 { return x * x * y * y }
	numintmvClose(t, "DoubleGL x2y2", DoubleGaussLegendre(x2y2, 0, 1, 0, 1, 5), 1.0/9, 1e-12)

	// mean of x+y over [0,2]x[0,2] = 2
	sum := func(x, y float64) float64 { return x + y }
	numintmvClose(t, "DoubleAverage", DoubleAverage(sum, 0, 2, 0, 2), 2.0, 1e-9)
}

func TestRectangleTripleIntegrals(t *testing.T) {
	unit := func(x, y, z float64) float64 { return 1 }
	numintmvClose(t, "TripleSimpson vol", TripleSimpson(unit, 0, 1, 0, 1, 0, 1, 4, 4, 4), 1.0, 1e-12)
	numintmvClose(t, "TripleTrapezoid vol", TripleTrapezoid(unit, 0, 1, 0, 1, 0, 1, 10, 10, 10), 1.0, 1e-12)

	prod := func(x, y, z float64) float64 { return x * y * z }
	numintmvClose(t, "TripleGaussLegendre", TripleGaussLegendre(prod, 0, 1, 0, 1, 0, 1, 4), 1.0/8, 1e-12)
	numintmvClose(t, "TripleSimpson prod", TripleSimpson(prod, 0, 1, 0, 1, 0, 1, 6, 6, 6), 1.0/8, 1e-12)
	numintmvClose(t, "TripleMidpoint prod", TripleMidpoint(prod, 0, 1, 0, 1, 0, 1, 100, 100, 100), 1.0/8, 1e-4)
	numintmvClose(t, "TripleIntegral prod", TripleIntegral(prod, 0, 1, 0, 1, 0, 1), 1.0/8, 1e-9)
	numintmvClose(t, "TripleAverage", TripleAverage(prod, 0, 1, 0, 1, 0, 1), 1.0/8, 1e-9)
}

func TestGeneralRegion(t *testing.T) {
	one := func(x, y float64) float64 { return 1 }
	zero := func(x float64) float64 { return 0 }
	line := func(x float64) float64 { return x }
	// area of triangle 0<=y<=x, 0<=x<=1 is 1/2
	numintmvClose(t, "RegionSimpson", DoubleIntegralRegion(one, 0, 1, zero, line, 16, 16), 0.5, 1e-10)
	numintmvClose(t, "RegionGL", DoubleIntegralRegionGL(one, 0, 1, zero, line, 8), 0.5, 1e-10)

	// area between y=x and y=x^2 over [0,1] = 1/6
	sq := func(x float64) float64 { return x * x }
	numintmvClose(t, "AreaBetween", AreaBetween(sq, line, 0, 1, 20), 1.0/6, 1e-10)

	// centroid of rectangle 0<=y<=1, 0<=x<=2 is (1, 0.5)
	oneF := func(x float64) float64 { return 1 }
	xb, yb := CentroidRegion(zero, oneF, 0, 2, 20)
	numintmvClose(t, "CentroidX", xb, 1.0, 1e-10)
	numintmvClose(t, "CentroidY", yb, 0.5, 1e-10)
}

func TestPlanarSpecialRegions(t *testing.T) {
	one := func(x, y float64) float64 { return 1 }
	// triangle area
	numintmvClose(t, "TriangleArea", IntegrateTriangle(one, 0, 0, 1, 0, 0, 1, 8), 0.5, 1e-10)
	// integral of x over unit right triangle = 1/6
	fx := func(x, y float64) float64 { return x }
	numintmvClose(t, "TriangleX", IntegrateTriangle(fx, 0, 0, 1, 0, 0, 1, 12), 1.0/6, 1e-10)

	// disk area of radius 2 = 4*pi
	numintmvClose(t, "DiskArea", IntegrateDisk(one, 0, 0, 2, 40, 40), 4*math.Pi, 1e-6)
	// integral of x^2+y^2 over unit disk = pi/2
	r2 := func(x, y float64) float64 { return x*x + y*y }
	numintmvClose(t, "DiskR2", IntegrateDisk(r2, 0, 0, 1, 40, 40), math.Pi/2, 1e-8)

	// polar: integral of r over full unit disk = pi
	gpolar := func(r, th float64) float64 { return 1 }
	numintmvClose(t, "Polar", IntegratePolar(gpolar, 0, 1, 0, 2*math.Pi, 20, 40), math.Pi, 1e-8)

	// annulus area between r=1 and r=2 = 3*pi
	numintmvClose(t, "Annulus", IntegrateAnnulus(one, 0, 0, 1, 2, 40, 40), 3*math.Pi, 1e-6)
}

func TestVolumeIntegrals(t *testing.T) {
	one := func(x, y, z float64) float64 { return 1 }
	numintmvClose(t, "SphereSurface", IntegrateSphereSurface(one, 0, 0, 0, 1, 60, 60), 4*math.Pi, 1e-4)
	numintmvClose(t, "BallVolume", IntegrateBall(one, 0, 0, 0, 1, 40, 40, 40), 4.0/3*math.Pi, 1e-4)
	// integral of x^2+y^2+z^2 over unit ball = 4*pi/5
	r2 := func(x, y, z float64) float64 { return x*x + y*y + z*z }
	numintmvClose(t, "BallR2", IntegrateBall(r2, 0, 0, 0, 1, 60, 60, 60), 4*math.Pi/5, 1e-4)
	// cylinder volume radius 1, height 2 = 2*pi
	numintmvClose(t, "CylinderVol", IntegrateCylinder(one, 0, 0, 1, 0, 2, 40, 40, 4), 2*math.Pi, 1e-6)
}

func TestGaussHermite(t *testing.T) {
	// integral of exp(-x^2) over R = sqrt(pi)
	one := func(x float64) float64 { return 1 }
	numintmvClose(t, "GH const", GaussHermite(one, 16), math.Sqrt(math.Pi), 1e-12)
	// integral of exp(-x^2) x^2 = sqrt(pi)/2
	x2 := func(x float64) float64 { return x * x }
	numintmvClose(t, "GH x2", GaussHermite(x2, 16), math.Sqrt(math.Pi)/2, 1e-12)
	// integral of exp(-x^2) x^4 = 3/4 sqrt(pi)
	x4 := func(x float64) float64 { return x * x * x * x }
	numintmvClose(t, "GH x4", GaussHermite(x4, 16), 0.75*math.Sqrt(math.Pi), 1e-11)

	nodes, weights := GaussHermiteNodes(10)
	ws := 0.0
	for _, w := range weights {
		ws += w
	}
	numintmvClose(t, "GH weight sum", ws, math.Sqrt(math.Pi), 1e-12)
	if len(nodes) != 10 {
		t.Fatalf("expected 10 nodes, got %d", len(nodes))
	}

	// E[X^2] for standard normal is 1; E[X] for N(3,2) is 3
	numintmvClose(t, "ExpNormal x2", ExpectationNormal(x2, 0, 1, 16), 1.0, 1e-12)
	id := func(x float64) float64 { return x }
	numintmvClose(t, "ExpNormal mean", ExpectationNormal(id, 3, 2, 16), 3.0, 1e-11)

	// integral of exp(-x^2/2) over R = sqrt(2*pi)
	numintmvClose(t, "GaussianWeighted", IntegrateGaussianWeighted(one, 16), math.Sqrt(2*math.Pi), 1e-12)
}

func TestGaussLaguerre(t *testing.T) {
	one := func(x float64) float64 { return 1 }
	id := func(x float64) float64 { return x }
	x2 := func(x float64) float64 { return x * x }
	x3 := func(x float64) float64 { return x * x * x }
	// moments of exp(-x): int x^k exp(-x) = k!
	numintmvClose(t, "GL 0!", GaussLaguerre(one, 12), 1.0, 1e-11)
	numintmvClose(t, "GL 1!", GaussLaguerre(id, 12), 1.0, 1e-11)
	numintmvClose(t, "GL 2!", GaussLaguerre(x2, 12), 2.0, 1e-10)
	numintmvClose(t, "GL 3!", GaussLaguerre(x3, 12), 6.0, 1e-10)

	// generalized: int x^alpha exp(-x) = Gamma(alpha+1)
	numintmvClose(t, "GLgen a=1", GaussLaguerreGen(one, 1, 12), math.Gamma(2), 1e-11)
	numintmvClose(t, "GLgen a=2", GaussLaguerreGen(one, 2, 12), math.Gamma(3), 1e-10)

	nodes, weights := GaussLaguerreNodes(8)
	ws := 0.0
	for _, w := range weights {
		ws += w
	}
	numintmvClose(t, "GL weight sum", ws, 1.0, 1e-12)
	if len(nodes) != 8 {
		t.Fatalf("expected 8 nodes, got %d", len(nodes))
	}

	// int_0^inf exp(-2x) dx = 1/2
	numintmvClose(t, "ExpWeighted", IntegrateExpWeighted(one, 2, 12), 0.5, 1e-11)
	// Laplace transform of f(t)=t at s=2 is 1/s^2 = 1/4
	numintmvClose(t, "Laplace", LaplaceTransform(id, 2, 12), 0.25, 1e-11)
}

func TestImproperIntegrals(t *testing.T) {
	// int_{-inf}^{inf} 1/(1+x^2) dx = pi
	cauchy := func(x float64) float64 { return 1 / (1 + x*x) }
	numintmvClose(t, "Infinite cauchy", IntegrateInfinite(cauchy, 64), math.Pi, 1e-9)
	// int_{-inf}^{inf} exp(-x^2) dx = sqrt(pi)
	gauss := func(x float64) float64 { return math.Exp(-x * x) }
	numintmvClose(t, "Infinite gauss", IntegrateInfinite(gauss, 120), math.Sqrt(math.Pi), 1e-4)

	// int_0^inf exp(-x) dx = 1
	numintmvClose(t, "UpperInfinite exp", IntegrateUpperInfinite(func(x float64) float64 { return math.Exp(-x) }, 0, 80), 1.0, 1e-6)
	// int_1^inf 1/x^2 dx = 1
	numintmvClose(t, "UpperInfinite inv", IntegrateUpperInfinite(func(x float64) float64 { return 1 / (x * x) }, 1, 32), 1.0, 1e-9)
	// int_{-inf}^0 exp(x) dx = 1
	numintmvClose(t, "LowerInfinite", IntegrateLowerInfinite(func(x float64) float64 { return math.Exp(x) }, 0, 80), 1.0, 1e-6)

	// int_0^1 1/sqrt(x) dx = 2
	numintmvClose(t, "Singular", IntegrateSingularEndpoint(func(x float64) float64 { return 1 / math.Sqrt(x) }, 0, 1, 20), 2.0, 1e-10)
}

func TestLineIntegrals(t *testing.T) {
	circle := func(t float64) (float64, float64) { return math.Cos(t), math.Sin(t) }
	// arc length of unit circle = 2*pi
	numintmvClose(t, "ArcLength2D", ArcLength2D(circle, 0, 2*math.Pi, 40), 2*math.Pi, 1e-6)
	// scalar line integral of 1 over unit circle = length = 2*pi
	one := func(x, y float64) float64 { return 1 }
	numintmvClose(t, "LineScalar2D", LineIntegralScalar2D(one, circle, 0, 2*math.Pi, 40), 2*math.Pi, 1e-6)

	// helix arc length over [0,2pi] = 2*pi*sqrt(2)
	helix := func(t float64) (float64, float64, float64) { return math.Cos(t), math.Sin(t), t }
	numintmvClose(t, "ArcLength3D", ArcLength3D(helix, 0, 2*math.Pi, 40), 2*math.Pi*math.Sqrt2, 1e-5)

	// circulation of (-y, x) around unit circle = 2*pi
	rot := func(x, y float64) (float64, float64) { return -y, x }
	numintmvClose(t, "Circulation2D", Circulation2D(rot, circle, 0, 2*math.Pi, 60), 2*math.Pi, 1e-5)
	numintmvClose(t, "LineVector2D", LineIntegralVector2D(rot, circle, 0, 2*math.Pi, 60), 2*math.Pi, 1e-5)

	// flux of (x, y) across unit circle = 2*pi
	radial := func(x, y float64) (float64, float64) { return x, y }
	numintmvClose(t, "Flux2D", Flux2D(radial, circle, 0, 2*math.Pi, 60), 2*math.Pi, 1e-5)

	// vertical field integrated along z-axis segment = 1
	up := func(x, y, z float64) (float64, float64, float64) { return 0, 0, 1 }
	seg := func(t float64) (float64, float64, float64) { return 0, 0, t }
	numintmvClose(t, "LineVector3D", LineIntegralVector3D(up, seg, 0, 1, 8), 1.0, 1e-9)
	numintmvClose(t, "LineScalar3D", LineIntegralScalar3D(func(x, y, z float64) float64 { return 1 }, seg, 0, 1, 8), 1.0, 1e-9)

	// polyline (0,0)->(1,0)->(1,1) has length 2
	xs := []float64{0, 1, 1}
	ys := []float64{0, 0, 1}
	numintmvClose(t, "PolylineLen", ArcLengthPolyline2D(xs, ys), 2.0, 1e-12)
	numintmvClose(t, "PolylineScalar", LineIntegralPolyline2D(one, xs, ys), 2.0, 1e-12)
}

func TestMonteCarlo(t *testing.T) {
	// int_0^1 x dx = 0.5 with error estimate
	res := MonteCarloWithError(func(x float64) float64 { return x }, 0, 1, 200000, 1)
	numintmvClose(t, "MC1D value", res.Value, 0.5, 5e-3)
	if res.Evals != 200000 || res.ErrEst <= 0 {
		t.Errorf("unexpected QuadResult %+v", res)
	}

	// int over [0,1]^2 of x+y = 1
	res2 := MonteCarlo2DWithError(func(x, y float64) float64 { return x + y }, 0, 1, 0, 1, 200000, 2)
	numintmvClose(t, "MC2D value", res2.Value, 1.0, 5e-3)

	// triple integral of constant 1 over [0,2]^3 = 8 (exact regardless of seed)
	numintmvClose(t, "MC3D", MonteCarlo3D(func(x, y, z float64) float64 { return 1 }, 0, 2, 0, 2, 0, 2, 1000, 5), 8.0, 1e-9)

	// antithetic is exact for linear integrands
	numintmvClose(t, "MCantithetic", MonteCarloAntithetic(func(x float64) float64 { return x }, 0, 1, 1000, 9), 0.5, 1e-12)

	// stratified estimate of int_0^1 x^2 dx = 1/3
	numintmvClose(t, "MCstratified", MonteCarloStratified(func(x float64) float64 { return x * x }, 0, 1, 200, 20, 11), 1.0/3, 5e-3)

	// importance sampling for int_0^1 x^2 dx = 1/3 with uniform proposal
	imp := MonteCarloImportance(
		func(x float64) float64 { return x * x },
		func(r *rand.Rand) float64 { return r.Float64() },
		func(x float64) float64 { return 1 },
		300000, 7)
	numintmvClose(t, "MCimportance", imp, 1.0/3, 5e-3)

	// disk area radius 2 exact for constant integrand
	numintmvClose(t, "MCdisk", MonteCarloDisk(func(x, y float64) float64 { return 1 }, 0, 0, 2, 1000, 3), 4*math.Pi, 1e-9)
	// ball volume radius 1 exact for constant integrand
	numintmvClose(t, "MCball", MonteCarloBall(func(x, y, z float64) float64 { return 1 }, 0, 0, 0, 1, 1000, 4), 4.0/3*math.Pi, 1e-9)

	// volume of unit disk inside [-1,1]^2 = pi
	vol := MonteCarloVolume(func(p []float64) bool { return p[0]*p[0]+p[1]*p[1] <= 1 },
		[]float64{-1, -1}, []float64{1, 1}, 400000, 3)
	numintmvClose(t, "MCvolume", vol, math.Pi, 1e-2)

	// integral of 1 over triangle region 0<=y<=x, 0<=x<=1 = 1/2
	reg := MonteCarloRegion(func(x, y float64) float64 { return 1 }, 0, 1,
		func(x float64) float64 { return 0 }, func(x float64) float64 { return x }, 300000, 6)
	numintmvClose(t, "MCregion", reg, 0.5, 5e-3)
}

func TestQuasiMonteCarloAndHalton(t *testing.T) {
	// van der Corput base-2 known values
	numintmvClose(t, "vdc1", VanDerCorput(1, 2), 0.5, 0)
	numintmvClose(t, "vdc2", VanDerCorput(2, 2), 0.25, 0)
	numintmvClose(t, "vdc3", VanDerCorput(3, 2), 0.75, 0)

	pt := HaltonPoint(1, []int{2, 3})
	numintmvClose(t, "halton x", pt[0], 0.5, 0)
	numintmvClose(t, "halton y", pt[1], 1.0/3, 1e-15)

	seq := HaltonSequence(4, []int{2})
	if len(seq) != 4 || seq[0][0] != 0.5 {
		t.Errorf("unexpected halton sequence %v", seq)
	}

	primes := HaltonPrimes(4)
	want := []int{2, 3, 5, 7}
	for i := range want {
		if primes[i] != want[i] {
			t.Errorf("HaltonPrimes[%d] = %d, want %d", i, primes[i], want[i])
		}
	}

	// quasi Monte Carlo integrals
	numintmvClose(t, "QMC x", QuasiMonteCarlo(func(x float64) float64 { return x }, 0, 1, 4096), 0.5, 1e-3)
	numintmvClose(t, "QMC2D xy", QuasiMonteCarlo2D(func(x, y float64) float64 { return x * y }, 0, 1, 0, 1, 8192), 0.25, 2e-3)
	numintmvClose(t, "QMCnd xyz", QuasiMonteCarloND(func(p []float64) float64 { return p[0] * p[1] * p[2] },
		[]float64{0, 0, 0}, []float64{1, 1, 1}, 16384), 0.125, 2e-3)
}

func BenchmarkIntegrateBall(b *testing.B) {
	f := func(x, y, z float64) float64 { return math.Exp(-(x*x + y*y + z*z)) }
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		IntegrateBall(f, 0, 0, 0, 1, 40, 40, 40)
	}
}
