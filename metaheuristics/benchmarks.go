package metaheuristics

import "math"

// Sphere is the sum-of-squares benchmark f(x) = sum x_i^2. It is convex,
// unimodal and separable with global minimum 0 at the origin.
func Sphere(x []float64) float64 {
	s := 0.0
	for _, v := range x {
		s += v * v
	}
	return s
}

// Rastrigin is a highly multimodal benchmark
// f(x) = 10n + sum(x_i^2 - 10 cos(2*pi*x_i)). Its global minimum is 0 at the
// origin; it is typically evaluated on [-5.12, 5.12]^n.
func Rastrigin(x []float64) float64 {
	s := 10.0 * float64(len(x))
	for _, v := range x {
		s += v*v - 10*math.Cos(2*math.Pi*v)
	}
	return s
}

// Ackley is the multimodal Ackley benchmark with a nearly flat outer region
// and a deep central basin. Its global minimum is 0 at the origin; it is
// typically evaluated on [-32.768, 32.768]^n.
func Ackley(x []float64) float64 {
	n := float64(len(x))
	if n == 0 {
		return 0
	}
	const (
		a = 20.0
		b = 0.2
		c = 2 * math.Pi
	)
	var sq, cs float64
	for _, v := range x {
		sq += v * v
		cs += math.Cos(c * v)
	}
	return -a*math.Exp(-b*math.Sqrt(sq/n)) - math.Exp(cs/n) + a + math.E
}

// Rosenbrock is the classic banana-valley benchmark
// f(x) = sum_{i<n-1} [100 (x_{i+1}-x_i^2)^2 + (1-x_i)^2]. Its global minimum
// is 0 at the all-ones vector.
func Rosenbrock(x []float64) float64 {
	s := 0.0
	for i := 0; i+1 < len(x); i++ {
		a := x[i+1] - x[i]*x[i]
		b := 1 - x[i]
		s += 100*a*a + b*b
	}
	return s
}

// Griewank is a multimodal benchmark
// f(x) = 1 + sum x_i^2/4000 - prod cos(x_i/sqrt(i+1)). Its global minimum is 0
// at the origin; it is typically evaluated on [-600, 600]^n.
func Griewank(x []float64) float64 {
	sum := 0.0
	prod := 1.0
	for i, v := range x {
		sum += v * v / 4000.0
		prod *= math.Cos(v / math.Sqrt(float64(i+1)))
	}
	return 1 + sum - prod
}

// Schwefel is a deceptive multimodal benchmark
// f(x) = 418.9829 n - sum x_i sin(sqrt(|x_i|)). Its global minimum is
// approximately 0 at x_i = 420.9687; it is evaluated on [-500, 500]^n.
func Schwefel(x []float64) float64 {
	s := 418.9828872724338 * float64(len(x))
	for _, v := range x {
		s -= v * math.Sin(math.Sqrt(math.Abs(v)))
	}
	return s
}

// Zakharov is a unimodal benchmark
// f(x) = sum x_i^2 + (sum 0.5 i x_i)^2 + (sum 0.5 i x_i)^4. Its global minimum
// is 0 at the origin.
func Zakharov(x []float64) float64 {
	var s1, s2 float64
	for i, v := range x {
		s1 += v * v
		s2 += 0.5 * float64(i+1) * v
	}
	return s1 + s2*s2 + s2*s2*s2*s2
}

// StyblinskiTang is a multimodal benchmark
// f(x) = 0.5 sum(x_i^4 - 16 x_i^2 + 5 x_i). Its global minimum is
// -39.16599 n at x_i = -2.903534.
func StyblinskiTang(x []float64) float64 {
	s := 0.0
	for _, v := range x {
		s += v*v*v*v - 16*v*v + 5*v
	}
	return 0.5 * s
}

// DixonPrice is a unimodal benchmark
// f(x) = (x_1-1)^2 + sum_{i>=2} i (2 x_i^2 - x_{i-1})^2. Its global minimum is
// 0 at x_i = 2^(-(2^i-2)/2^i).
func DixonPrice(x []float64) float64 {
	if len(x) == 0 {
		return 0
	}
	s := (x[0] - 1) * (x[0] - 1)
	for i := 1; i < len(x); i++ {
		t := 2*x[i]*x[i] - x[i-1]
		s += float64(i+1) * t * t
	}
	return s
}

// SumSquares is the weighted sum-of-squares benchmark
// f(x) = sum (i+1) x_i^2. Its global minimum is 0 at the origin.
func SumSquares(x []float64) float64 {
	s := 0.0
	for i, v := range x {
		s += float64(i+1) * v * v
	}
	return s
}

// SumOfDifferentPowers is the benchmark f(x) = sum |x_i|^(i+2). Its global
// minimum is 0 at the origin.
func SumOfDifferentPowers(x []float64) float64 {
	s := 0.0
	for i, v := range x {
		s += math.Pow(math.Abs(v), float64(i+2))
	}
	return s
}

// Trid is the Trid benchmark
// f(x) = sum (x_i-1)^2 - sum_{i>=2} x_i x_{i-1}. Its global minimum is
// -n(n+4)(n-1)/6 at x_i = (i+1)(n-i).
func Trid(x []float64) float64 {
	var s1, s2 float64
	for i, v := range x {
		d := v - 1
		s1 += d * d
		if i > 0 {
			s2 += v * x[i-1]
		}
	}
	return s1 - s2
}

// Sphere shifted variants and 2-D functions follow.

// Booth is the 2-D Booth benchmark
// f(x,y) = (x+2y-7)^2 + (2x+y-5)^2 with global minimum 0 at (1, 3).
func Booth(x []float64) float64 {
	a := x[0] + 2*x[1] - 7
	b := 2*x[0] + x[1] - 5
	return a*a + b*b
}

// Matyas is the 2-D Matyas benchmark
// f(x,y) = 0.26(x^2+y^2) - 0.48 x y with global minimum 0 at the origin.
func Matyas(x []float64) float64 {
	return 0.26*(x[0]*x[0]+x[1]*x[1]) - 0.48*x[0]*x[1]
}

// Beale is the 2-D Beale benchmark with global minimum 0 at (3, 0.5).
func Beale(x []float64) float64 {
	a := 1.5 - x[0] + x[0]*x[1]
	b := 2.25 - x[0] + x[0]*x[1]*x[1]
	c := 2.625 - x[0] + x[0]*x[1]*x[1]*x[1]
	return a*a + b*b + c*c
}

// Himmelblau is the 2-D Himmelblau benchmark
// f(x,y) = (x^2+y-11)^2 + (x+y^2-7)^2. It has four identical global minima of
// value 0, one of which is (3, 2).
func Himmelblau(x []float64) float64 {
	a := x[0]*x[0] + x[1] - 11
	b := x[0] + x[1]*x[1] - 7
	return a*a + b*b
}

// ThreeHumpCamel is the 2-D three-hump camel benchmark with global minimum 0
// at the origin.
func ThreeHumpCamel(x []float64) float64 {
	xx := x[0] * x[0]
	return 2*xx - 1.05*xx*xx + xx*xx*xx/6 + x[0]*x[1] + x[1]*x[1]
}

// SixHumpCamel is the 2-D six-hump camel benchmark with two global minima of
// value approximately -1.0316.
func SixHumpCamel(x []float64) float64 {
	xx := x[0] * x[0]
	yy := x[1] * x[1]
	return (4-2.1*xx+xx*xx/3)*xx + x[0]*x[1] + (-4+4*yy)*yy
}

// Easom is the 2-D Easom benchmark with a single narrow global minimum of -1
// at (pi, pi).
func Easom(x []float64) float64 {
	return -math.Cos(x[0]) * math.Cos(x[1]) *
		math.Exp(-((x[0]-math.Pi)*(x[0]-math.Pi) + (x[1]-math.Pi)*(x[1]-math.Pi)))
}

// LeviN13 is the 2-D Levi function N.13 with global minimum 0 at (1, 1).
func LeviN13(x []float64) float64 {
	p := math.Pi
	a := math.Sin(3 * p * x[0])
	b := (x[0] - 1) * (x[0] - 1) * (1 + math.Sin(3*p*x[1])*math.Sin(3*p*x[1]))
	c := (x[1] - 1) * (x[1] - 1) * (1 + math.Sin(2*p*x[1])*math.Sin(2*p*x[1]))
	return a*a + b + c
}

// GoldsteinPrice is the 2-D Goldstein-Price benchmark with global minimum 3 at
// (0, -1).
func GoldsteinPrice(x []float64) float64 {
	a := x[0]
	b := x[1]
	t1 := 1 + (a+b+1)*(a+b+1)*(19-14*a+3*a*a-14*b+6*a*b+3*b*b)
	t2 := 30 + (2*a-3*b)*(2*a-3*b)*(18-32*a+12*a*a+48*b-36*a*b+27*b*b)
	return t1 * t2
}

// McCormick is the 2-D McCormick benchmark with global minimum approximately
// -1.9133 at (-0.54719, -1.54719).
func McCormick(x []float64) float64 {
	return math.Sin(x[0]+x[1]) + (x[0]-x[1])*(x[0]-x[1]) - 1.5*x[0] + 2.5*x[1] + 1
}

// SchafferN2 is the 2-D Schaffer function N.2 with global minimum 0 at the
// origin.
func SchafferN2(x []float64) float64 {
	num := math.Sin(x[0]*x[0]-x[1]*x[1])*math.Sin(x[0]*x[0]-x[1]*x[1]) - 0.5
	den := 1 + 0.001*(x[0]*x[0]+x[1]*x[1])
	return 0.5 + num/(den*den)
}

// SchafferN4 is the 2-D Schaffer function N.4 with global minimum approximately
// 0.292579 at (0, 1.25313).
func SchafferN4(x []float64) float64 {
	num := math.Cos(math.Sin(math.Abs(x[0]*x[0]-x[1]*x[1])))*math.Cos(math.Sin(math.Abs(x[0]*x[0]-x[1]*x[1]))) - 0.5
	den := 1 + 0.001*(x[0]*x[0]+x[1]*x[1])
	return 0.5 + num/(den*den)
}

// DropWave is the 2-D drop-wave benchmark with global minimum -1 at the origin.
func DropWave(x []float64) float64 {
	r := x[0]*x[0] + x[1]*x[1]
	return -(1 + math.Cos(12*math.Sqrt(r))) / (0.5*r + 2)
}

// Bukin6 is the 2-D Bukin function N.6 with global minimum 0 at (-10, 1).
func Bukin6(x []float64) float64 {
	return 100*math.Sqrt(math.Abs(x[1]-0.01*x[0]*x[0])) + 0.01*math.Abs(x[0]+10)
}

// Eggholder is the 2-D Eggholder benchmark with global minimum approximately
// -959.6407 at (512, 404.2319).
func Eggholder(x []float64) float64 {
	return -(x[1]+47)*math.Sin(math.Sqrt(math.Abs(x[0]/2+x[1]+47))) -
		x[0]*math.Sin(math.Sqrt(math.Abs(x[0]-x[1]-47)))
}

// HolderTable is the 2-D Holder table benchmark with four global minima of
// value approximately -19.2085.
func HolderTable(x []float64) float64 {
	return -math.Abs(math.Sin(x[0]) * math.Cos(x[1]) *
		math.Exp(math.Abs(1-math.Sqrt(x[0]*x[0]+x[1]*x[1])/math.Pi)))
}

// CrossInTray is the 2-D cross-in-tray benchmark with four global minima of
// value approximately -2.06261.
func CrossInTray(x []float64) float64 {
	e := math.Abs(100 - math.Sqrt(x[0]*x[0]+x[1]*x[1])/math.Pi)
	t := math.Abs(math.Sin(x[0])*math.Sin(x[1])*math.Exp(e)) + 1
	return -0.0001 * math.Pow(t, 0.1)
}

// Michalewicz is the multimodal Michalewicz benchmark with steepness parameter
// m = 10: f(x) = -sum sin(x_i) sin^{2m}((i+1) x_i^2 / pi). It is evaluated on
// [0, pi]^n and has many local minima.
func Michalewicz(x []float64) float64 {
	const m = 10.0
	s := 0.0
	for i, v := range x {
		s += math.Sin(v) * math.Pow(math.Sin(float64(i+1)*v*v/math.Pi), 2*m)
	}
	return -s
}

// Powell is the 4k-dimensional Powell benchmark with global minimum 0 at the
// origin. It requires the dimension to be a multiple of 4; extra trailing
// coordinates are ignored.
func Powell(x []float64) float64 {
	s := 0.0
	for i := 0; i+3 < len(x); i += 4 {
		a := x[i] + 10*x[i+1]
		b := x[i+2] - x[i+3]
		c := x[i+1] - 2*x[i+2]
		d := x[i] - x[i+3]
		s += a*a + 5*b*b + c*c*c*c + 10*d*d*d*d
	}
	return s
}

// Benchmark bundles a benchmark objective with metadata: a human-readable
// name, a recommended search box for a given dimension, and its known global
// minimum value.
type Benchmark struct {
	Name string
	Func ObjectiveFunc
	// BoundsFor returns a recommended search box of the given dimension.
	BoundsFor func(dim int) Bounds
	// GlobalMin returns the global minimum value for the given dimension.
	GlobalMin func(dim int) float64
	// FixedDim, when positive, is the only dimension the function accepts.
	FixedDim int
}

// StandardBenchmarks returns a slice of the scalable n-dimensional benchmark
// functions (Sphere, Rastrigin, Ackley, Rosenbrock, Griewank, Schwefel,
// Zakharov, Styblinski-Tang and Sum-of-Squares) with their metadata.
func StandardBenchmarks() []Benchmark {
	return []Benchmark{
		{
			Name:      "Sphere",
			Func:      Sphere,
			BoundsFor: func(d int) Bounds { return UniformBounds(d, -5.12, 5.12) },
			GlobalMin: func(int) float64 { return 0 },
		},
		{
			Name:      "Rastrigin",
			Func:      Rastrigin,
			BoundsFor: func(d int) Bounds { return UniformBounds(d, -5.12, 5.12) },
			GlobalMin: func(int) float64 { return 0 },
		},
		{
			Name:      "Ackley",
			Func:      Ackley,
			BoundsFor: func(d int) Bounds { return UniformBounds(d, -32.768, 32.768) },
			GlobalMin: func(int) float64 { return 0 },
		},
		{
			Name:      "Rosenbrock",
			Func:      Rosenbrock,
			BoundsFor: func(d int) Bounds { return UniformBounds(d, -5, 10) },
			GlobalMin: func(int) float64 { return 0 },
		},
		{
			Name:      "Griewank",
			Func:      Griewank,
			BoundsFor: func(d int) Bounds { return UniformBounds(d, -600, 600) },
			GlobalMin: func(int) float64 { return 0 },
		},
		{
			Name:      "Schwefel",
			Func:      Schwefel,
			BoundsFor: func(d int) Bounds { return UniformBounds(d, -500, 500) },
			GlobalMin: func(int) float64 { return 0 },
		},
		{
			Name:      "Zakharov",
			Func:      Zakharov,
			BoundsFor: func(d int) Bounds { return UniformBounds(d, -5, 10) },
			GlobalMin: func(int) float64 { return 0 },
		},
		{
			Name:      "StyblinskiTang",
			Func:      StyblinskiTang,
			BoundsFor: func(d int) Bounds { return UniformBounds(d, -5, 5) },
			GlobalMin: func(d int) float64 { return -39.16599 * float64(d) },
		},
		{
			Name:      "SumSquares",
			Func:      SumSquares,
			BoundsFor: func(d int) Bounds { return UniformBounds(d, -10, 10) },
			GlobalMin: func(int) float64 { return 0 },
		},
	}
}
