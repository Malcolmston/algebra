package fractal

import "math/rand"

// AffineMap is a 2-D affine transformation of the form
//
//	x' = A*x + B*y + E
//	y' = C*x + D*y + F
//
// Such contraction maps are the building blocks of iterated function systems.
type AffineMap struct {
	A, B, C, D, E, F float64
}

// NewAffineMap constructs an AffineMap from its six coefficients in the order
// a, b, c, d, e, f matching the field layout x' = a*x+b*y+e, y' = c*x+d*y+f.
func NewAffineMap(a, b, c, d, e, f float64) AffineMap {
	return AffineMap{a, b, c, d, e, f}
}

// Apply returns the image of p under the affine map.
func (m AffineMap) Apply(p Point2D) Point2D {
	return Point2D{
		X: m.A*p.X + m.B*p.Y + m.E,
		Y: m.C*p.X + m.D*p.Y + m.F,
	}
}

// IFS is an iterated function system: a list of affine contraction Maps with an
// equal-length list of selection Weights. When run as a chaos game the map at
// index i is chosen with probability Weights[i] / sum(Weights). If Weights is
// nil or empty the maps are chosen uniformly.
type IFS struct {
	Maps    []AffineMap
	Weights []float64
}

// fractalCumulative returns the normalized cumulative distribution of the
// system's weights, or nil when the maps should be selected uniformly.
func (s IFS) fractalCumulative() []float64 {
	if len(s.Weights) == 0 {
		return nil
	}
	cum := make([]float64, len(s.Maps))
	var total float64
	for i := range s.Maps {
		w := 0.0
		if i < len(s.Weights) {
			w = s.Weights[i]
		}
		total += w
		cum[i] = total
	}
	if total <= 0 {
		return nil
	}
	for i := range cum {
		cum[i] /= total
	}
	return cum
}

// ChaosGameFrom runs the random-iteration ("chaos game") algorithm starting
// from the point start: it repeatedly picks a map according to the system's
// weights and applies it, discarding the first transient points before
// recording n points. Randomness is driven solely by seed, so the output is
// fully deterministic for a given seed. It returns n points and panics if the
// system has no maps or n is negative.
func (s IFS) ChaosGameFrom(start Point2D, n, transient int, seed int64) []Point2D {
	if len(s.Maps) == 0 {
		panic("fractal: ChaosGame needs at least one map")
	}
	if n < 0 {
		panic("fractal: ChaosGame needs non-negative n")
	}
	if transient < 0 {
		transient = 0
	}
	rng := rand.New(rand.NewSource(seed))
	cum := s.fractalCumulative()
	p := start
	pick := func() AffineMap {
		if cum == nil {
			return s.Maps[rng.Intn(len(s.Maps))]
		}
		u := rng.Float64()
		for i, c := range cum {
			if u <= c {
				return s.Maps[i]
			}
		}
		return s.Maps[len(s.Maps)-1]
	}
	for i := 0; i < transient; i++ {
		p = pick().Apply(p)
	}
	out := make([]Point2D, n)
	for i := 0; i < n; i++ {
		p = pick().Apply(p)
		out[i] = p
	}
	return out
}

// ChaosGame runs the chaos game starting from the origin, discarding 20
// transient points, and returns n attractor points. It is a convenience wrapper
// around [IFS.ChaosGameFrom]. The output is deterministic for a given seed.
func (s IFS) ChaosGame(n int, seed int64) []Point2D {
	return s.ChaosGameFrom(Point2D{0, 0}, n, 20, seed)
}

// SierpinskiTriangleIFS returns the three-map iterated function system whose
// attractor is the Sierpinski triangle with vertices (0,0), (1,0) and
// (0.5, 0.5). Each map scales the plane by 1/2 toward one vertex; the maps are
// chosen with equal probability.
func SierpinskiTriangleIFS() IFS {
	return IFS{
		Maps: []AffineMap{
			{0.5, 0, 0, 0.5, 0, 0},
			{0.5, 0, 0, 0.5, 0.5, 0},
			{0.5, 0, 0, 0.5, 0.25, 0.25},
		},
		Weights: []float64{1, 1, 1},
	}
}

// SierpinskiCarpetIFS returns the eight-map iterated function system whose
// attractor is the Sierpinski carpet: the unit square is divided into a 3x3
// grid and every cell except the center maps back to the whole via a 1/3-scale
// contraction. The maps are chosen with equal probability.
func SierpinskiCarpetIFS() IFS {
	var maps []AffineMap
	for gy := 0; gy < 3; gy++ {
		for gx := 0; gx < 3; gx++ {
			if gx == 1 && gy == 1 {
				continue
			}
			maps = append(maps, AffineMap{
				1.0 / 3, 0, 0, 1.0 / 3,
				float64(gx) / 3, float64(gy) / 3,
			})
		}
	}
	w := make([]float64, len(maps))
	for i := range w {
		w[i] = 1
	}
	return IFS{Maps: maps, Weights: w}
}

// BarnsleyFern returns the four-map iterated function system that generates
// Michael Barnsley's fern. The maps and their selection probabilities are the
// classical values; the attractor fits within roughly x in [-2.2, 2.7] and
// y in [0, 10].
func BarnsleyFern() IFS {
	return IFS{
		Maps: []AffineMap{
			{0, 0, 0, 0.16, 0, 0},
			{0.85, 0.04, -0.04, 0.85, 0, 1.6},
			{0.2, -0.26, 0.23, 0.22, 0, 1.6},
			{-0.15, 0.28, 0.26, 0.24, 0, 0.44},
		},
		Weights: []float64{0.01, 0.85, 0.07, 0.07},
	}
}
