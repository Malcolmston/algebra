package fem

import "math"

// QuadRule1D is a quadrature rule on the reference interval [-1, 1]. The nodes
// and weights integrate polynomials exactly up to a rule-dependent degree.
type QuadRule1D struct {
	Nodes   []float64
	Weights []float64
}

// GaussLegendre1D returns the n-point Gauss–Legendre rule on [-1, 1], exact for
// polynomials of degree up to 2n-1. It panics if n < 1.
func GaussLegendre1D(n int) QuadRule1D {
	if n < 1 {
		panic("fem: GaussLegendre1D requires n >= 1")
	}
	nodes := make([]float64, n)
	weights := make([]float64, n)
	for i := 0; i < n; i++ {
		// Initial guess for the i-th root of the Legendre polynomial.
		x := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		var pp float64
		for iter := 0; iter < 100; iter++ {
			p0, p1 := 1.0, 0.0
			for k := 0; k < n; k++ {
				p2 := p1
				p1 = p0
				p0 = ((2*float64(k)+1)*x*p1 - float64(k)*p2) / (float64(k) + 1)
			}
			// p0 = P_n(x); derivative pp.
			pp = float64(n) * (x*p0 - p1) / (x*x - 1)
			dx := p0 / pp
			x -= dx
			if math.Abs(dx) < 1e-15 {
				break
			}
		}
		nodes[i] = x
		weights[i] = 2 / ((1 - x*x) * pp * pp)
	}
	return QuadRule1D{Nodes: nodes, Weights: weights}
}

// MapToInterval returns nodes and weights of the rule mapped from [-1,1] to the
// interval [a, b].
func (q QuadRule1D) MapToInterval(a, b float64) (nodes, weights []float64) {
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	nodes = make([]float64, len(q.Nodes))
	weights = make([]float64, len(q.Weights))
	for i := range q.Nodes {
		nodes[i] = mid + half*q.Nodes[i]
		weights[i] = half * q.Weights[i]
	}
	return nodes, weights
}

// IntegrateInterval approximates the integral of f over [a, b] using the
// n-point Gauss–Legendre rule.
func IntegrateInterval(f func(float64) float64, a, b float64, n int) float64 {
	q := GaussLegendre1D(n)
	nodes, weights := q.MapToInterval(a, b)
	var s float64
	for i := range nodes {
		s += weights[i] * f(nodes[i])
	}
	return s
}

// TriQuadRule is a symmetric quadrature rule on a triangle. Bary holds the
// barycentric coordinates of each point and Weights holds the corresponding
// weights, normalised to sum to one. The physical integral over a triangle of
// area A is A times the weighted sum of the integrand at the mapped points.
type TriQuadRule struct {
	Bary    [][3]float64
	Weights []float64
}

// TriangleQuadrature returns a symmetric quadrature rule on the reference
// triangle exact for polynomials up to the requested degree (1, 2, 3 or 5; a
// request for degree 4 is served by the degree-5 rule). It panics on an
// unsupported degree.
func TriangleQuadrature(degree int) TriQuadRule {
	switch degree {
	case 1:
		return TriQuadRule{
			Bary:    [][3]float64{{1.0 / 3, 1.0 / 3, 1.0 / 3}},
			Weights: []float64{1},
		}
	case 2:
		return TriQuadRule{
			Bary: [][3]float64{
				{2.0 / 3, 1.0 / 6, 1.0 / 6},
				{1.0 / 6, 2.0 / 3, 1.0 / 6},
				{1.0 / 6, 1.0 / 6, 2.0 / 3},
			},
			Weights: []float64{1.0 / 3, 1.0 / 3, 1.0 / 3},
		}
	case 3:
		w0 := -27.0 / 48.0
		w1 := 25.0 / 48.0
		return TriQuadRule{
			Bary: [][3]float64{
				{1.0 / 3, 1.0 / 3, 1.0 / 3},
				{0.6, 0.2, 0.2},
				{0.2, 0.6, 0.2},
				{0.2, 0.2, 0.6},
			},
			Weights: []float64{w0, w1, w1, w1},
		}
	case 4, 5:
		a1 := 0.059715871789770
		b1 := 0.470142064105115
		a2 := 0.797426985353087
		b2 := 0.101286507323456
		w0 := 0.225
		w1 := 0.132394152788506
		w2 := 0.125939180544827
		return TriQuadRule{
			Bary: [][3]float64{
				{1.0 / 3, 1.0 / 3, 1.0 / 3},
				{a1, b1, b1},
				{b1, a1, b1},
				{b1, b1, a1},
				{a2, b2, b2},
				{b2, a2, b2},
				{b2, b2, a2},
			},
			Weights: []float64{w0, w1, w1, w1, w2, w2, w2},
		}
	default:
		panic("fem: unsupported triangle quadrature degree")
	}
}

// IntegrateTriangle approximates the integral of f over the triangle with the
// given vertices using a rule of the requested degree. f is evaluated at
// physical coordinates.
func IntegrateTriangle(f func(x, y float64) float64, v1, v2, v3 [2]float64, degree int) float64 {
	rule := TriangleQuadrature(degree)
	area := TriangleArea(v1, v2, v3)
	var s float64
	for i, bc := range rule.Bary {
		x := bc[0]*v1[0] + bc[1]*v2[0] + bc[2]*v3[0]
		y := bc[0]*v1[1] + bc[1]*v2[1] + bc[2]*v3[1]
		s += rule.Weights[i] * f(x, y)
	}
	return area * s
}

// TriangleArea returns the (unsigned) area of the triangle with the given
// vertices.
func TriangleArea(v1, v2, v3 [2]float64) float64 {
	return 0.5 * math.Abs((v2[0]-v1[0])*(v3[1]-v1[1])-(v3[0]-v1[0])*(v2[1]-v1[1]))
}

// TriangleSignedArea returns the signed area of the triangle, positive when the
// vertices are counter-clockwise.
func TriangleSignedArea(v1, v2, v3 [2]float64) float64 {
	return 0.5 * ((v2[0]-v1[0])*(v3[1]-v1[1]) - (v3[0]-v1[0])*(v2[1]-v1[1]))
}
