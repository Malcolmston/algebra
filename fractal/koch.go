package fractal

import "math"

// fractalKochPeak returns the apex of the equilateral bump raised on the middle
// third of the segment from a to b, rotated 60 degrees counterclockwise so the
// bump lies to the left of the a->b direction.
func fractalKochPeak(c, d Point2D) Point2D {
	v := d.Sub(c)
	// Rotate v by +60 degrees.
	const cos60, sin60 = 0.5, 0.8660254037844386
	rot := Point2D{v.X*cos60 - v.Y*sin60, v.X*sin60 + v.Y*cos60}
	return c.Add(rot)
}

// KochCurve returns the polyline vertices of the Koch curve constructed between
// the endpoints p0 and p1 after the given number of iterations. Iteration 0
// returns the two endpoints [p0, p1]; each further iteration replaces every
// segment by four segments (a bump raised on its middle third), so the result
// has 4^iterations + 1 points. Bumps are raised to the left of the p0->p1
// direction. It panics if iterations is negative.
func KochCurve(p0, p1 Point2D, iterations int) []Point2D {
	if iterations < 0 {
		panic("fractal: KochCurve needs non-negative iterations")
	}
	pts := []Point2D{p0, p1}
	for it := 0; it < iterations; it++ {
		next := make([]Point2D, 0, (len(pts)-1)*4+1)
		for i := 0; i < len(pts)-1; i++ {
			a, b := pts[i], pts[i+1]
			c := a.Lerp(b, 1.0/3.0)
			d := a.Lerp(b, 2.0/3.0)
			peak := fractalKochPeak(c, d)
			next = append(next, a, c, peak, d)
		}
		next = append(next, pts[len(pts)-1])
		pts = next
	}
	return pts
}

// KochSnowflake returns the closed polyline of the Koch snowflake inscribed in
// the circle of the given radius centered at center, after the given number of
// iterations. The base is an equilateral triangle whose corners lie on the
// circle; each edge is replaced by a Koch curve with bumps pointing outward.
// The returned slice is closed (its last point equals its first). It panics if
// iterations is negative.
func KochSnowflake(center Point2D, radius float64, iterations int) []Point2D {
	if iterations < 0 {
		panic("fractal: KochSnowflake needs non-negative iterations")
	}
	// Three vertices at 90, 210, 330 degrees. Traverse clockwise so the
	// left-hand Koch bumps point outward.
	angles := []float64{90, 330, 210}
	verts := make([]Point2D, 3)
	for i, deg := range angles {
		rad := deg * math.Pi / 180
		verts[i] = Point2D{center.X + radius*math.Cos(rad), center.Y + radius*math.Sin(rad)}
	}
	var out []Point2D
	for i := 0; i < 3; i++ {
		a := verts[i]
		b := verts[(i+1)%3]
		edge := KochCurve(a, b, iterations)
		// Drop the last point of each edge to avoid duplicating the shared
		// vertex; the closing point is appended at the end.
		out = append(out, edge[:len(edge)-1]...)
	}
	out = append(out, out[0])
	return out
}

// Triangle is a triangle with vertices A, B and C.
type Triangle struct {
	A, B, C Point2D
}

// SierpinskiTriangle returns the list of filled sub-triangles produced by
// recursively removing the central triangle from the triangle (a, b, c) to the
// given depth. Depth 0 returns the single input triangle; each level replaces
// every triangle by its three corner sub-triangles, so the result contains
// 3^depth triangles. It panics if depth is negative.
func SierpinskiTriangle(a, b, c Point2D, depth int) []Triangle {
	if depth < 0 {
		panic("fractal: SierpinskiTriangle needs non-negative depth")
	}
	tris := []Triangle{{a, b, c}}
	for d := 0; d < depth; d++ {
		next := make([]Triangle, 0, len(tris)*3)
		for _, t := range tris {
			ab := t.A.Midpoint(t.B)
			bc := t.B.Midpoint(t.C)
			ca := t.C.Midpoint(t.A)
			next = append(next,
				Triangle{t.A, ab, ca},
				Triangle{ab, t.B, bc},
				Triangle{ca, bc, t.C},
			)
		}
		tris = next
	}
	return tris
}

// Interval is a closed real interval [Start, End].
type Interval struct {
	Start, End float64
}

// Length returns End - Start.
func (iv Interval) Length() float64 { return iv.End - iv.Start }

// CantorSet returns the intervals remaining after iterations steps of the
// middle-thirds Cantor construction applied to [start, end]. Iteration 0
// returns the single interval [start, end]; each step replaces every interval
// by its first and last thirds, so the result contains 2^iterations intervals
// with total length (end-start)*(2/3)^iterations. Intervals are returned in
// increasing order. It panics if iterations is negative.
func CantorSet(start, end float64, iterations int) []Interval {
	if iterations < 0 {
		panic("fractal: CantorSet needs non-negative iterations")
	}
	ivs := []Interval{{start, end}}
	for it := 0; it < iterations; it++ {
		next := make([]Interval, 0, len(ivs)*2)
		for _, iv := range ivs {
			third := iv.Length() / 3
			next = append(next,
				Interval{iv.Start, iv.Start + third},
				Interval{iv.End - third, iv.End},
			)
		}
		ivs = next
	}
	return ivs
}
