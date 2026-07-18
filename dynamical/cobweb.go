package dynamical

// CobwebSegment is one line segment of a cobweb (staircase) diagram, joining
// the point (X0, Y0) to (X1, Y1). Cobweb diagrams alternate vertical segments
// (from the diagonal up or down to the graph of f) with horizontal segments
// (from the graph across to the diagonal), visualizing the iteration of a
// one-dimensional map.
type CobwebSegment struct {
	X0, Y0, X1, Y1 float64
}

// Cobweb returns the segments of the cobweb diagram for n iterations of the
// map f starting at x0. The construction begins on the diagonal at (x0, x0)
// and, for each step, adds a vertical segment up to the graph (x, f(x))
// followed by a horizontal segment across to the diagonal (f(x), f(x)). The
// result therefore contains 2*n segments (for n >= 0).
func Cobweb(f Map1D, x0 float64, n int) []CobwebSegment {
	if n < 0 {
		n = 0
	}
	segs := make([]CobwebSegment, 0, 2*n)
	x := x0
	for i := 0; i < n; i++ {
		y := f(x)
		// Vertical: from (x, x) on the diagonal up/down to (x, y) on the graph.
		segs = append(segs, CobwebSegment{X0: x, Y0: x, X1: x, Y1: y})
		// Horizontal: from (x, y) across to (y, y) on the diagonal.
		segs = append(segs, CobwebSegment{X0: x, Y0: y, X1: y, Y1: y})
		x = y
	}
	return segs
}
