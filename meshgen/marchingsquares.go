package meshgen

import "errors"

// ScalarGrid2 is a scalar field sampled on a regular rectangular lattice of
// Nx by Ny nodes. Node (i, j) sits at (X0 + i*Dx, Y0 + j*Dy) and its value is
// stored row-major at Values[j*Nx+i].
type ScalarGrid2 struct {
	Nx, Ny         int
	X0, Y0, Dx, Dy float64
	Values         []float64
}

// NewScalarGrid2 returns a zero-initialised grid of nx by ny nodes with the
// given origin and spacing. It panics if nx or ny is less than two or the
// spacing is not positive.
func NewScalarGrid2(nx, ny int, x0, y0, dx, dy float64) *ScalarGrid2 {
	if nx < 2 || ny < 2 {
		panic("meshgen: grid needs at least 2x2 nodes")
	}
	if dx <= 0 || dy <= 0 {
		panic("meshgen: grid spacing must be positive")
	}
	return &ScalarGrid2{
		Nx: nx, Ny: ny, X0: x0, Y0: y0, Dx: dx, Dy: dy,
		Values: make([]float64, nx*ny),
	}
}

// SampleGrid2 builds a scalar grid by evaluating f at every node of the
// nx by ny lattice.
func SampleGrid2(nx, ny int, x0, y0, dx, dy float64, f func(x, y float64) float64) *ScalarGrid2 {
	g := NewScalarGrid2(nx, ny, x0, y0, dx, dy)
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			g.Set(i, j, f(g.X0+float64(i)*g.Dx, g.Y0+float64(j)*g.Dy))
		}
	}
	return g
}

// At returns the value at node (i, j).
func (g *ScalarGrid2) At(i, j int) float64 { return g.Values[j*g.Nx+i] }

// Set assigns v to node (i, j).
func (g *ScalarGrid2) Set(i, j int, v float64) { g.Values[j*g.Nx+i] = v }

// Node returns the position of node (i, j).
func (g *ScalarGrid2) Node(i, j int) Vec2 {
	return Vec2{g.X0 + float64(i)*g.Dx, g.Y0 + float64(j)*g.Dy}
}

// MarchingSquares extracts the iso-contour of the grid at the given isovalue as
// a set of line segments. Saddle cells are resolved with the asymptotic
// (mean-value) decider so contours never cross inside a cell.
func MarchingSquares(g *ScalarGrid2, iso float64) []Segment {
	var out []Segment
	for j := 0; j < g.Ny-1; j++ {
		for i := 0; i < g.Nx-1; i++ {
			out = append(out, marchingSquaresCell(g, i, j, iso)...)
		}
	}
	return out
}

// marchingSquaresCell returns the contour segments inside a single cell whose
// lower-left corner is node (i, j).
func marchingSquaresCell(g *ScalarGrid2, i, j int, iso float64) []Segment {
	// Corners: 0=bl, 1=br, 2=tr, 3=tl.
	pos := [4]Vec2{
		g.Node(i, j),
		g.Node(i+1, j),
		g.Node(i+1, j+1),
		g.Node(i, j+1),
	}
	val := [4]float64{
		g.At(i, j),
		g.At(i+1, j),
		g.At(i+1, j+1),
		g.At(i, j+1),
	}
	inside := func(k int) bool { return val[k] >= iso }
	// Edge k joins corner k and corner (k+1)%4.
	crossPoint := func(k int) Vec2 {
		a := k
		b := (k + 1) % 4
		return isoInterp(pos[a], pos[b], val[a], val[b], iso)
	}
	var crossing []int
	for k := 0; k < 4; k++ {
		if inside(k) != inside((k+1)%4) {
			crossing = append(crossing, k)
		}
	}
	switch len(crossing) {
	case 2:
		return []Segment{{crossPoint(crossing[0]), crossPoint(crossing[1])}}
	case 4:
		center := (val[0] + val[1] + val[2] + val[3]) / 4
		p0 := crossPoint(0)
		p1 := crossPoint(1)
		p2 := crossPoint(2)
		p3 := crossPoint(3)
		if (center >= iso) == inside(0) {
			return []Segment{{p3, p0}, {p1, p2}}
		}
		return []Segment{{p0, p1}, {p2, p3}}
	default:
		return nil
	}
}

// isoInterp returns the point on the segment a-b where the linearly
// interpolated field equals iso. It falls back to the midpoint when the two
// endpoint values are equal.
func isoInterp(a, b Vec2, va, vb, iso float64) Vec2 {
	den := vb - va
	if den == 0 {
		return a.Midpoint(b)
	}
	t := (iso - va) / den
	return a.Lerp(b, t)
}

// IsoContour is a convenience wrapper that samples f on an nx by ny lattice and
// returns the marching-squares contour at the given isovalue.
func IsoContour(nx, ny int, x0, y0, dx, dy float64, f func(x, y float64) float64, iso float64) []Segment {
	return MarchingSquares(SampleGrid2(nx, ny, x0, y0, dx, dy, f), iso)
}

// TotalSegmentLength returns the sum of the lengths of the given segments.
func TotalSegmentLength(segs []Segment) float64 {
	var s float64
	for _, seg := range segs {
		s += seg.Length()
	}
	return s
}

// ContourError is returned by contour helpers on invalid input.
var ContourError = errors.New("meshgen: invalid contour input")
