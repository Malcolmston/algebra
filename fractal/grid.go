package fractal

import "math"

// Grid is a dense Width×Height array of scalar values stored in row-major order
// (index = y*Width + x). It is the natural output of rasterizing an
// escape-time fractal over a rectangular region.
type Grid struct {
	Width, Height int
	Data          []float64
}

// NewGrid allocates a zero-filled Grid with the given dimensions. It panics if
// width or height is negative.
func NewGrid(width, height int) *Grid {
	if width < 0 || height < 0 {
		panic("fractal: negative grid dimension")
	}
	return &Grid{Width: width, Height: height, Data: make([]float64, width*height)}
}

// At returns the value at column x, row y. It panics if the coordinates are out
// of range.
func (g *Grid) At(x, y int) float64 {
	if x < 0 || x >= g.Width || y < 0 || y >= g.Height {
		panic("fractal: grid index out of range")
	}
	return g.Data[y*g.Width+x]
}

// Set stores v at column x, row y. It panics if the coordinates are out of
// range.
func (g *Grid) Set(x, y int, v float64) {
	if x < 0 || x >= g.Width || y < 0 || y >= g.Height {
		panic("fractal: grid index out of range")
	}
	g.Data[y*g.Width+x] = v
}

// Row returns the y-th row of the grid as a slice that aliases the underlying
// storage. Mutating it mutates the grid.
func (g *Grid) Row(y int) []float64 {
	if y < 0 || y >= g.Height {
		panic("fractal: grid row out of range")
	}
	return g.Data[y*g.Width : y*g.Width+g.Width]
}

// MinMax returns the minimum and maximum values in the grid. For an empty grid
// it returns (0, 0).
func (g *Grid) MinMax() (min, max float64) {
	if len(g.Data) == 0 {
		return 0, 0
	}
	min, max = g.Data[0], g.Data[0]
	for _, v := range g.Data[1:] {
		if v < min {
			min = v
		}
		if v > max {
			max = v
		}
	}
	return min, max
}

// Min returns the smallest value in the grid, or 0 for an empty grid.
func (g *Grid) Min() float64 {
	min, _ := g.MinMax()
	return min
}

// Max returns the largest value in the grid, or 0 for an empty grid.
func (g *Grid) Max() float64 {
	_, max := g.MinMax()
	return max
}

// Normalize returns a new Grid whose values are linearly rescaled so that the
// original minimum maps to 0 and the original maximum maps to 1. If all values
// are equal the result is all zeros. The receiver is not modified.
func (g *Grid) Normalize() *Grid {
	out := NewGrid(g.Width, g.Height)
	min, max := g.MinMax()
	span := max - min
	if span == 0 {
		return out
	}
	inv := 1.0 / span
	for i, v := range g.Data {
		out.Data[i] = (v - min) * inv
	}
	return out
}

// Clone returns a deep copy of the grid.
func (g *Grid) Clone() *Grid {
	out := NewGrid(g.Width, g.Height)
	copy(out.Data, g.Data)
	return out
}

// Mean returns the arithmetic mean of all grid values, or 0 for an empty grid.
func (g *Grid) Mean() float64 {
	if len(g.Data) == 0 {
		return 0
	}
	var s float64
	for _, v := range g.Data {
		s += v
	}
	return s / float64(len(g.Data))
}

// CountFinite returns the number of grid cells whose value is finite (neither
// NaN nor infinite).
func (g *Grid) CountFinite() int {
	n := 0
	for _, v := range g.Data {
		if !math.IsNaN(v) && !math.IsInf(v, 0) {
			n++
		}
	}
	return n
}
