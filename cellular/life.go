package cellular

import (
	"fmt"
	"sort"
	"strings"
)

// Grid is a rectangular two-dimensional binary configuration stored in
// row-major order. It is the state type for life-like automata.
type Grid struct {
	rows, cols int
	cells      []int
}

// NewGrid returns a rows x cols grid of zeros. It returns nil for non-positive
// dimensions.
func NewGrid(rows, cols int) *Grid {
	if rows <= 0 || cols <= 0 {
		return nil
	}
	return &Grid{rows: rows, cols: cols, cells: make([]int, rows*cols)}
}

// NewGridFrom builds a grid from a rectangular slice of rows, mapping non-zero
// entries to 1. It returns an error for empty or ragged input.
func NewGridFrom(cells [][]int) (*Grid, error) {
	if len(cells) == 0 || len(cells[0]) == 0 {
		return nil, fmt.Errorf("cellular: NewGridFrom needs a non-empty grid")
	}
	cols := len(cells[0])
	g := NewGrid(len(cells), cols)
	for r, row := range cells {
		if len(row) != cols {
			return nil, fmt.Errorf("cellular: NewGridFrom row %d has width %d, want %d", r, len(row), cols)
		}
		for c, v := range row {
			if v != 0 {
				g.cells[r*cols+c] = 1
			}
		}
	}
	return g, nil
}

// GridFromStrings builds a grid from a slice of equal-length strings, mapping
// each rune equal to on to 1 and all others to 0.
func GridFromStrings(lines []string, on rune) (*Grid, error) {
	if len(lines) == 0 {
		return nil, fmt.Errorf("cellular: GridFromStrings needs at least one line")
	}
	cols := len([]rune(lines[0]))
	g := NewGrid(len(lines), cols)
	for r, ln := range lines {
		rs := []rune(ln)
		if len(rs) != cols {
			return nil, fmt.Errorf("cellular: GridFromStrings line %d has width %d, want %d", r, len(rs), cols)
		}
		for c, ch := range rs {
			if ch == on {
				g.cells[r*cols+c] = 1
			}
		}
	}
	return g, nil
}

// Rows returns the number of rows.
func (g *Grid) Rows() int { return g.rows }

// Cols returns the number of columns.
func (g *Grid) Cols() int { return g.cols }

// at reads a cell without bounds checking beyond the slice; callers pass indices
// already reduced into range.
func (g *Grid) at(r, c int) int { return g.cells[r*g.cols+c] }

// set writes a cell in range.
func (g *Grid) set(r, c, v int) { g.cells[r*g.cols+c] = v }

// At returns the value of cell (r, c). Out-of-range coordinates return 0.
func (g *Grid) At(r, c int) int {
	if r < 0 || r >= g.rows || c < 0 || c >= g.cols {
		return 0
	}
	return g.at(r, c)
}

// Set writes v (any non-zero value becomes 1) into cell (r, c). Out-of-range
// coordinates are ignored.
func (g *Grid) Set(r, c, v int) {
	if r < 0 || r >= g.rows || c < 0 || c >= g.cols {
		return
	}
	if v != 0 {
		v = 1
	}
	g.set(r, c, v)
}

// Clone returns an independent copy of the grid.
func (g *Grid) Clone() *Grid {
	return &Grid{rows: g.rows, cols: g.cols, cells: append([]int(nil), g.cells...)}
}

// Population returns the number of live cells.
func (g *Grid) Population() int {
	c := 0
	for _, v := range g.cells {
		if v != 0 {
			c++
		}
	}
	return c
}

// Equal reports whether two grids have identical dimensions and contents.
func (g *Grid) Equal(other *Grid) bool {
	if other == nil || g.rows != other.rows || g.cols != other.cols {
		return false
	}
	for i := range g.cells {
		if g.cells[i] != other.cells[i] {
			return false
		}
	}
	return true
}

// Cells returns a fresh row-major two-dimensional copy of the grid contents.
func (g *Grid) Cells() [][]int {
	out := make([][]int, g.rows)
	for r := 0; r < g.rows; r++ {
		row := make([]int, g.cols)
		copy(row, g.cells[r*g.cols:(r+1)*g.cols])
		out[r] = row
	}
	return out
}

// String renders the grid with '.' for dead and '#' for live cells.
func (g *Grid) String() string {
	var b strings.Builder
	for r := 0; r < g.rows; r++ {
		if r > 0 {
			b.WriteByte('\n')
		}
		for c := 0; c < g.cols; c++ {
			if g.at(r, c) != 0 {
				b.WriteByte('#')
			} else {
				b.WriteByte('.')
			}
		}
	}
	return b.String()
}

// Neighbourhood selects the cell adjacency used by a life-like rule.
type Neighbourhood int

const (
	// Moore counts the 8 surrounding cells (king moves).
	Moore Neighbourhood = iota
	// VonNeumann counts only the 4 orthogonally adjacent cells.
	VonNeumann
)

// String names the neighbourhood.
func (n Neighbourhood) String() string {
	switch n {
	case Moore:
		return "Moore"
	case VonNeumann:
		return "von Neumann"
	default:
		return fmt.Sprintf("Neighbourhood(%d)", int(n))
	}
}

// mooreOffsets and vonNeumannOffsets list the relative neighbour coordinates.
var mooreOffsets = [][2]int{
	{-1, -1}, {-1, 0}, {-1, 1},
	{0, -1}, {0, 1},
	{1, -1}, {1, 0}, {1, 1},
}

var vonNeumannOffsets = [][2]int{
	{-1, 0}, {0, -1}, {0, 1}, {1, 0},
}

// CountNeighbours returns the number of live cells adjacent to (r, c) under the
// given neighbourhood and boundary condition. For Periodic the grid is a torus;
// FixedZero, FixedOne and Reflect apply componentwise to each axis.
func CountNeighbours(g *Grid, r, c int, nb Neighbourhood, bc Boundary) int {
	offs := mooreOffsets
	if nb == VonNeumann {
		offs = vonNeumannOffsets
	}
	count := 0
	for _, o := range offs {
		count += gridCell(g, r+o[0], c+o[1], bc)
	}
	return count
}

// gridCell reads a possibly out-of-range grid cell under a boundary condition.
func gridCell(g *Grid, r, c int, bc Boundary) int {
	rr := axisIndex(r, g.rows, bc)
	cc := axisIndex(c, g.cols, bc)
	if rr < 0 || cc < 0 {
		if bc == FixedOne {
			return 1
		}
		return 0
	}
	return g.at(rr, cc)
}

// axisIndex maps an index along one axis into range, or returns -1 to signal an
// out-of-range fixed-value read.
func axisIndex(i, n int, bc Boundary) int {
	if i >= 0 && i < n {
		return i
	}
	switch bc {
	case Periodic:
		return ((i % n) + n) % n
	case Reflect:
		if n == 1 {
			return 0
		}
		m := ((i % (2 * n)) + 2*n) % (2 * n)
		if m >= n {
			m = 2*n - 1 - m
		}
		return m
	default: // FixedZero or FixedOne handled by caller
		return -1
	}
}

// LifeRule is a life-like (outer-totalistic binary two-dimensional) rule given
// by its birth and survival neighbour counts. Born[n] is true if a dead cell
// with n live neighbours becomes alive; Survive[n] is true if a live cell with n
// live neighbours stays alive.
type LifeRule struct {
	Born    [9]bool
	Survive [9]bool
}

// NewLifeRule builds a life-like rule from explicit birth and survival counts.
// It returns an error if any count is outside 0..8.
func NewLifeRule(born, survive []int) (*LifeRule, error) {
	l := &LifeRule{}
	for _, n := range born {
		if n < 0 || n > 8 {
			return nil, fmt.Errorf("cellular: birth count %d out of range [0,8]", n)
		}
		l.Born[n] = true
	}
	for _, n := range survive {
		if n < 0 || n > 8 {
			return nil, fmt.Errorf("cellular: survival count %d out of range [0,8]", n)
		}
		l.Survive[n] = true
	}
	return l, nil
}

// ParseRuleString parses a life-like rulestring in B/S form (for example
// "B3/S23") or in the survival/birth Golly form ("23/3"). Digits must lie in
// 0..8.
func ParseRuleString(s string) (*LifeRule, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil, fmt.Errorf("cellular: empty rulestring")
	}
	var bornStr, survStr string
	if s[0] == 'B' || s[0] == 'b' {
		parts := strings.Split(s, "/")
		if len(parts) != 2 {
			return nil, fmt.Errorf("cellular: rulestring %q must have one '/'", s)
		}
		b := strings.TrimSpace(parts[0])
		sv := strings.TrimSpace(parts[1])
		if len(b) == 0 || (b[0] != 'B' && b[0] != 'b') {
			return nil, fmt.Errorf("cellular: rulestring %q missing B section", s)
		}
		if len(sv) == 0 || (sv[0] != 'S' && sv[0] != 's') {
			return nil, fmt.Errorf("cellular: rulestring %q missing S section", s)
		}
		bornStr = b[1:]
		survStr = sv[1:]
	} else {
		parts := strings.Split(s, "/")
		if len(parts) != 2 {
			return nil, fmt.Errorf("cellular: rulestring %q must have one '/'", s)
		}
		survStr = strings.TrimSpace(parts[0])
		bornStr = strings.TrimSpace(parts[1])
	}
	born, err := parseCounts(bornStr)
	if err != nil {
		return nil, err
	}
	surv, err := parseCounts(survStr)
	if err != nil {
		return nil, err
	}
	return NewLifeRule(born, surv)
}

// parseCounts converts a run of digits into a slice of counts.
func parseCounts(s string) ([]int, error) {
	var out []int
	for _, ch := range s {
		if ch < '0' || ch > '8' {
			return nil, fmt.Errorf("cellular: invalid rulestring digit %q", ch)
		}
		out = append(out, int(ch-'0'))
	}
	return out, nil
}

// String returns the rule in canonical B/S form, for example "B3/S23".
func (l *LifeRule) String() string {
	var b strings.Builder
	b.WriteByte('B')
	for n := 0; n <= 8; n++ {
		if l.Born[n] {
			b.WriteByte(byte('0' + n))
		}
	}
	b.WriteString("/S")
	for n := 0; n <= 8; n++ {
		if l.Survive[n] {
			b.WriteByte(byte('0' + n))
		}
	}
	return b.String()
}

// Next returns the next value of a cell given its current value and its live
// neighbour count.
func (l *LifeRule) Next(alive, neighbours int) int {
	if alive != 0 {
		if l.Survive[neighbours] {
			return 1
		}
		return 0
	}
	if l.Born[neighbours] {
		return 1
	}
	return 0
}

// LifeStep advances a grid one generation under a life-like rule, neighbourhood
// and boundary condition, returning a new grid. The input is not modified.
func LifeStep(l *LifeRule, g *Grid, nb Neighbourhood, bc Boundary) *Grid {
	out := NewGrid(g.rows, g.cols)
	for r := 0; r < g.rows; r++ {
		for c := 0; c < g.cols; c++ {
			n := CountNeighbours(g, r, c, nb, bc)
			out.set(r, c, l.Next(g.at(r, c), n))
		}
	}
	return out
}

// LifeEvolve runs a life-like rule for steps generations and returns steps+1
// grids beginning with a copy of g.
func LifeEvolve(l *LifeRule, g *Grid, steps int, nb Neighbourhood, bc Boundary) []*Grid {
	frames := make([]*Grid, 0, steps+1)
	frames = append(frames, g.Clone())
	cur := g
	for t := 0; t < steps; t++ {
		cur = LifeStep(l, cur, nb, bc)
		frames = append(frames, cur)
	}
	return frames
}

// GridPopulations returns the population of each grid in a sequence of frames.
func GridPopulations(frames []*Grid) []int {
	out := make([]int, len(frames))
	for i, f := range frames {
		out[i] = f.Population()
	}
	return out
}

// FindPeriod runs a life-like rule from g under a fixed boundary and reports the
// smallest period p in 1..maxPeriod for which the configuration recurs, together
// with the number of preperiod steps before the cycle is entered. It returns
// period 0 if no cycle is found within maxSteps steps. Only exact grid equality
// (no translation) is considered, so translating patterns such as gliders are
// reported with period 0.
func FindPeriod(l *LifeRule, g *Grid, nb Neighbourhood, bc Boundary, maxSteps int) (period, preperiod int) {
	seen := map[string]int{}
	cur := g.Clone()
	for t := 0; t <= maxSteps; t++ {
		key := cur.String()
		if prev, ok := seen[key]; ok {
			return t - prev, prev
		}
		seen[key] = t
		cur = LifeStep(l, cur, nb, bc)
	}
	return 0, 0
}

// LiveCells returns the sorted coordinates of the live cells of a grid as [row,
// col] pairs, ordered by row then column.
func LiveCells(g *Grid) [][2]int {
	var out [][2]int
	for r := 0; r < g.rows; r++ {
		for c := 0; c < g.cols; c++ {
			if g.at(r, c) != 0 {
				out = append(out, [2]int{r, c})
			}
		}
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i][0] != out[j][0] {
			return out[i][0] < out[j][0]
		}
		return out[i][1] < out[j][1]
	})
	return out
}

// BoundingBox returns the smallest rectangle containing every live cell as
// (minRow, minCol, maxRow, maxCol). If the grid is empty it returns all zeros
// and ok=false.
func BoundingBox(g *Grid) (minRow, minCol, maxRow, maxCol int, ok bool) {
	minRow, minCol = g.rows, g.cols
	maxRow, maxCol = -1, -1
	for r := 0; r < g.rows; r++ {
		for c := 0; c < g.cols; c++ {
			if g.at(r, c) != 0 {
				if r < minRow {
					minRow = r
				}
				if r > maxRow {
					maxRow = r
				}
				if c < minCol {
					minCol = c
				}
				if c > maxCol {
					maxCol = c
				}
				ok = true
			}
		}
	}
	if !ok {
		return 0, 0, 0, 0, false
	}
	return minRow, minCol, maxRow, maxCol, true
}

// Stamp copies the live cells of pattern into g with its top-left corner at (top,
// left), leaving the rest of g unchanged. Cells that fall outside g are ignored.
func Stamp(g, pattern *Grid, top, left int) {
	for r := 0; r < pattern.rows; r++ {
		for c := 0; c < pattern.cols; c++ {
			if pattern.at(r, c) != 0 {
				g.Set(top+r, left+c, 1)
			}
		}
	}
}
