package cellular

import (
	"fmt"
	"strings"
)

// RenderASCII renders a binary spacetime diagram as a multi-line string, using
// on for non-zero cells and off for zero cells, one row per line.
func RenderASCII(rows [][]int, off, on rune) string {
	var b strings.Builder
	for i, row := range rows {
		if i > 0 {
			b.WriteByte('\n')
		}
		b.WriteString(StateToString(row, off, on))
	}
	return b.String()
}

// RenderStates renders a spacetime diagram over an arbitrary alphabet using the
// runes in glyphs indexed by cell value. Cell values outside the glyph range are
// rendered with '?'.
func RenderStates(rows [][]int, glyphs []rune) string {
	var b strings.Builder
	for i, row := range rows {
		if i > 0 {
			b.WriteByte('\n')
		}
		for _, v := range row {
			if v >= 0 && v < len(glyphs) {
				b.WriteRune(glyphs[v])
			} else {
				b.WriteRune('?')
			}
		}
	}
	return b.String()
}

// TransposeRows returns the transpose of a rectangular spacetime diagram, so
// that column j of the input becomes row j of the output. It returns nil for an
// empty or ragged input.
func TransposeRows(rows [][]int) [][]int {
	if len(rows) == 0 {
		return nil
	}
	w := len(rows[0])
	for _, r := range rows {
		if len(r) != w {
			return nil
		}
	}
	out := make([][]int, w)
	for j := 0; j < w; j++ {
		col := make([]int, len(rows))
		for i := range rows {
			col[i] = rows[i][j]
		}
		out[j] = col
	}
	return out
}

// RowPopulations returns the population (non-zero count) of each row of a
// spacetime diagram.
func RowPopulations(rows [][]int) []int {
	out := make([]int, len(rows))
	for i, r := range rows {
		out[i] = Population(r)
	}
	return out
}

// RowDensities returns the density of each row of a spacetime diagram.
func RowDensities(rows [][]int) []float64 {
	out := make([]float64, len(rows))
	for i, r := range rows {
		out[i] = Density(r)
	}
	return out
}

// LightConeWidth returns, for each row of a diagram grown from a central seed,
// the distance from the centre to the outermost non-zero cell, i.e. the observed
// half-width of the light cone at that time step.
func LightConeWidth(rows [][]int) []int {
	out := make([]int, len(rows))
	for i, row := range rows {
		centre := len(row) / 2
		w := 0
		for j, v := range row {
			if v != 0 {
				d := j - centre
				if d < 0 {
					d = -d
				}
				if d > w {
					w = d
				}
			}
		}
		out[i] = w
	}
	return out
}

// DiagramString is a convenience that renders a binary diagram with '.' for 0
// and '#' for 1.
func DiagramString(rows [][]int) string {
	return RenderASCII(rows, '.', '#')
}

// ParseDiagram parses a multi-line ASCII binary diagram (as produced by
// DiagramString) back into a spacetime diagram, mapping on to 1 and every other
// rune to 0. Blank trailing lines are ignored.
func ParseDiagram(s string, on rune) [][]int {
	lines := strings.Split(s, "\n")
	var out [][]int
	for _, ln := range lines {
		if ln == "" {
			continue
		}
		out = append(out, StateFromString(ln, on))
	}
	return out
}

// CountOnes returns the total number of non-zero cells across a spacetime
// diagram.
func CountOnes(rows [][]int) int {
	total := 0
	for _, r := range rows {
		total += Population(r)
	}
	return total
}

// FinalRow returns the last row of a spacetime diagram, or nil if it is empty.
func FinalRow(rows [][]int) []int {
	if len(rows) == 0 {
		return nil
	}
	return rows[len(rows)-1]
}

// DiagramDimensions returns the number of rows and the (first-row) width of a
// diagram.
func DiagramDimensions(rows [][]int) (height, width int) {
	height = len(rows)
	if height > 0 {
		width = len(rows[0])
	}
	return height, width
}

// FormatRuleSummary returns a one-line summary of an elementary rule combining
// its number, symmetry conjugates, lambda parameter and heuristic class.
func FormatRuleSummary(r ElementaryRule) string {
	return fmt.Sprintf("rule %d: mirror=%d complement=%d lambda=%.3f class=%d",
		int(r), int(r.MirrorRule()), int(r.ComplementRule()), r.LambdaParameter(), r.Class())
}
