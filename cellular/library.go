package cellular

// mustRule builds a LifeRule from a rulestring and panics on error; it is used
// only for the compile-time-constant named rules below.
func mustRule(s string) *LifeRule {
	l, err := ParseRuleString(s)
	if err != nil {
		panic(err)
	}
	return l
}

// Conway returns Conway's Game of Life, B3/S23.
func Conway() *LifeRule { return mustRule("B3/S23") }

// HighLife returns the HighLife rule B36/S23, notable for its replicator.
func HighLife() *LifeRule { return mustRule("B36/S23") }

// DayAndNight returns the Day & Night rule B3678/S34678, which is symmetric
// under exchange of live and dead cells.
func DayAndNight() *LifeRule { return mustRule("B3678/S34678") }

// Seeds returns the Seeds rule B2/S, in which every live cell dies each
// generation.
func Seeds() *LifeRule { return mustRule("B2/S") }

// LifeWithoutDeath returns B3/S012345678, in which cells never die (Life without
// Death).
func LifeWithoutDeath() *LifeRule { return mustRule("B3/S012345678") }

// Replicator returns the Replicator rule B1357/S1357, which copies any pattern.
func Replicator() *LifeRule { return mustRule("B1357/S1357") }

// Diamoeba returns the Diamoeba rule B35678/S5678.
func Diamoeba() *LifeRule { return mustRule("B35678/S5678") }

// Maze returns the Maze rule B3/S12345.
func Maze() *LifeRule { return mustRule("B3/S12345") }

// Move returns the Move rule B368/S245, which supports slow-moving spaceships.
func Move() *LifeRule { return mustRule("B368/S245") }

// Anneal returns the Anneal (twisted majority) rule B4678/S35678.
func Anneal() *LifeRule { return mustRule("B4678/S35678") }

// Coral returns the Coral rule B3/S45678.
func Coral() *LifeRule { return mustRule("B3/S45678") }

// Amoeba returns the Amoeba rule B357/S1358.
func Amoeba() *LifeRule { return mustRule("B357/S1358") }

// mustPattern builds a Grid from ASCII rows using '#' for live cells and panics
// on error; used for the named starting patterns below.
func mustPattern(lines []string) *Grid {
	g, err := GridFromStrings(lines, '#')
	if err != nil {
		panic(err)
	}
	return g
}

// Block returns the 2x2 still life (Conway period 1).
func Block() *Grid {
	return mustPattern([]string{
		"##",
		"##",
	})
}

// Beehive returns the beehive still life.
func Beehive() *Grid {
	return mustPattern([]string{
		".##.",
		"#..#",
		".##.",
	})
}

// Blinker returns the horizontal blinker, a Conway period-2 oscillator.
func Blinker() *Grid {
	return mustPattern([]string{
		".....",
		".....",
		".###.",
		".....",
		".....",
	})
}

// Toad returns the toad, a Conway period-2 oscillator.
func Toad() *Grid {
	return mustPattern([]string{
		"......",
		"..###.",
		".###..",
		"......",
	})
}

// Beacon returns the beacon, a Conway period-2 oscillator.
func Beacon() *Grid {
	return mustPattern([]string{
		"......",
		".##...",
		".##...",
		"...##.",
		"...##.",
		"......",
	})
}

// Glider returns the glider, the smallest Conway spaceship, travelling diagonally
// with period 4.
func Glider() *Grid {
	return mustPattern([]string{
		".#.",
		"..#",
		"###",
	})
}

// LWSS returns the lightweight spaceship, a Conway period-4 orthogonal
// spaceship.
func LWSS() *Grid {
	return mustPattern([]string{
		"#..#.",
		"....#",
		"#...#",
		".####",
	})
}

// RPentomino returns the R-pentomino, a famous long-lived Conway methuselah.
func RPentomino() *Grid {
	return mustPattern([]string{
		".##",
		"##.",
		".#.",
	})
}

// Diehard returns the diehard methuselah, which vanishes after 130 Conway
// generations.
func Diehard() *Grid {
	return mustPattern([]string{
		"......#.",
		"##......",
		".#...###",
	})
}

// Acorn returns the acorn methuselah, which stabilises only after thousands of
// Conway generations.
func Acorn() *Grid {
	return mustPattern([]string{
		".#.....",
		"...#...",
		"##..###",
	})
}

// GliderOn returns a rows x cols grid containing a single glider stamped near the
// top-left corner. It returns nil if the grid is too small to hold the glider.
func GliderOn(rows, cols int) *Grid {
	if rows < 3 || cols < 3 {
		return nil
	}
	g := NewGrid(rows, cols)
	Stamp(g, Glider(), 0, 0)
	return g
}

// CenteredPattern returns a rows x cols grid with pattern stamped so its bounding
// box is centred. It returns nil if the pattern does not fit.
func CenteredPattern(pattern *Grid, rows, cols int) *Grid {
	if pattern.rows > rows || pattern.cols > cols {
		return nil
	}
	g := NewGrid(rows, cols)
	top := (rows - pattern.rows) / 2
	left := (cols - pattern.cols) / 2
	Stamp(g, pattern, top, left)
	return g
}
