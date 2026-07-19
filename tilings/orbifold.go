package tilings

import (
	"errors"
	"math"
)

// ErrInvalidOrbifold is returned when an orbifold symbol cannot be parsed.
var ErrInvalidOrbifold = errors.New("tilings: invalid orbifold symbol")

// OrbifoldCost returns the total cost of a Conway orbifold symbol, the sum of
// the individual feature costs used in Conway's magic theorem. A digit n before
// any '*' (a gyration point) costs (n-1)/n; a digit n after '*' (a kaleidoscope
// corner) costs (n-1)/(2n); '*' (a mirror boundary) and 'x' (a cross-cap) each
// cost 1; 'o' (a handle) costs 2. The infinity symbol '∞' is treated as the
// limit n -> infinity, costing 1 before '*' and 1/2 after. A total cost of 2
// characterises the plane groups (wallpaper and frieze). It returns an error for
// unrecognised characters.
func OrbifoldCost(symbol string) (float64, error) {
	var cost float64
	afterStar := false
	for _, r := range symbol {
		switch {
		case r == 'o' || r == 'O':
			cost += 2
		case r == 'x' || r == 'X' || r == '×':
			cost += 1
		case r == '*':
			cost += 1
			afterStar = true
		case r == '∞':
			if afterStar {
				cost += 0.5
			} else {
				cost += 1
			}
		case r >= '2' && r <= '9':
			n := float64(r - '0')
			if afterStar {
				cost += (n - 1) / (2 * n)
			} else {
				cost += (n - 1) / n
			}
		case r == '1' || r == ' ':
			// A 1-fold centre contributes nothing.
		default:
			return 0, ErrInvalidOrbifold
		}
	}
	return cost, nil
}

// OrbifoldEulerCharacteristic returns the orbifold Euler characteristic
// 2 - OrbifoldCost(symbol). It is 0 for every wallpaper and frieze group.
func OrbifoldEulerCharacteristic(symbol string) (float64, error) {
	c, err := OrbifoldCost(symbol)
	if err != nil {
		return 0, err
	}
	return 2 - c, nil
}

// IsWallpaperSignature reports whether the orbifold symbol has total cost 2 (to
// within eps), the signature of a two-dimensional space group.
func IsWallpaperSignature(symbol string, eps float64) bool {
	c, err := OrbifoldCost(symbol)
	if err != nil {
		return false
	}
	return math.Abs(c-2) <= eps
}

// GyrationOrders returns the orders of the cone (gyration) points listed before
// any '*' in the orbifold symbol; '∞' is reported as 0.
func GyrationOrders(symbol string) []int {
	var out []int
	for _, r := range symbol {
		if r == '*' {
			break
		}
		switch {
		case r == '∞':
			out = append(out, 0)
		case r >= '2' && r <= '9':
			out = append(out, int(r-'0'))
		}
	}
	return out
}

// KaleidoscopeCorners returns the orders of the corner (kaleidoscope) points
// listed after a '*' in the orbifold symbol; '∞' is reported as 0.
func KaleidoscopeCorners(symbol string) []int {
	var out []int
	afterStar := false
	for _, r := range symbol {
		if r == '*' {
			afterStar = true
			continue
		}
		if !afterStar {
			continue
		}
		switch {
		case r == '∞':
			out = append(out, 0)
		case r >= '2' && r <= '9':
			out = append(out, int(r-'0'))
		}
	}
	return out
}

// HasMirror reports whether the orbifold symbol contains a mirror boundary '*'.
func HasMirror(symbol string) bool {
	for _, r := range symbol {
		if r == '*' {
			return true
		}
	}
	return false
}

// HasCrossCap reports whether the orbifold symbol contains a cross-cap ('x').
func HasCrossCap(symbol string) bool {
	for _, r := range symbol {
		if r == 'x' || r == 'X' || r == '×' {
			return true
		}
	}
	return false
}
