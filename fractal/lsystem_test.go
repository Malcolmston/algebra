package fractal

import (
	"math"
	"testing"
)

func TestLSystemExpand(t *testing.T) {
	k := KochLSystem()
	if got := k.Expand(0); got != "F" {
		t.Errorf("koch expand 0: got %q want F", got)
	}
	if got := k.Expand(1); got != "F+F--F+F" {
		t.Errorf("koch expand 1: got %q", got)
	}
	// Each iteration multiplies the number of F symbols by 4.
	countF := func(s string) int {
		n := 0
		for _, r := range s {
			if r == 'F' {
				n++
			}
		}
		return n
	}
	if got := countF(k.Expand(3)); got != 64 {
		t.Errorf("koch F count at depth 3: got %d want 64", got)
	}
}

func TestLSystemDragon(t *testing.T) {
	d := DragonLSystem()
	if got := d.Expand(0); got != "FX" {
		t.Errorf("dragon expand 0: got %q", got)
	}
	// One rewrite: F stays, X -> X+YF+.
	if got := d.Expand(1); got != "FX+YF+" {
		t.Errorf("dragon expand 1: got %q want FX+YF+", got)
	}
}

func TestTurtleSquare(t *testing.T) {
	cfg := TurtleConfig{Step: 1, AngleDeg: 90, StartAngleDeg: 0, Start: Point2D{0, 0}}
	segs := Turtle("F+F+F+F", cfg)
	if len(segs) != 4 {
		t.Fatalf("expected 4 segments, got %d", len(segs))
	}
	// A unit square traced counterclockwise returns to the origin.
	last := segs[len(segs)-1]
	approx(t, last.X1, 0, 1e-12, "square close X")
	approx(t, last.Y1, 0, 1e-12, "square close Y")
	// First segment goes east.
	approx(t, segs[0].X1, 1, 1e-12, "first seg X")
	approx(t, segs[0].Y1, 0, 1e-12, "first seg Y")
	// Second segment goes north from (1,0) to (1,1).
	approx(t, segs[1].X1, 1, 1e-12, "second seg X")
	approx(t, segs[1].Y1, 1, 1e-12, "second seg Y")
}

func TestTurtleBranchingCount(t *testing.T) {
	// Bracketed branch: F[+F]F draws 3 forward moves.
	cfg := TurtleConfig{Step: 1, AngleDeg: 25, StartAngleDeg: 90}
	segs := Turtle("F[+F]F", cfg)
	if len(segs) != 3 {
		t.Fatalf("expected 3 segments, got %d", len(segs))
	}
	// After the closing bracket the turtle resumes from the pushed state:
	// the third segment starts where the first ended.
	if segs[2].X0 != segs[0].X1 || segs[2].Y0 != segs[0].Y1 {
		t.Errorf("branch pop did not restore state")
	}
}

func TestTurtlePathLengthKoch(t *testing.T) {
	// The Koch L-system rendered with a 60-degree turn has total drawn length
	// equal to (number of F symbols) * step.
	k := KochLSystem()
	cmds := k.Expand(2)
	cfg := TurtleConfig{Step: 1, AngleDeg: 60, StartAngleDeg: 0}
	segs := Turtle(cmds, cfg)
	total := 0.0
	for _, s := range segs {
		total += math.Hypot(s.X1-s.X0, s.Y1-s.Y0)
	}
	approx(t, total, 16, 1e-9, "koch depth-2 length") // 4^2 = 16 unit steps
}
