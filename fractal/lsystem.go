package fractal

import (
	"math"
	"strings"
)

// LSystem is a deterministic context-free Lindenmayer system: an Axiom string
// and a set of production Rules mapping a single symbol (rune) to a replacement
// string. Symbols with no rule are left unchanged (they act as constants).
type LSystem struct {
	Axiom string
	Rules map[rune]string
}

// NewLSystem constructs an LSystem from an axiom and a rule map. The rule map is
// copied so later mutation of the argument does not affect the system.
func NewLSystem(axiom string, rules map[rune]string) LSystem {
	cp := make(map[rune]string, len(rules))
	for k, v := range rules {
		cp[k] = v
	}
	return LSystem{Axiom: axiom, Rules: cp}
}

// Step applies every production rule once to s in parallel, replacing each
// symbol by its rule's right-hand side (or leaving it unchanged if it has no
// rule) and returning the rewritten string.
func (l LSystem) Step(s string) string {
	var b strings.Builder
	for _, r := range s {
		if rep, ok := l.Rules[r]; ok {
			b.WriteString(rep)
		} else {
			b.WriteRune(r)
		}
	}
	return b.String()
}

// Expand applies the production rules iterations times starting from the axiom
// and returns the resulting string. Zero iterations returns the axiom itself.
// It panics if iterations is negative.
func (l LSystem) Expand(iterations int) string {
	if iterations < 0 {
		panic("fractal: Expand needs non-negative iterations")
	}
	s := l.Axiom
	for i := 0; i < iterations; i++ {
		s = l.Step(s)
	}
	return s
}

// TurtleConfig parameterizes the turtle-graphics interpreter. Step is the
// forward move length; AngleDeg is the turn angle in degrees applied by the
// '+' and '-' commands; StartAngleDeg is the initial heading in degrees
// measured counterclockwise from the positive x-axis; Start is the initial pen
// position.
type TurtleConfig struct {
	Step          float64
	AngleDeg      float64
	StartAngleDeg float64
	Start         Point2D
}

// TurtleState is the pose of the turtle: its position and its heading in
// radians measured counterclockwise from the positive x-axis.
type TurtleState struct {
	Pos     Point2D
	HeadDeg float64
}

// Segment is a directed line segment from (X0,Y0) to (X1,Y1), the drawing
// primitive produced by the turtle interpreter.
type Segment struct {
	X0, Y0, X1, Y1 float64
}

// Turtle interprets a command string as turtle graphics and returns the drawn
// line segments. The recognized commands are:
//
//	F, G : move forward by Step, drawing a segment
//	f, g : move forward by Step without drawing
//	+    : turn left  (counterclockwise) by AngleDeg
//	-    : turn right (clockwise)        by AngleDeg
//	|    : turn around (180 degrees)
//	[    : push the current state onto a stack
//	]    : pop and restore the most recently pushed state
//
// All other symbols are ignored, so L-system variables such as 'X' and 'Y' act
// as no-ops. Bracketed sections produce disconnected segment groups, allowing
// branching structures.
func Turtle(commands string, cfg TurtleConfig) []Segment {
	var segs []Segment
	pos := cfg.Start
	head := cfg.StartAngleDeg * math.Pi / 180
	turn := cfg.AngleDeg * math.Pi / 180
	var stack []TurtleState
	for _, r := range commands {
		switch r {
		case 'F', 'G':
			sn, cs := math.Sincos(head)
			np := Point2D{pos.X + cfg.Step*cs, pos.Y + cfg.Step*sn}
			segs = append(segs, Segment{pos.X, pos.Y, np.X, np.Y})
			pos = np
		case 'f', 'g':
			sn, cs := math.Sincos(head)
			pos = Point2D{pos.X + cfg.Step*cs, pos.Y + cfg.Step*sn}
		case '+':
			head += turn
		case '-':
			head -= turn
		case '|':
			head += math.Pi
		case '[':
			stack = append(stack, TurtleState{pos, head * 180 / math.Pi})
		case ']':
			if len(stack) > 0 {
				top := stack[len(stack)-1]
				stack = stack[:len(stack)-1]
				pos = top.Pos
				head = top.HeadDeg * math.Pi / 180
			}
		}
	}
	return segs
}

// TurtlePath interprets a command string and returns the ordered vertices
// visited by the turtle: the start position followed by the endpoint of every
// forward move (both drawing and non-drawing). It is intended for non-branching
// L-systems; when brackets are present the returned polyline includes the jumps
// created by popping the stack. See [Turtle] for the command set.
func TurtlePath(commands string, cfg TurtleConfig) []Point2D {
	pts := []Point2D{cfg.Start}
	pos := cfg.Start
	head := cfg.StartAngleDeg * math.Pi / 180
	turn := cfg.AngleDeg * math.Pi / 180
	var stack []TurtleState
	for _, r := range commands {
		switch r {
		case 'F', 'G', 'f', 'g':
			sn, cs := math.Sincos(head)
			pos = Point2D{pos.X + cfg.Step*cs, pos.Y + cfg.Step*sn}
			pts = append(pts, pos)
		case '+':
			head += turn
		case '-':
			head -= turn
		case '|':
			head += math.Pi
		case '[':
			stack = append(stack, TurtleState{pos, head * 180 / math.Pi})
		case ']':
			if len(stack) > 0 {
				top := stack[len(stack)-1]
				stack = stack[:len(stack)-1]
				pos = top.Pos
				head = top.HeadDeg * math.Pi / 180
				pts = append(pts, pos)
			}
		}
	}
	return pts
}

// KochLSystem returns the classic Koch curve L-system: axiom "F" with the rule
// F -> F+F--F+F. Rendered with a turn angle of 60 degrees it produces the Koch
// curve.
func KochLSystem() LSystem {
	return NewLSystem("F", map[rune]string{'F': "F+F--F+F"})
}

// SierpinskiLSystem returns an L-system whose turtle rendering (turn angle 60
// degrees) draws the Sierpinski triangle. Axiom "F-G-G" with rules
// F -> F-G+F+G-F and G -> GG.
func SierpinskiLSystem() LSystem {
	return NewLSystem("F-G-G", map[rune]string{
		'F': "F-G+F+G-F",
		'G': "GG",
	})
}

// DragonLSystem returns the Heighway dragon curve L-system: axiom "F" with
// rules X -> X+YF+ and Y -> -FX-Y, rendered with a turn angle of 90 degrees.
// The symbols X and Y control the recursion and draw nothing themselves.
func DragonLSystem() LSystem {
	return NewLSystem("FX", map[rune]string{
		'X': "X+YF+",
		'Y': "-FX-Y",
	})
}

// PlantLSystem returns a branching plant L-system: axiom "X" with rules
// X -> F+[[X]-X]-F[-FX]+X and F -> FF, rendered with a turn angle of about 25
// degrees. It uses the bracket commands to create branches.
func PlantLSystem() LSystem {
	return NewLSystem("X", map[rune]string{
		'X': "F+[[X]-X]-F[-FX]+X",
		'F': "FF",
	})
}
