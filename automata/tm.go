package automata

import "fmt"

// Direction is the head-movement direction of a Turing machine or the tape
// motion of a two-way automaton.
type Direction int

const (
	// Left moves the head one cell towards lower indices.
	Left Direction = -1
	// Stay leaves the head in place.
	Stay Direction = 0
	// Right moves the head one cell towards higher indices.
	Right Direction = 1
)

// String renders a Direction as "L", "S" or "R".
func (d Direction) String() string {
	switch d {
	case Left:
		return "L"
	case Right:
		return "R"
	default:
		return "S"
	}
}

// TMRule is the action of a deterministic Turing-machine transition: the state
// to enter, the symbol to write, and the direction to move the head.
type TMRule struct {
	NextState string
	Write     rune
	Move      Direction
}

// tmKey identifies a (state, symbol) pair in the transition table.
type tmKey struct {
	state string
	sym   rune
}

// TM is a deterministic single-tape Turing machine. The tape is unbounded in
// both directions and filled with the Blank symbol outside the written region.
type TM struct {
	// Start is the initial control state.
	Start string
	// Accept is the accepting halting state.
	Accept string
	// Reject is the rejecting halting state (may be empty if unused).
	Reject string
	// Blank is the tape's blank symbol.
	Blank rune
	rules map[tmKey]TMRule
}

// NewTM constructs a Turing machine with the given start, accept and reject
// states and blank symbol.
func NewTM(start, accept, reject string, blank rune) *TM {
	return &TM{
		Start:  start,
		Accept: accept,
		Reject: reject,
		Blank:  blank,
		rules:  make(map[tmKey]TMRule),
	}
}

// SetRule installs the transition δ(state, read) = (next, write, move).
func (m *TM) SetRule(state string, read rune, next string, write rune, move Direction) {
	m.rules[tmKey{state, read}] = TMRule{NextState: next, Write: write, Move: move}
}

// Rule returns the transition for (state, read) and whether one is defined.
func (m *TM) Rule(state string, read rune) (TMRule, bool) {
	r, ok := m.rules[tmKey{state, read}]
	return r, ok
}

// TMConfig is an instantaneous description of a Turing-machine computation: the
// current state, the tape contents and the head position.
type TMConfig struct {
	State string
	Tape  []rune
	Head  int
	blank rune
}

// Read returns the symbol under the head, extending the tape view with blanks
// as needed (without mutating the configuration).
func (c *TMConfig) Read() rune {
	if c.Head < 0 || c.Head >= len(c.Tape) {
		return c.blank
	}
	return c.Tape[c.Head]
}

// TapeString returns the tape contents as a string with leading and trailing
// blank symbols trimmed.
func (c *TMConfig) TapeString() string {
	lo, hi := 0, len(c.Tape)
	for lo < hi && c.Tape[lo] == c.blank {
		lo++
	}
	for hi > lo && c.Tape[hi-1] == c.blank {
		hi--
	}
	return string(c.Tape[lo:hi])
}

// InitialConfig returns the starting configuration for the given input, with the
// head at the leftmost input cell.
func (m *TM) InitialConfig(input string) TMConfig {
	tape := []rune(input)
	if len(tape) == 0 {
		tape = []rune{m.Blank}
	}
	return TMConfig{State: m.Start, Tape: tape, Head: 0, blank: m.Blank}
}

// Step applies one transition to c and returns the successor configuration and
// true, or c unchanged and false if the machine has halted (no applicable rule
// or an accept/reject state).
func (m *TM) Step(c TMConfig) (TMConfig, bool) {
	if c.State == m.Accept || c.State == m.Reject {
		return c, false
	}
	read := c.Read()
	rule, ok := m.rules[tmKey{c.State, read}]
	if !ok {
		return c, false
	}
	tape := append([]rune{}, c.Tape...)
	head := c.Head
	// Grow the tape if the head is at or past an edge.
	if head < 0 {
		grow := make([]rune, -head)
		for i := range grow {
			grow[i] = m.Blank
		}
		tape = append(grow, tape...)
		head = 0
	}
	for head >= len(tape) {
		tape = append(tape, m.Blank)
	}
	tape[head] = rule.Write
	head += int(rule.Move)
	next := TMConfig{State: rule.NextState, Tape: tape, Head: head, blank: m.Blank}
	// Normalise a negative head by prepending a blank.
	if next.Head < 0 {
		next.Tape = append([]rune{m.Blank}, next.Tape...)
		next.Head = 0
	}
	return next, true
}

// TMResult reports the outcome of running a Turing machine.
type TMResult struct {
	// Accepted is true when the machine halted in the accept state.
	Accepted bool
	// Halted is true when the machine halted within the step budget.
	Halted bool
	// Steps is the number of transitions taken.
	Steps int
	// Final is the configuration at which the run stopped.
	Final TMConfig
}

// Run simulates the machine on input for up to maxSteps transitions and reports
// the outcome. A non-positive maxSteps means run until the machine halts, which
// may not terminate for some machines and inputs.
func (m *TM) Run(input string, maxSteps int) TMResult {
	c := m.InitialConfig(input)
	steps := 0
	for {
		if c.State == m.Accept {
			return TMResult{Accepted: true, Halted: true, Steps: steps, Final: c}
		}
		if c.State == m.Reject {
			return TMResult{Accepted: false, Halted: true, Steps: steps, Final: c}
		}
		if maxSteps > 0 && steps >= maxSteps {
			return TMResult{Accepted: false, Halted: false, Steps: steps, Final: c}
		}
		next, ok := m.Step(c)
		if !ok {
			// Halted with no rule and not in accept: reject by halting.
			return TMResult{Accepted: false, Halted: true, Steps: steps, Final: c}
		}
		c = next
		steps++
	}
}

// Accepts reports whether the machine accepts input within maxSteps transitions
// and returns an error if the step budget was exhausted before halting.
func (m *TM) Accepts(input string, maxSteps int) (bool, error) {
	res := m.Run(input, maxSteps)
	if !res.Halted {
		return false, fmt.Errorf("automata: TM did not halt within %d steps", maxSteps)
	}
	return res.Accepted, nil
}

// Compute runs the machine on input and returns the trimmed tape contents when
// it halts in the accept state, modelling a transducer. It errors if the machine
// rejects or fails to halt within maxSteps.
func (m *TM) Compute(input string, maxSteps int) (string, error) {
	res := m.Run(input, maxSteps)
	if !res.Halted {
		return "", fmt.Errorf("automata: TM did not halt within %d steps", maxSteps)
	}
	if !res.Accepted {
		return "", fmt.Errorf("automata: TM rejected input %q", input)
	}
	return res.Final.TapeString(), nil
}
