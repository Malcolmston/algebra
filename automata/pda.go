package automata

// AcceptMode selects a pushdown automaton's acceptance condition.
type AcceptMode int

const (
	// AcceptByFinalState accepts when the whole input is consumed in an
	// accepting state.
	AcceptByFinalState AcceptMode = iota
	// AcceptByEmptyStack accepts when the whole input is consumed with an empty
	// stack.
	AcceptByEmptyStack
)

// PDATransition is one nondeterministic move of a pushdown automaton: enter
// NextState and replace the popped stack top with the symbols in Push, where
// Push[0] becomes the new top. An empty Push pops without pushing.
type PDATransition struct {
	NextState string
	Push      []rune
}

// pdaKey identifies a (state, input, stack-top) triple. When inputEps is true
// the move is an epsilon move that consumes no input.
type pdaKey struct {
	state    string
	input    rune
	inputEps bool
	top      rune
}

// PDA is a nondeterministic pushdown automaton with epsilon moves over rune
// input and rune stack alphabets.
type PDA struct {
	// Start is the initial control state.
	Start string
	// InitialStack is the symbol on the stack before any input is read.
	InitialStack rune
	// Accept is the set of accepting states (used in AcceptByFinalState mode).
	Accept map[string]bool
	// Mode selects the acceptance condition.
	Mode AcceptMode
	// MaxSteps bounds the number of configurations explored during a run to
	// guarantee termination in the presence of epsilon-push loops; zero selects
	// a generous default.
	MaxSteps int
	rules    map[pdaKey][]PDATransition
}

// NewPDA constructs a pushdown automaton with the given start state, initial
// stack symbol and acceptance mode.
func NewPDA(start string, initialStack rune, mode AcceptMode) *PDA {
	return &PDA{
		Start:        start,
		InitialStack: initialStack,
		Accept:       make(map[string]bool),
		Mode:         mode,
		rules:        make(map[pdaKey][]PDATransition),
	}
}

// AddAccept marks each of the given states as accepting.
func (p *PDA) AddAccept(states ...string) {
	for _, s := range states {
		p.Accept[s] = true
	}
}

// AddTransition adds a move on input symbol in (popping top): from state, read
// in and pop top, go to next and push the given symbols.
func (p *PDA) AddTransition(state string, in rune, top rune, next string, push ...rune) {
	k := pdaKey{state: state, input: in, top: top}
	p.rules[k] = append(p.rules[k], PDATransition{NextState: next, Push: append([]rune{}, push...)})
}

// AddEpsilonTransition adds an epsilon move (consuming no input) from state,
// popping top, entering next and pushing the given symbols.
func (p *PDA) AddEpsilonTransition(state string, top rune, next string, push ...rune) {
	k := pdaKey{state: state, inputEps: true, top: top}
	p.rules[k] = append(p.rules[k], PDATransition{NextState: next, Push: append([]rune{}, push...)})
}

// pdaConfig is an instantaneous description of a PDA computation.
type pdaConfig struct {
	state string
	pos   int    // index into the input runes
	stack string // stack contents, top at the end
}

// applyPush returns the stack (top at end) obtained by popping the current top
// and pushing push with push[0] on top.
func applyPush(stackWithoutTop string, push []rune) string {
	// push[0] must end on top, so append in reverse order.
	buf := []rune(stackWithoutTop)
	for i := len(push) - 1; i >= 0; i-- {
		buf = append(buf, push[i])
	}
	return string(buf)
}

// Accepts reports whether the PDA accepts input under its acceptance mode. The
// search explores the configuration graph breadth-first with a visited set and
// a step budget (PDA.MaxSteps or a default derived from the input length) to
// ensure termination.
func (p *PDA) Accepts(input string) bool {
	runes := []rune(input)
	budget := p.MaxSteps
	if budget <= 0 {
		budget = 10000 + 200*(len(runes)+1)
	}
	start := pdaConfig{state: p.Start, pos: 0, stack: string([]rune{p.InitialStack})}
	visited := map[pdaConfig]bool{start: true}
	queue := []pdaConfig{start}
	steps := 0
	for len(queue) > 0 {
		if steps >= budget {
			return false
		}
		cur := queue[0]
		queue = queue[1:]
		steps++

		if p.isAccepting(cur, len(runes)) {
			return true
		}

		stackRunes := []rune(cur.stack)
		if len(stackRunes) == 0 {
			// No top to pop; no moves possible.
			continue
		}
		top := stackRunes[len(stackRunes)-1]
		base := string(stackRunes[:len(stackRunes)-1])

		// Epsilon moves.
		for _, tr := range p.rules[pdaKey{state: cur.state, inputEps: true, top: top}] {
			nc := pdaConfig{state: tr.NextState, pos: cur.pos, stack: applyPush(base, tr.Push)}
			if !visited[nc] {
				visited[nc] = true
				queue = append(queue, nc)
			}
		}
		// Input-consuming moves.
		if cur.pos < len(runes) {
			in := runes[cur.pos]
			for _, tr := range p.rules[pdaKey{state: cur.state, input: in, top: top}] {
				nc := pdaConfig{state: tr.NextState, pos: cur.pos + 1, stack: applyPush(base, tr.Push)}
				if !visited[nc] {
					visited[nc] = true
					queue = append(queue, nc)
				}
			}
		}
	}
	return false
}

// isAccepting reports whether cur is an accepting configuration: all input read
// and either an accepting state or an empty stack depending on the mode.
func (p *PDA) isAccepting(cur pdaConfig, inputLen int) bool {
	if cur.pos != inputLen {
		return false
	}
	switch p.Mode {
	case AcceptByEmptyStack:
		return len([]rune(cur.stack)) == 0
	default:
		return p.Accept[cur.state]
	}
}
