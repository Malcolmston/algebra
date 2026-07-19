package proofsystems

import (
	"fmt"
	"strings"
)

// HornClause is a definite Horn clause head :- b1, ..., bn, where the head and
// body atoms are positive first-order literals. A clause with an empty body is
// a fact.
type HornClause struct {
	Head FOLiteral
	Body []FOLiteral
}

// Fact builds a Horn clause with the given head and no body.
func Fact(head FOLiteral) HornClause { return HornClause{Head: head} }

// Rule builds a Horn clause head :- body.
func Rule(head FOLiteral, body ...FOLiteral) HornClause {
	cp := make([]FOLiteral, len(body))
	copy(cp, body)
	return HornClause{Head: head, Body: cp}
}

// IsFact reports whether the clause has an empty body.
func (h HornClause) IsFact() bool { return len(h.Body) == 0 }

// String renders the clause in Prolog-like syntax.
func (h HornClause) String() string {
	if h.IsFact() {
		return h.Head.String() + "."
	}
	parts := make([]string, len(h.Body))
	for i, b := range h.Body {
		parts[i] = b.String()
	}
	return h.Head.String() + " :- " + strings.Join(parts, ", ") + "."
}

// vars returns the variable names in the clause.
func (h HornClause) vars() []string {
	set := map[string]bool{}
	for _, v := range h.Head.Vars() {
		set[v] = true
	}
	for _, b := range h.Body {
		for _, v := range b.Vars() {
			set[v] = true
		}
	}
	return sortedKeys(set)
}

func (h HornClause) rename(suffix string) HornClause {
	r := map[string]string{}
	for _, v := range h.vars() {
		r[v] = v + suffix
	}
	head := renameLiteral(h.Head, r)
	body := make([]FOLiteral, len(h.Body))
	for i, b := range h.Body {
		body[i] = renameLiteral(b, r)
	}
	return HornClause{Head: head, Body: body}
}

func renameLiteral(l FOLiteral, r map[string]string) FOLiteral {
	args := make([]Term, len(l.Args))
	for i, a := range l.Args {
		args[i] = a.Rename(r)
	}
	return FOLiteral{Neg: l.Neg, Pred: l.Pred, Args: args}
}

// Program is an ordered collection of Horn clauses forming a tiny Prolog-style
// logic program. Clauses are tried top to bottom during SLD resolution.
type Program struct {
	Clauses []HornClause
}

// NewProgram builds a program from the given clauses.
func NewProgram(clauses ...HornClause) *Program {
	cp := make([]HornClause, len(clauses))
	copy(cp, clauses)
	return &Program{Clauses: cp}
}

// Add appends a clause to the program.
func (p *Program) Add(c HornClause) { p.Clauses = append(p.Clauses, c) }

// AddFact appends a fact to the program.
func (p *Program) AddFact(head FOLiteral) { p.Add(Fact(head)) }

// AddRule appends a rule to the program.
func (p *Program) AddRule(head FOLiteral, body ...FOLiteral) { p.Add(Rule(head, body...)) }

// String renders the whole program, one clause per line.
func (p *Program) String() string {
	parts := make([]string, len(p.Clauses))
	for i, c := range p.Clauses {
		parts[i] = c.String()
	}
	return strings.Join(parts, "\n")
}

// SolveOptions bounds an SLD search.
type SolveOptions struct {
	// MaxDepth limits the SLD derivation length to avoid non-termination.
	MaxDepth int
	// MaxSolutions limits how many answer substitutions are collected; zero or
	// negative means unbounded (subject to MaxDepth).
	MaxSolutions int
}

// DefaultSolveOptions returns reasonable search bounds.
func DefaultSolveOptions() SolveOptions {
	return SolveOptions{MaxDepth: 1000, MaxSolutions: 100}
}

// Solve performs SLD resolution of the goal list against the program and returns
// the answer substitutions restricted to the variables of the query. The
// clauses are tried in program order with depth-first backtracking; fresh
// variants of each clause are used so no variable capture occurs.
func (p *Program) Solve(goals []FOLiteral, opts SolveOptions) []Substitution {
	if opts.MaxDepth <= 0 {
		opts.MaxDepth = 1000
	}
	qVars := queryVars(goals)
	var answers []Substitution
	counter := 0
	var step func(gs []FOLiteral, sub Substitution, depth int)
	step = func(gs []FOLiteral, sub Substitution, depth int) {
		if opts.MaxSolutions > 0 && len(answers) >= opts.MaxSolutions {
			return
		}
		if len(gs) == 0 {
			answers = append(answers, sub.Restrict(qVars))
			return
		}
		if depth >= opts.MaxDepth {
			return
		}
		goal := sub.ApplyLiteral(gs[0])
		for _, clause := range p.Clauses {
			counter++
			fresh := clause.rename(fmt.Sprintf("_%d", counter))
			s, err := UnifyLiterals(goal, fresh.Head)
			if err != nil {
				continue
			}
			newGoals := append(append([]FOLiteral{}, fresh.Body...), gs[1:]...)
			step(newGoals, sub.Compose(s), depth+1)
			if opts.MaxSolutions > 0 && len(answers) >= opts.MaxSolutions {
				return
			}
		}
	}
	step(goals, NewSubstitution(), 0)
	return answers
}

// Query reports whether the goal list is derivable from the program (there is
// at least one SLD refutation) within the default search bounds.
func (p *Program) Query(goals ...FOLiteral) bool {
	return len(p.Solve(goals, SolveOptions{MaxDepth: 1000, MaxSolutions: 1})) > 0
}

// SolveOne returns the first answer substitution for the goal list and whether
// one was found.
func (p *Program) SolveOne(goals ...FOLiteral) (Substitution, bool) {
	ans := p.Solve(goals, SolveOptions{MaxDepth: 1000, MaxSolutions: 1})
	if len(ans) == 0 {
		return Substitution{}, false
	}
	return ans[0], true
}

// Prove is a synonym for Query used when the goals encode a yes/no question.
func (p *Program) Prove(goals ...FOLiteral) bool { return p.Query(goals...) }

func queryVars(goals []FOLiteral) []string {
	set := map[string]bool{}
	for _, g := range goals {
		for _, v := range g.Vars() {
			set[v] = true
		}
	}
	return sortedKeys(set)
}

// SLDResolve is a package-level convenience that builds a program from clauses
// and solves the goals with the default options.
func SLDResolve(clauses []HornClause, goals []FOLiteral) []Substitution {
	return NewProgram(clauses...).Solve(goals, DefaultSolveOptions())
}
