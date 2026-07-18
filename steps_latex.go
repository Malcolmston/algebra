package algebra

import (
	"strconv"
	"strings"
)

// This file renders a worked [Solution] (and its individual [Step]s) as LaTeX,
// reusing the package's expression renderer [LaTeX] for every piece of
// mathematics so that the derivation looks exactly like the rest of the
// package's output. The whole solution is emitted as a single aligned
// environment.

// stepArrow is the transformation arrow placed between a step's before and
// after expressions.
const stepArrow = " \\;\\longrightarrow\\; "

// latexText escapes s so it is safe to place inside a LaTeX \text{...} group.
// The LaTeX special characters are neutralised so an explanation containing,
// for example, "^", "_" or "&" cannot break the surrounding math environment.
func latexText(s string) string {
	repl := strings.NewReplacer(
		"\\", "\\textbackslash{}",
		"{", "\\{",
		"}", "\\}",
		"$", "\\$",
		"&", "\\&",
		"%", "\\%",
		"#", "\\#",
		"_", "\\_",
		"^", "\\textasciicircum{}",
		"~", "\\textasciitilde{}",
	)
	return repl.Replace(s)
}

// stepMath renders the before (and, when it differs, after) mathematics of a
// step, joined by the transformation arrow.
func stepMath(s Step) string {
	math := LaTeX(s.Before)
	if s.After != nil && !s.Before.Equal(s.After) {
		math += stepArrow + LaTeX(s.After)
	}
	return math
}

// LaTeX renders the step as a self-contained LaTeX fragment combining its
// plain-language explanation (as a \text{...} group) with its before/after
// mathematics rendered by the package's expression renderer.
func (s Step) LaTeX() string {
	return "\\text{" + latexText(s.Explanation) + "}\\quad " + stepMath(s)
}

// LaTeX renders the whole solution as a complete LaTeX derivation inside a
// single aligned environment. Each step contributes an explanation line
// (numbered, as a \text{...} group) followed by its rendered before/after
// mathematics; the final line shows the result. The output is deterministic and
// depends only on the standard library.
func (sol *Solution) LaTeX() string {
	var b strings.Builder
	b.WriteString("\\begin{aligned}\n")
	if sol.Title != "" {
		b.WriteString("&\\text{" + latexText(sol.Title) + "} \\\\\n")
	}
	for i, s := range sol.Steps {
		label := strconv.Itoa(i+1) + ". "
		if s.Rule != "" {
			label += "(" + s.Rule + ") "
		}
		b.WriteString("&\\text{" + latexText(label+s.Explanation) + "} \\\\\n")
		b.WriteString("&\\quad " + stepMath(s) + " \\\\\n")
	}
	if sol.Result != nil {
		b.WriteString("&\\text{Result: }\\quad " + LaTeX(sol.Result) + "\n")
	}
	b.WriteString("\\end{aligned}")
	return b.String()
}

// GenerateLaTeX returns the complete LaTeX derivation of the solution sol,
// identical to sol.LaTeX(). It is the explicit package-level entry point for
// turning a worked [Solution] into a renderable LaTeX string.
func GenerateLaTeX(sol *Solution) string { return sol.LaTeX() }

// SolutionLaTeX is an alias for [GenerateLaTeX]: it returns the complete LaTeX
// derivation of the solution sol.
func SolutionLaTeX(sol *Solution) string { return sol.LaTeX() }
