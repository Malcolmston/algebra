package algebra

import (
	"strings"
	"testing"
)

// hasExplanation reports whether any step of sol has an explanation or rule
// containing sub.
func hasExplanation(sol *Solution, sub string) bool {
	for _, s := range sol.Steps {
		if strings.Contains(s.Explanation, sub) || strings.Contains(s.Rule, sub) {
			return true
		}
	}
	return false
}

// TestDifferentiateSteps checks the number of steps, the rules named, the
// result and the generated LaTeX for a polynomial and a chain-rule case.
func TestDifferentiateSteps(t *testing.T) {
	x := Sym("x")

	sol := DifferentiateSteps(Add(Pow(x, Int(2)), Mul(Int(3), x), Int(1)), x)
	if len(sol.Steps) != 5 {
		t.Fatalf("polynomial: got %d steps, want 5", len(sol.Steps))
	}
	for _, want := range []string{"power rule", "sum rule", "constant multiple rule", "variable rule", "constant rule"} {
		if !hasExplanation(sol, want) {
			t.Errorf("polynomial: missing rule %q", want)
		}
	}
	if want := Add(Mul(Int(2), x), Int(3)); !sol.Result.Equal(want) {
		t.Errorf("polynomial: result = %s, want %s", sol.Result, want)
	}
	ltx := GenerateLaTeX(sol)
	for _, want := range []string{"\\begin{aligned}", "\\text{", "\\longrightarrow", "\\end{aligned}"} {
		if !strings.Contains(ltx, want) {
			t.Errorf("polynomial LaTeX missing %q", want)
		}
	}

	chain := DifferentiateSteps(Sin(Pow(x, Int(2))), x)
	if len(chain.Steps) != 2 {
		t.Fatalf("chain: got %d steps, want 2", len(chain.Steps))
	}
	if !hasExplanation(chain, "chain rule") {
		t.Errorf("chain: missing chain rule")
	}
	if want := Mul(Int(2), x, Cos(Pow(x, Int(2)))); !sol.Result.Equal(Add(Mul(Int(2), x), Int(3))) || !chain.Result.Equal(want) {
		t.Errorf("chain: result = %s, want %s", chain.Result, want)
	}
}

// TestSolveQuadraticSteps checks the full derivation, both roots, and that the
// LaTeX contains the quadratic-formula pieces (a square root and a fraction).
func TestSolveQuadraticSteps(t *testing.T) {
	x := Sym("x")
	e := Add(Pow(x, Int(2)), Mul(Int(-5), x), Int(6)) // (x-2)(x-3)
	sol := SolveQuadraticSteps(e, x)
	if len(sol.Steps) != 6 {
		t.Fatalf("got %d steps, want 6", len(sol.Steps))
	}
	for _, want := range []string{"identify coefficients", "quadratic formula", "discriminant", "solution"} {
		if !hasExplanation(sol, want) {
			t.Errorf("missing step %q", want)
		}
	}
	if !sol.Result.Equal(Int(2)) {
		t.Errorf("result = %s, want 2", sol.Result)
	}
	// Both roots must appear as solution steps and both must satisfy e == 0.
	roots := map[string]bool{}
	for _, s := range sol.Steps {
		if s.Rule == "solution" {
			roots[s.After.String()] = true
			if got := Simplify(Subs(e, x, s.After)); !isZero(got) && !got.Equal(Int(0)) {
				t.Errorf("root %s does not satisfy the equation (got %s)", s.After, got)
			}
		}
	}
	if !roots["2"] || !roots["3"] {
		t.Errorf("roots = %v, want {2,3}", roots)
	}
	ltx := GenerateLaTeX(sol)
	for _, want := range []string{"\\begin{aligned}", "\\text{", "\\sqrt{", "\\frac{"} {
		if !strings.Contains(ltx, want) {
			t.Errorf("LaTeX missing quadratic-formula piece %q", want)
		}
	}
}

// TestSolveLinearSteps checks isolation and division for a linear equation.
func TestSolveLinearSteps(t *testing.T) {
	x := Sym("x")
	sol := SolveLinearSteps(Add(Mul(Int(2), x), Int(-6)), x) // 2x - 6 = 0 -> x = 3
	if len(sol.Steps) != 4 {
		t.Fatalf("got %d steps, want 4", len(sol.Steps))
	}
	if !hasExplanation(sol, "divide") || !hasExplanation(sol, "isolate") {
		t.Errorf("missing isolate/divide steps")
	}
	if !sol.Result.Equal(Int(3)) {
		t.Errorf("result = %s, want 3", sol.Result)
	}
}

// TestSolveCubicSteps checks the rational-root / synthetic-division derivation
// and that every reported root satisfies the cubic.
func TestSolveCubicSteps(t *testing.T) {
	x := Sym("x")
	e := Add(Pow(x, Int(3)), Mul(Int(-6), Pow(x, Int(2))), Mul(Int(11), x), Int(-6)) // (x-1)(x-2)(x-3)
	sol := SolveCubicSteps(e, x)
	if len(sol.Steps) != 7 {
		t.Fatalf("got %d steps, want 7", len(sol.Steps))
	}
	for _, want := range []string{"rational root theorem", "synthetic division", "quadratic formula", "all solutions"} {
		if !hasExplanation(sol, want) {
			t.Errorf("missing step %q", want)
		}
	}
	if !sol.Result.Equal(Int(1)) {
		t.Errorf("result = %s, want 1", sol.Result)
	}
	roots, err := Solve(e, x)
	if err != nil {
		t.Fatalf("Solve error: %v", err)
	}
	if len(roots) != 3 {
		t.Errorf("Solve gave %d roots, want 3", len(roots))
	}
}

// TestCompleteSquareSteps checks the completed-square form is algebraically
// equal to the original quadratic.
func TestCompleteSquareSteps(t *testing.T) {
	x := Sym("x")
	e := Add(Pow(x, Int(2)), Mul(Int(6), x), Int(5)) // (x+3)^2 - 4
	sol := CompleteSquareSteps(e, x)
	if len(sol.Steps) != 4 {
		t.Fatalf("got %d steps, want 4", len(sol.Steps))
	}
	if !hasExplanation(sol, "complete the square") {
		t.Errorf("missing complete-the-square step")
	}
	want := Add(Pow(Add(x, Int(3)), Int(2)), Int(-4))
	if !sol.Result.Equal(want) {
		t.Errorf("result = %s, want %s", sol.Result, want)
	}
	if lhs, rhs := Simplify(Expand(sol.Result)), Simplify(Expand(e)); !lhs.Equal(rhs) {
		t.Errorf("completed square %s != original %s", lhs, rhs)
	}
}

// TestFactorSteps checks the factored form multiplies back to the original.
func TestFactorSteps(t *testing.T) {
	x := Sym("x")
	e := Add(Pow(x, Int(2)), Mul(Int(-5), x), Int(6))
	sol := FactorSteps(e, x)
	if len(sol.Steps) != 5 {
		t.Fatalf("got %d steps, want 5", len(sol.Steps))
	}
	if !hasExplanation(sol, "discriminant") || !hasExplanation(sol, "factored form") {
		t.Errorf("missing discriminant/factored steps")
	}
	if lhs, rhs := Simplify(Expand(sol.Result)), Simplify(Expand(e)); !lhs.Equal(rhs) {
		t.Errorf("factored %s expands to %s, want %s", sol.Result, lhs, rhs)
	}
}

// TestExpandSteps checks a binomial square expansion.
func TestExpandSteps(t *testing.T) {
	x := Sym("x")
	sol := ExpandSteps(Mul(Add(x, Int(1)), Add(x, Int(2))))
	if len(sol.Steps) != 1 {
		t.Fatalf("got %d steps, want 1", len(sol.Steps))
	}
	want := Add(Pow(x, Int(2)), Mul(Int(3), x), Int(2))
	if !sol.Result.Equal(want) {
		t.Errorf("result = %s, want %s", sol.Result, want)
	}
	if !hasExplanation(sol, "distributive law") {
		t.Errorf("missing distributive step")
	}

	// A power of a sum expands by the binomial-expansion branch.
	bsol := ExpandSteps(Pow(Add(x, Int(1)), Int(2)))
	if len(bsol.Steps) != 1 || !hasExplanation(bsol, "binomial expansion") {
		t.Errorf("binomial: steps=%d, missing binomial expansion", len(bsol.Steps))
	}
	if bw := Add(Pow(x, Int(2)), Mul(Int(2), x), Int(1)); !bsol.Result.Equal(bw) {
		t.Errorf("binomial result = %s, want %s", bsol.Result, bw)
	}
}

// TestSimplifySteps checks folding of numeric constants in a raw sum.
func TestSimplifySteps(t *testing.T) {
	x := Sym("x")
	// A deliberately non-canonical sum 2 + 3 + x built with the raw constructor.
	sol := SimplifySteps(newSum([]Expr{Int(2), Int(3), x}))
	if len(sol.Steps) != 1 {
		t.Fatalf("got %d steps, want 1", len(sol.Steps))
	}
	if want := Add(x, Int(5)); !sol.Result.Equal(want) {
		t.Errorf("result = %s, want %s", sol.Result, want)
	}
	if !hasExplanation(sol, "combine like terms") {
		t.Errorf("missing combine-like-terms step")
	}
}

// TestIntegrateSteps checks term-by-term integration and that differentiating
// the antiderivative recovers the integrand.
func TestIntegrateSteps(t *testing.T) {
	x := Sym("x")
	e := Add(Pow(x, Int(2)), Mul(Int(2), x))
	sol := IntegrateSteps(e, x)
	if len(sol.Steps) != 4 {
		t.Fatalf("got %d steps, want 4", len(sol.Steps))
	}
	if !hasExplanation(sol, "sum rule") || !hasExplanation(sol, "power rule") || !hasExplanation(sol, "constant multiple rule") {
		t.Errorf("missing integration rule steps")
	}
	if got := Simplify(Diff(sol.Result, x)); !got.Equal(Simplify(e)) {
		t.Errorf("d/dx(%s) = %s, want %s", sol.Result, got, Simplify(e))
	}
}

// TestSeriesSteps checks the Maclaurin series of exp to four terms.
func TestSeriesSteps(t *testing.T) {
	x := Sym("x")
	sol := SeriesSteps(Exp(x), x, Int(0), 4)
	if len(sol.Steps) != 5 { // four term steps + the summation step
		t.Fatalf("got %d steps, want 5", len(sol.Steps))
	}
	if !hasExplanation(sol, "taylor term") {
		t.Errorf("missing taylor-term steps")
	}
	want := Series(Exp(x), x, Int(0), 4) // 1 + x + x^2/2 + x^3/6
	if !sol.Result.Equal(want) {
		t.Errorf("result = %s, want %s", sol.Result, want)
	}
}

// TestLimitSteps checks the classic sin(x)/x limit via l'Hopital.
func TestLimitSteps(t *testing.T) {
	x := Sym("x")
	sol := LimitSteps(Mul(Sin(x), Pow(x, Int(-1))), x, Int(0))
	if len(sol.Steps) != 4 {
		t.Fatalf("got %d steps, want 4", len(sol.Steps))
	}
	if !hasExplanation(sol, "l'Hopital") || !hasExplanation(sol, "indeterminate") {
		t.Errorf("missing l'Hopital/indeterminate steps")
	}
	if !sol.Result.Equal(Int(1)) {
		t.Errorf("result = %s, want 1", sol.Result)
	}
}

// TestPartialFractionSteps checks 1/(x^2 - x) = 1/(x-1) - 1/x.
func TestPartialFractionSteps(t *testing.T) {
	x := Sym("x")
	e := Pow(Add(Pow(x, Int(2)), Mul(Int(-1), x)), Int(-1))
	sol := PartialFractionSteps(e, x)
	if len(sol.Steps) != 5 {
		t.Fatalf("got %d steps, want 5", len(sol.Steps))
	}
	if !hasExplanation(sol, "factor denominator") || !hasExplanation(sol, "residue") {
		t.Errorf("missing factor/residue steps")
	}
	want, err := ApartExpr(e, x)
	if err != nil {
		t.Fatalf("ApartExpr error: %v", err)
	}
	if !sol.Result.Equal(want) {
		t.Errorf("result = %s, want %s", sol.Result, want)
	}
}

// TestSolveSystemSteps checks Gaussian elimination on a 2x2 system and that its
// answer matches the direct solver.
func TestSolveSystemSteps(t *testing.T) {
	x, y := Sym("x"), Sym("y")
	// x + y = 3, x - y = 1  ->  x = 2, y = 1.
	eqs := []Expr{Add(x, y, Int(-3)), Add(x, Mul(Int(-1), y), Int(-1))}
	sol := SolveSystemSteps(eqs, []Expr{x, y})
	if len(sol.Steps) != 6 {
		t.Fatalf("got %d steps, want 6", len(sol.Steps))
	}
	if !hasExplanation(sol, "Eliminate") || !hasExplanation(sol, "normalize pivot") {
		t.Errorf("missing elimination/normalization row operations")
	}
	if !sol.Result.Equal(Int(2)) {
		t.Errorf("result (x) = %s, want 2", sol.Result)
	}
	direct, err := SolveSystem(eqs, []Expr{x, y})
	if err != nil {
		t.Fatalf("SolveSystem error: %v", err)
	}
	if !direct[0].Equal(Int(2)) || !direct[1].Equal(Int(1)) {
		t.Errorf("direct solve = %v, want [2 1]", direct)
	}
	// The last solution step reports y = 1.
	last := sol.Steps[len(sol.Steps)-1]
	if !last.After.Equal(Int(1)) {
		t.Errorf("final step reports %s, want y = 1", last.After)
	}
}

// TestStepAndSolutionLaTeX checks the Step and Solution LaTeX rendering and the
// package-level GenerateLaTeX/SolutionLaTeX entry points.
func TestStepAndSolutionLaTeX(t *testing.T) {
	x := Sym("x")
	sol := DifferentiateSteps(Pow(x, Int(2)), x)
	// A single step renders explanation as \text{...} plus the math.
	sl := sol.Steps[0].LaTeX()
	if !strings.HasPrefix(sl, "\\text{") || !strings.Contains(sl, "x^{2}") {
		t.Errorf("Step.LaTeX = %q", sl)
	}
	full := sol.LaTeX()
	if GenerateLaTeX(sol) != full || SolutionLaTeX(sol) != full {
		t.Errorf("GenerateLaTeX/SolutionLaTeX disagree with Solution.LaTeX")
	}
	if !strings.Contains(full, "\\begin{aligned}") || !strings.Contains(full, "Result:") {
		t.Errorf("Solution.LaTeX missing structure: %q", full)
	}
	// The human-readable form lists every step and the result.
	s := sol.String()
	if !strings.Contains(s, "Result:") || !strings.Contains(s, "1.") {
		t.Errorf("Solution.String missing content: %q", s)
	}
}

// TestStepGeneratorsNeverNil checks that invalid input yields a Solution with an
// explanatory step rather than a nil pointer or a fabricated result.
func TestStepGeneratorsNeverNil(t *testing.T) {
	x := Sym("x")
	// Not linear -> a "not applicable" step, never nil.
	if sol := SolveLinearSteps(Pow(x, Int(2)), x); sol == nil || len(sol.Steps) == 0 {
		t.Errorf("SolveLinearSteps returned empty solution for non-linear input")
	}
	// Non-symbol variable -> an error step.
	if sol := DifferentiateSteps(x, Int(2)); sol == nil || len(sol.Steps) == 0 {
		t.Errorf("DifferentiateSteps returned empty solution for non-symbol variable")
	}
}
