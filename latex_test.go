package algebra

import "testing"

// TestLaTeX exercises the LaTeX renderer with known-answer symbolic cases.
func TestLaTeX(t *testing.T) {
	cases := []struct {
		name string
		in   Expr
		want string
	}{
		{"integer", Int(2), "2"},
		{"rational", Rat(1, 2), "\\frac{1}{2}"},
		{"neg-rational", Rat(-3, 4), "-\\frac{3}{4}"},
		{"symbol", Sym("x"), "x"},
		{"pi", Pi, "\\pi"},
		{"e", E, "e"},
		{"imag", I, "i"},
		{"inf", Inf, "\\infty"},
		{"power", Pow(Sym("x"), Int(2)), "x^{2}"},
		{"reciprocal", Pow(Sym("x"), Int(-1)), "\\frac{1}{x}"},
		{"neg-power", Pow(Sym("x"), Int(-2)), "\\frac{1}{x^{2}}"},
		{"sqrt-fn", Sqrt(Sym("x")), "\\sqrt{x}"},
		{"half-power", Pow(Sym("x"), Rat(1, 2)), "\\sqrt{x}"},
		{"sum", Add(Sym("x"), Int(1)), "x + 1"},
		{"poly", Add(Pow(Sym("x"), Int(2)), Int(1)), "x^{2} + 1"},
		{"coeff-symbol", Mul(Int(2), Sym("x")), "2 x"},
		{"neg-symbol", Mul(Int(-1), Sym("x")), "-x"},
		{"rational-coeff", Mul(Rat(1, 2), Sym("x")), "\\frac{x}{2}"},
		{"sin", Sin(Sym("x")), "\\sin\\left(x\\right)"},
		{"gamma", Gamma(Sym("x")), "\\Gamma\\left(x\\right)"},
		{"abs", Abs(Sym("x")), "\\left|x\\right|"},
		{"factorial", Factorial(Sym("x")), "x!"},
		{"integral", newIntegral(Pow(Sym("x"), Int(2)), Sym("x")), "\\int x^{2} \\, d x"},
		{"sum-op", newBigOp("Sum", Sym("k"), "k", Int(1), Sym("n")), "\\sum_{k=1}^{n} k"},
		{"limit", newBigOp("Limit", Pow(Sym("x"), Int(2)), "x", Int(0), Int(0)), "\\lim_{x \\to 0} x^{2}"},
	}
	for _, c := range cases {
		if got := LaTeX(c.in); got != c.want {
			t.Errorf("%s: LaTeX = %q, want %q", c.name, got, c.want)
		}
	}
}

// TestLaTeXEq checks the two-sided equation renderer.
func TestLaTeXEq(t *testing.T) {
	got := LaTeXEq(Sym("y"), Add(Sym("x"), Int(1)))
	if want := "y = x + 1"; got != want {
		t.Errorf("LaTeXEq = %q, want %q", got, want)
	}
}

// TestPretty exercises the 2D Unicode pretty printer with hand-verified
// layouts.
func TestPretty(t *testing.T) {
	cases := []struct {
		name string
		in   Expr
		want string
	}{
		{"integer", Int(5), "5"},
		{"symbol", Sym("x"), "x"},
		{"fraction", Rat(1, 2), "1\n─\n2"},
		{"power", Pow(Sym("x"), Int(2)), " 2\nx "},
		{"sum", Add(Sym("x"), Int(1)), "x + 1"},
		{"poly", Add(Pow(Sym("x"), Int(2)), Int(1)), " 2    \nx  + 1"},
		{"sqrt", Sqrt(Sym("x")), "  __\n/  x"},
	}
	for _, c := range cases {
		if got := Pretty(c.in); got != c.want {
			t.Errorf("%s: Pretty =\n%q\nwant\n%q", c.name, got, c.want)
		}
	}
}

// TestMathML exercises the presentation MathML renderer.
func TestMathML(t *testing.T) {
	const open = "<math xmlns=\"http://www.w3.org/1998/Math/MathML\">"
	const close = "</math>"
	cases := []struct {
		name string
		in   Expr
		want string
	}{
		{"symbol", Sym("x"), "<mi>x</mi>"},
		{"integer", Int(2), "<mn>2</mn>"},
		{"rational", Rat(1, 2), "<mfrac><mn>1</mn><mn>2</mn></mfrac>"},
		{"power", Pow(Sym("x"), Int(2)), "<msup><mi>x</mi><mn>2</mn></msup>"},
		{"sqrt", Sqrt(Sym("x")), "<msqrt><mi>x</mi></msqrt>"},
		{"sum", Add(Sym("x"), Int(1)), "<mrow><mi>x</mi><mo>+</mo><mn>1</mn></mrow>"},
		{"reciprocal", Pow(Sym("x"), Int(-1)), "<mfrac><mn>1</mn><mrow><mi>x</mi></mrow></mfrac>"},
	}
	for _, c := range cases {
		want := open + c.want + close
		if got := MathML(c.in); got != want {
			t.Errorf("%s: MathML = %q, want %q", c.name, got, want)
		}
	}
}
